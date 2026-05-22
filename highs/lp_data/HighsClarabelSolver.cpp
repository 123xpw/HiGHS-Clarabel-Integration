/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsClarabelSolver.cpp
 * @brief M4 — Clarabel LP solver integration (solveLpClarabel implementation).
 *
 * Compiled only when HIGHS_USE_CLARABEL is defined.
 *
 * M2 (LP→Clarabel data mapping) and M3 (Clarabel solution→HiGHS write-back)
 * are implemented inline here, keeping this translation unit self-contained
 * and avoiding duplicate-symbol issues with the standalone src/lp_data library.
 *
 * Row expansion rules (M2):
 *   HiGHS row type   Condition                   Clarabel row     Cone
 *   Equality         lo == hi                    Ã =  A, b̃ = lo  ZeroConeT(1)
 *   Upper-bound      isinf(lo) && !isinf(hi)     Ã =  A, b̃ = hi  NonnegConeT(1)
 *   Lower-bound      !isinf(lo) && isinf(hi)     Ã = -A, b̃ =-lo  NonnegConeT(1)
 *   Double-sided     both finite, lo != hi        UB row + LB row  two NonnegConeT(1)
 *   Free             both infinite               (skipped)
 * Column bounds add one extra row per finite bound.
 *
 * Dual sign convention (M3):
 *   es = (kMinimize ? +1 : -1)
 *   kRowEq:  row_dual[i]  = -z[k]*es
 *   kRowUB:  row_dual[i] += -z[k]*es
 *   kRowLB:  row_dual[i] += +z[k]*es
 *   kColUB:  col_dual[j] += -z[k]*es
 *   kColLB:  col_dual[j] += +z[k]*es
 */

#ifdef HIGHS_USE_CLARABEL

#include "lp_data/HighsClarabelSolver.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/Sparse>
#include <clarabel.hpp>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "model/HighsHessian.h"
#include "presolve/ICrashX.h"

// ── internal helpers ─────────────────────────────────────────────────────────

namespace {

// True when v represents +∞ in HiGHS (kHighsInf ≈ 1e30 or IEEE infinity).
static inline bool isHInf(double v) { return std::isinf(v) || v >= kHighsInf; }
static inline bool isNInf(double v) { return std::isinf(v) || v <= -kHighsInf; }

// Tags identifying the origin of each Clarabel constraint row.
enum class RowKind : uint8_t { kRowEq, kRowUB, kRowLB, kColUB, kColLB };
struct RowEntry { RowKind kind; int highs_idx; };

// ── buildClarabelInput ───────────────────────────────────────────────────────
//
// Converts a HiGHS LP to the Clarabel input format.
// On return: P (zero for LP), q, A, b, cones, row_map, settings are populated.
struct ClarabelInput {
  Eigen::SparseMatrix<double>                     P;
  Eigen::VectorXd                                 q;
  Eigen::SparseMatrix<double>                     A;
  Eigen::VectorXd                                 b;
  std::vector<clarabel::SupportedConeT<double>>   cones;
  clarabel::DefaultSettings<double>               settings;
  std::vector<RowEntry>                           row_map;
  std::vector<double>                             col_shift;  // x = x' + shift
  int      num_highs_cols;
  int      num_highs_rows;
  ObjSense sense;
  double   obj_offset;
};

static ClarabelInput buildClarabelInput(const HighsLp&        lp,
                                         const HighsOptions&   opts,
                                         const HighsHessian*   hessian_ptr = nullptr) {
  const_cast<HighsLp&>(lp).ensureColwise();
  const int n       = static_cast<int>(lp.num_col_);
  const int num_row = static_cast<int>(lp.num_row_);
  const auto& rl = lp.row_lower_;
  const auto& ru = lp.row_upper_;
  const auto& cl = lp.col_lower_;
  const auto& cu = lp.col_upper_;

  ClarabelInput in;
  in.P.resize(n, n);  // zero by default; filled below for QP

  // ── Hessian (QP mode) ─────────────────────────────────────────────────────
  // HiGHS stores Q in upper-triangular CSC (kTriangular format).
  // Clarabel expects P upper-triangular.  Direct copy; negate for maximization.
  if (hessian_ptr != nullptr && hessian_ptr->dim_ > 0 &&
      hessian_ptr->numNz() > 0) {
    const double sign = (lp.sense_ == ObjSense::kMaximize) ? -1.0 : 1.0;
    const int dim = static_cast<int>(hessian_ptr->dim_);
    std::vector<Eigen::Triplet<double>> p_trips;
    p_trips.reserve(static_cast<size_t>(hessian_ptr->numNz()));
    for (int j = 0; j < dim; ++j) {
      for (int k = static_cast<int>(hessian_ptr->start_[j]);
           k < static_cast<int>(hessian_ptr->start_[j + 1]); ++k) {
        const int i = static_cast<int>(hessian_ptr->index_[k]);
        // HiGHS Hessian data may be stored as either triangular orientation
        // depending on the source format. Clarabel expects the upper triangle.
        p_trips.emplace_back(std::min(i, j), std::max(i, j),
                             sign * hessian_ptr->value_[k]);
      }
    }
    in.P.setFromTriplets(p_trips.begin(), p_trips.end());
    in.P.makeCompressed();
  }
  in.num_highs_cols = n;
  in.num_highs_rows = num_row;
  in.sense          = lp.sense_;
  in.obj_offset     = lp.offset_;

  // ── Variable shift: x = x' + shift (LP only) ───────────────────────────
  //
  // Shifting variables so their effective lower bound is 0 ensures that all
  // kColLB Clarabel rows have b = 0 (NonnegCone with b < 0 triggers a known
  // Clarabel crash on some instances such as 25fv47).  Skipped for QP because
  // the quadratic term would also need a q and offset correction.
  const auto& mat = lp.a_matrix_;
  in.col_shift.assign(n, 0.0);
  std::vector<double> row_b_adj(num_row, 0.0);
  if (hessian_ptr == nullptr || hessian_ptr->numNz() == 0) {
    for (int j = 0; j < n; ++j)
      if (!isNInf(cl[j])) in.col_shift[j] = cl[j];

    // Adjust objective offset: c^T (x' + shift) + offset = c^T x' + (c^T shift + offset)
    for (int j = 0; j < n; ++j)
      in.obj_offset += lp.col_cost_[j] * in.col_shift[j];

    // Compute A*shift for RHS adjustments: new b_i = b_i - A[i,:]*shift
    for (int j = 0; j < n; ++j) {
      if (in.col_shift[j] == 0.0) continue;
      for (int k = static_cast<int>(mat.start_[j]);
           k < static_cast<int>(mat.start_[j + 1]); ++k)
        row_b_adj[static_cast<int>(mat.index_[k])] +=
            mat.value_[k] * in.col_shift[j];
    }
  }

  // ── cost vector q ────────────────────────────────────────────────────────
  {
    const double sign = (lp.sense_ == ObjSense::kMaximize) ? -1.0 : 1.0;
    in.q.resize(n);
    for (int j = 0; j < n; ++j) in.q[j] = sign * lp.col_cost_[j];
  }

  // ── Clarabel row index assignment ─────────────────────────────────────────
  //
  // Equalities first (ZeroCone), then inequalities (NonnegCone).
  // Each HiGHS row maps to 0, 1, or 2 Clarabel rows.
  struct CRef { int idx; double sign; RowKind kind; };
  std::vector<std::vector<CRef>> r2c(num_row);

  int n_eq = 0;
  for (int i = 0; i < num_row; ++i)
    if (!isNInf(rl[i]) && !isHInf(ru[i]) && rl[i] == ru[i]) ++n_eq;

  int eq_cur   = 0;
  int ineq_cur = n_eq;
  for (int i = 0; i < num_row; ++i) {
    bool lo_fin = !isNInf(rl[i]);
    bool hi_fin = !isHInf(ru[i]);
    if (lo_fin && hi_fin && rl[i] == ru[i]) {
      r2c[i].push_back({eq_cur++,   +1.0, RowKind::kRowEq});
    } else if (!lo_fin && hi_fin) {
      r2c[i].push_back({ineq_cur++, +1.0, RowKind::kRowUB});
    } else if (lo_fin && !hi_fin) {
      r2c[i].push_back({ineq_cur++, -1.0, RowKind::kRowLB});
    } else if (lo_fin && hi_fin) {
      r2c[i].push_back({ineq_cur++, +1.0, RowKind::kRowUB});
      r2c[i].push_back({ineq_cur++, -1.0, RowKind::kRowLB});
    }
    // free row: no Clarabel rows
  }
  const int n_lp_rows = ineq_cur;

  // Count column-bound rows.
  int n_col_rows = 0;
  for (int j = 0; j < n; ++j) {
    if (!isHInf(cu[j])) ++n_col_rows;
    if (!isNInf(cl[j])) ++n_col_rows;
  }
  const int total_rows = n_lp_rows + n_col_rows;

  in.b.resize(total_rows);
  in.cones.reserve(total_rows);
  in.row_map.reserve(total_rows);

  // ── b, cones, row_map: equality rows first ────────────────────────────────
  for (int i = 0; i < num_row; ++i)
    for (const auto& ref : r2c[i])
      if (ref.kind == RowKind::kRowEq) {
        in.b[ref.idx] = rl[i] - row_b_adj[i];
        in.cones.push_back(clarabel::ZeroConeT<double>(1));
        in.row_map.push_back({RowKind::kRowEq, i});
      }
  // inequality rows
  for (int i = 0; i < num_row; ++i)
    for (const auto& ref : r2c[i])
      if (ref.kind != RowKind::kRowEq) {
        if (ref.kind == RowKind::kRowUB)
          in.b[ref.idx] = ru[i] - row_b_adj[i];
        else  // kRowLB: -A*x' <= -lo  →  b = -lo + A*shift
          in.b[ref.idx] = -rl[i] + row_b_adj[i];
        in.cones.push_back(clarabel::NonnegativeConeT<double>(1));
        in.row_map.push_back({ref.kind, i});
      }

  // ── Triplets for the LP portion of Ã ─────────────────────────────────────
  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(mat.value_.size() * 2 + n_col_rows);

  for (int j = 0; j < n; ++j)
    for (int k = static_cast<int>(mat.start_[j]);
         k < static_cast<int>(mat.start_[j + 1]); ++k) {
      const int    hi = static_cast<int>(mat.index_[k]);
      const double v  = mat.value_[k];
      for (const auto& ref : r2c[hi])
        trips.emplace_back(ref.idx, j, ref.sign * v);
    }

  // ── Column-bound rows ─────────────────────────────────────────────────────
  int col_row = n_lp_rows;
  for (int j = 0; j < n; ++j) {
    if (!isHInf(cu[j])) {
      trips.emplace_back(col_row, j, +1.0);
      // After shift: x'[j] <= cu[j] - shift[j]
      in.b[col_row] = cu[j] - in.col_shift[j];
      in.cones.push_back(clarabel::NonnegativeConeT<double>(1));
      in.row_map.push_back({RowKind::kColUB, j});
      ++col_row;
    }
    if (!isNInf(cl[j])) {
      trips.emplace_back(col_row, j, -1.0);
      // After shift: x'[j] >= 0  (shift makes lb = 0)
      in.b[col_row] = 0.0;
      in.cones.push_back(clarabel::NonnegativeConeT<double>(1));
      in.row_map.push_back({RowKind::kColLB, j});
      ++col_row;
    }
  }
  assert(col_row == total_rows);
  assert(static_cast<int>(in.cones.size())   == total_rows);
  assert(static_cast<int>(in.row_map.size()) == total_rows);

  in.A.resize(total_rows, n);
  in.A.setFromTriplets(trips.begin(), trips.end());
  in.A.makeCompressed();

  // ── Clarabel settings ────────────────────────────────────────────────────
  auto& s = in.settings;
  s = clarabel::DefaultSettings<double>::default_settings();
  if (opts.time_limit > 0.0 && !std::isinf(opts.time_limit))
    s.time_limit = opts.time_limit;
  if (opts.ipm_iteration_limit > 0)
    s.max_iter = static_cast<uint32_t>(opts.ipm_iteration_limit);
  if (opts.ipm_optimality_tolerance > 0.0) {
    s.tol_gap_abs = opts.ipm_optimality_tolerance;
    s.tol_gap_rel = opts.ipm_optimality_tolerance;
  }
  if (opts.primal_feasibility_tolerance > 0.0)
    s.tol_feas = opts.primal_feasibility_tolerance;
  s.verbose = opts.output_flag;

  return in;
}

// ── mapStatus ────────────────────────────────────────────────────────────────

static HighsModelStatus mapClarabelStatus(clarabel::SolverStatus cs) {
  using S = clarabel::SolverStatus;
  switch (cs) {
    case S::Solved:
    case S::AlmostSolved:
      return HighsModelStatus::kOptimal;
    case S::PrimalInfeasible:
    case S::AlmostPrimalInfeasible:
      return HighsModelStatus::kInfeasible;
    case S::DualInfeasible:
    case S::AlmostDualInfeasible:
      return HighsModelStatus::kUnbounded;
    case S::MaxIterations:
      return HighsModelStatus::kIterationLimit;
    case S::MaxTime:
      return HighsModelStatus::kTimeLimit;
    case S::NumericalError:
    case S::InsufficientProgress:
    case S::CallbackTerminated:
    case S::Unsolved:
    default:
      return HighsModelStatus::kSolveError;
  }
}

static bool isClarabelFeasible(clarabel::SolverStatus cs) {
  return cs == clarabel::SolverStatus::Solved ||
         cs == clarabel::SolverStatus::AlmostSolved;
}

// ── writeClarabelSolution ────────────────────────────────────────────────────
//
// Fills highs_solution from Clarabel output.  Returns the HiGHS objective.

static double writeClarabelSolution(
    const clarabel::DefaultSolution<double>& csol,
    const ClarabelInput&                     in,
    const HighsLp&                           lp,
    HighsSolution&                           highs_solution) {

  const int nc = in.num_highs_cols;
  const int nr = in.num_highs_rows;
  highs_solution.col_value.assign(nc, 0.0);
  highs_solution.row_value.assign(nr, 0.0);
  highs_solution.col_dual.assign(nc, 0.0);
  highs_solution.row_dual.assign(nr, 0.0);
  highs_solution.value_valid = true;
  highs_solution.dual_valid  = true;

  // Objective: es * obj_val + offset
  const double es  = (in.sense == ObjSense::kMinimize) ? 1.0 : -1.0;
  const double obj = es * csol.obj_val + in.obj_offset;

  // Primal solution: unshift x' back to original variables
  for (int j = 0; j < nc; ++j)
    highs_solution.col_value[j] = csol.x[j] + in.col_shift[j];

  // Row activities: A*x (CSC traversal using unshifted col_value)
  const auto& mat = lp.a_matrix_;
  for (int j = 0; j < nc; ++j)
    for (int k = static_cast<int>(mat.start_[j]);
         k < static_cast<int>(mat.start_[j + 1]); ++k)
      highs_solution.row_value[static_cast<int>(mat.index_[k])] +=
          mat.value_[k] * highs_solution.col_value[j];

  // Duals
  const int total_rows = static_cast<int>(in.row_map.size());
  assert(static_cast<int>(csol.z.size()) >= total_rows);
  for (int k = 0; k < total_rows; ++k) {
    const auto& entry = in.row_map[k];
    const double zk   = csol.z[k];
    switch (entry.kind) {
      case RowKind::kRowEq:
        highs_solution.row_dual[entry.highs_idx]  = -zk * es;
        break;
      case RowKind::kRowUB:
        highs_solution.row_dual[entry.highs_idx] += -zk * es;
        break;
      case RowKind::kRowLB:
        highs_solution.row_dual[entry.highs_idx] += +zk * es;
        break;
      case RowKind::kColUB:
        highs_solution.col_dual[entry.highs_idx] += -zk * es;
        break;
      case RowKind::kColLB:
        highs_solution.col_dual[entry.highs_idx] += +zk * es;
        break;
    }
  }
  return obj;
}

}  // namespace

// ── solveLpClarabel ──────────────────────────────────────────────────────────

HighsStatus solveLpClarabel(HighsLpSolverObject& solver_object,
                             const HighsHessian*  hessian_ptr) {
  HighsLp&       lp      = solver_object.lp_;
  HighsOptions&  opts    = solver_object.options_;
  HighsSolution& sol     = solver_object.solution_;
  HighsBasis&    basis   = solver_object.basis_;
  HighsInfo&     info    = solver_object.highs_info_;

  // Reset outputs
  basis.valid       = false;
  sol.value_valid   = false;
  sol.dual_valid    = false;
  resetModelStatusAndHighsInfo(solver_object);

  const bool is_qp = (hessian_ptr != nullptr && hessian_ptr->dim_ > 0 &&
                       hessian_ptr->numNz() > 0);
  highsLogUser(opts.log_options, HighsLogType::kInfo,
               "Solving %s with Clarabel (interior-point)\n",
               is_qp ? "QP" : "LP");

  // ── M2: build Clarabel input ──────────────────────────────────────────────
  ClarabelInput in;
  try {
    in = buildClarabelInput(lp, opts, hessian_ptr);
  } catch (const std::exception& e) {
    highsLogDev(opts.log_options, HighsLogType::kError,
                "Exception in Clarabel input construction: %s\n", e.what());
    solver_object.model_status_ = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  if (in.A.rows() == 0) {
    // No constraints after mapping — degenerate; let the unconstrained solver handle it
    solver_object.model_status_ = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // ── Clarabel solve + M3 recovery ────────────────────────────────────────
  //
  // DefaultSolution has no default constructor, so we must initialise it
  // directly from solver.solution() inside the try block.  State is propagated
  // via the variables below.
  HighsModelStatus clarabel_model_status = HighsModelStatus::kSolveError;
  bool             feasible              = false;
  double           obj                   = 0.0;

  try {
    clarabel::DefaultSolver<double> cla_solver(in.P, in.q, in.A, in.b,
                                               in.cones, in.settings);
    cla_solver.solve();
    auto csol = cla_solver.solution();

    clarabel_model_status = mapClarabelStatus(csol.status);
    feasible              = isClarabelFeasible(csol.status);

    if (feasible) {
      obj = writeClarabelSolution(csol, in, lp, sol);

      // ── Primal solution quality gate ─────────────────────────────────────
      //
      // Clarabel may report Solved while the HiGHS-side residuals still
      // exceed tolerance (e.g. greenbea), or return numerically extreme
      // values that crash IPX crossover (e.g. 25fv47 heap corruption).
      //
      // Check each col_value / row_value for:
      //   (a) finiteness and magnitude  < 1e20  (crossover safety)
      //   (b) primal feasibility within 100 × primal_tolerance for finite bounds
      //
      // On failure: mark kUnknown so solveLp() routes to simplex cleanup,
      // and set feasible=false to skip crossover entirely.
      const double ptol    = 100.0 * opts.primal_feasibility_tolerance;
      const double max_val = 1e20;
      bool quality_ok = true;

      for (int j = 0; quality_ok && j < static_cast<int>(lp.num_col_); ++j) {
        const double v = sol.col_value[j];
        if (!std::isfinite(v) || std::fabs(v) > max_val) {
          quality_ok = false;
        } else {
          if (!isNInf(lp.col_lower_[j]) && v < lp.col_lower_[j] - ptol)
            quality_ok = false;
          if (!isHInf(lp.col_upper_[j]) && v > lp.col_upper_[j] + ptol)
            quality_ok = false;
        }
      }
      for (int i = 0; quality_ok && i < static_cast<int>(lp.num_row_); ++i) {
        const double v = sol.row_value[i];
        if (!std::isfinite(v) || std::fabs(v) > max_val) {
          quality_ok = false;
        } else {
          if (!isNInf(lp.row_lower_[i]) && v < lp.row_lower_[i] - ptol)
            quality_ok = false;
          if (!isHInf(lp.row_upper_[i]) && v > lp.row_upper_[i] + ptol)
            quality_ok = false;
        }
      }

      if (!quality_ok) {
        highsLogUser(opts.log_options, HighsLogType::kWarning,
                     "Clarabel: primal solution quality insufficient "
                     "(non-finite, magnitude > 1e20, or bound violation > tol x100) "
                     "— routing to simplex cleanup\n");
        clarabel_model_status = HighsModelStatus::kUnknown;
        feasible = false;  // skip crossover to avoid potential crash
      }
    }
  } catch (const std::exception& e) {
    highsLogDev(opts.log_options, HighsLogType::kError,
                "Exception in Clarabel solve: %s\n", e.what());
    solver_object.model_status_ = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  solver_object.model_status_ = clarabel_model_status;

  if (!feasible) {
    sol.value_valid = false;
    sol.dual_valid  = false;
    return HighsStatus::kOk;
  }

  info.objective_function_value = obj;

  // ── Crossover ─────────────────────────────────────────────────────────────
  //
  // Clarabel is an IPM solver — it produces no simplex basis.
  // Run IPX crossover if the user hasn't disabled it.
  // Only reached when quality_ok == true (guard above skips this on bad sol).
  const bool run_crossover =
      (opts.run_crossover != kHighsOffString);
  if (run_crossover) {
    HighsStatus co_status =
        callCrossover(opts, lp, basis, sol, solver_object.model_status_,
                      info, solver_object.callback_);
    if (co_status == HighsStatus::kError) {
      highsLogUser(opts.log_options, HighsLogType::kWarning,
                   "Clarabel: crossover failed; returning IPM solution without basis\n");
      basis.valid = false;
      // Solution is still valid; keep kOptimal
    }
    // co_status == kWarning is acceptable (imprecise duals; simplex cleans up)
  }

  return HighsStatus::kOk;
}

#endif  // HIGHS_USE_CLARABEL
