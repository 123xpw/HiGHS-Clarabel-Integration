// src/lp_data/HighsClarabelInterface.cpp
//
// Implementation of HighsClarabelInterface (M2).
// See HighsClarabelInterface.h for the problem correspondence, row expansion
// rules, and dual sign conventions.

#include "HighsClarabelInterface.h"

#include <cassert>
#include <vector>

// ── buildFromLp ──────────────────────────────────────────────────────────────

HighsClarabelInterface::ClarabelInput
HighsClarabelInterface::buildFromLp(const HighsLp&      lp,
                                    const HighsOptions& opts) {
  const_cast<HighsLp&>(lp).ensureColwise();
  const int n = static_cast<int>(lp.num_col_);

  ClarabelInput in;
  in.P.resize(n, n);  // zero Hessian for LP
  in.q              = buildCostVector(lp);
  buildConstraints(lp, in.A, in.b, in.cones, in.row_map);
  in.settings       = makeSettings(opts);
  in.num_highs_cols = n;
  in.num_highs_rows = static_cast<int>(lp.num_row_);
  in.sense          = lp.sense_;
  in.obj_offset     = lp.offset_;
  return in;
}

// ── buildFromModel ───────────────────────────────────────────────────────────

HighsClarabelInterface::ClarabelInput
HighsClarabelInterface::buildFromModel(const HighsModel&   model,
                                        const HighsOptions& opts) {
  const HighsLp& lp = model.lp_;
  const_cast<HighsLp&>(lp).ensureColwise();
  const int n = static_cast<int>(lp.num_col_);

  ClarabelInput in;
  if (model.isQp()) {
    fillEigenHessian(model.hessian_, in.P);
  } else {
    in.P.resize(n, n);
  }
  in.q              = buildCostVector(lp);
  buildConstraints(lp, in.A, in.b, in.cones, in.row_map);
  in.settings       = makeSettings(opts);
  in.num_highs_cols = n;
  in.num_highs_rows = static_cast<int>(lp.num_row_);
  in.sense          = lp.sense_;
  in.obj_offset     = lp.offset_;
  return in;
}

// ── fillEigenHessian ─────────────────────────────────────────────────────────
//
// HighsHessian::kTriangular stores the upper triangle: column j contains
// row indices i ≤ j.  That is exactly what Clarabel expects (upper-tri CSC).
// kSquare is the full symmetric matrix; we keep only the upper triangle.

void HighsClarabelInterface::fillEigenHessian(
    const HighsHessian&          H,
    Eigen::SparseMatrix<double>& P) {
  const int n = static_cast<int>(H.dim_);
  P.resize(n, n);
  if (H.start_.empty() || n == 0) return;

  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(H.index_.size());

  if (H.format_ == HessianFormat::kTriangular) {
    for (int j = 0; j < n; ++j) {
      for (int k = static_cast<int>(H.start_[j]);
           k < static_cast<int>(H.start_[j + 1]); ++k) {
        trips.emplace_back(static_cast<int>(H.index_[k]), j, H.value_[k]);
      }
    }
  } else {
    // kSquare: full symmetric — keep i ≤ j (upper triangle)
    for (int j = 0; j < n; ++j) {
      for (int k = static_cast<int>(H.start_[j]);
           k < static_cast<int>(H.start_[j + 1]); ++k) {
        const int i = static_cast<int>(H.index_[k]);
        if (i <= j) trips.emplace_back(i, j, H.value_[k]);
      }
    }
  }

  P.setFromTriplets(trips.begin(), trips.end());
  P.makeCompressed();
}

// ── buildCostVector ──────────────────────────────────────────────────────────

Eigen::VectorXd
HighsClarabelInterface::buildCostVector(const HighsLp& lp) {
  const int    n    = static_cast<int>(lp.num_col_);
  const double sign = (lp.sense_ == ObjSense::kMaximize) ? -1.0 : 1.0;
  Eigen::VectorXd q(n);
  for (int j = 0; j < n; ++j) q[j] = sign * lp.col_cost_[j];
  return q;
}

// ── buildConstraints ─────────────────────────────────────────────────────────
//
// Builds Ã, b̃, cones, row_map per the row expansion rules in the header.
//
// Clarabel row ordering:
//   [0,        n_eq)         — equality LP rows     (ZeroConeT(1))
//   [n_eq,     n_lp_rows)    — inequality LP rows   (NonnegativeConeT(1))
//   [n_lp_rows, total_rows)  — column bound rows    (NonnegativeConeT(1))
//
// The "equalities first" ordering keeps the ZeroCone block contiguous while
// preserving a 1-to-1 correspondence between cones[k] and Clarabel row k.

void HighsClarabelInterface::buildConstraints(
    const HighsLp&                                   lp,
    Eigen::SparseMatrix<double>&                     A,
    Eigen::VectorXd&                                 b,
    std::vector<clarabel::SupportedConeT<double>>&   cones,
    std::vector<RowEntry>&                           row_map) {

  const int num_row = static_cast<int>(lp.num_row_);
  const int num_col = static_cast<int>(lp.num_col_);
  const auto& rl = lp.row_lower_;
  const auto& ru = lp.row_upper_;
  const auto& cl = lp.col_lower_;
  const auto& cu = lp.col_upper_;

  // ── Step 1: pre-assign Clarabel row indices for every HiGHS LP row.
  //
  // For each HiGHS row we record the list of Clarabel rows it maps to
  // (0, 1, or 2 entries).  We assign equalities first (indices 0..n_eq-1)
  // then inequalities, so two cursors advance independently.

  struct ClarabelRef { int idx; double sign; RowKind kind; };
  std::vector<std::vector<ClarabelRef>> row_to_clarabel(num_row);

  // Count equality rows to determine where the inequality block starts.
  int n_eq = 0;
  for (int i = 0; i < num_row; ++i) {
    bool lo_fin = !isInf(rl[i]);
    bool hi_fin = !isInf(ru[i]);
    if (lo_fin && hi_fin && rl[i] == ru[i]) ++n_eq;
  }

  int eq_cursor   = 0;      // fills Clarabel rows [0, n_eq)
  int ineq_cursor = n_eq;   // fills Clarabel rows [n_eq, ...)

  for (int i = 0; i < num_row; ++i) {
    bool lo_fin = !isInf(rl[i]);
    bool hi_fin = !isInf(ru[i]);

    if (lo_fin && hi_fin && rl[i] == ru[i]) {
      row_to_clarabel[i].push_back({eq_cursor++, +1.0, RowKind::kRowEq});
    } else if (!lo_fin && hi_fin) {
      row_to_clarabel[i].push_back({ineq_cursor++, +1.0, RowKind::kRowUB});
    } else if (lo_fin && !hi_fin) {
      row_to_clarabel[i].push_back({ineq_cursor++, -1.0, RowKind::kRowLB});
    } else if (lo_fin && hi_fin) {
      // double-sided: UB row (sign +1) then LB row (sign -1)
      row_to_clarabel[i].push_back({ineq_cursor++, +1.0, RowKind::kRowUB});
      row_to_clarabel[i].push_back({ineq_cursor++, -1.0, RowKind::kRowLB});
    }
    // free row (both infinite): no Clarabel rows emitted
  }

  const int n_lp_rows = ineq_cursor;

  // ── Step 2: count column-bound rows
  int n_col_rows = 0;
  for (int j = 0; j < num_col; ++j) {
    if (!isInf(cu[j])) ++n_col_rows;
    if (!isInf(cl[j])) ++n_col_rows;
  }

  const int total_rows = n_lp_rows + n_col_rows;

  // ── Step 3: allocate output containers
  b.resize(total_rows);
  cones.reserve(total_rows);
  row_map.reserve(total_rows);

  // ── Step 4: fill b, cones, row_map for LP rows
  //
  // Two sweeps keep equalities contiguous in the cones list.
  // Sweep A — equality rows
  for (int i = 0; i < num_row; ++i) {
    for (const auto& ref : row_to_clarabel[i]) {
      if (ref.kind != RowKind::kRowEq) continue;
      b[ref.idx] = rl[i];  // equality: sign = +1, value = lo == hi
      cones.push_back(clarabel::ZeroConeT<double>(1));
      row_map.push_back({ref.kind, i});
    }
  }
  // Sweep B — inequality rows
  for (int i = 0; i < num_row; ++i) {
    for (const auto& ref : row_to_clarabel[i]) {
      if (ref.kind == RowKind::kRowEq) continue;
      if (ref.kind == RowKind::kRowUB)
        b[ref.idx] = ru[i];
      else
        b[ref.idx] = -rl[i];  // kRowLB: row was negated
      cones.push_back(clarabel::NonnegativeConeT<double>(1));
      row_map.push_back({ref.kind, i});
    }
  }

  // ── Step 5: collect triplets for the LP part of Ã
  const auto& mat = lp.a_matrix_;
  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(mat.value_.size() * 2 + n_col_rows);

  for (int j = 0; j < num_col; ++j) {
    for (int k = static_cast<int>(mat.start_[j]);
         k < static_cast<int>(mat.start_[j + 1]); ++k) {
      const int    hi_row = static_cast<int>(mat.index_[k]);
      const double val    = mat.value_[k];
      for (const auto& ref : row_to_clarabel[hi_row]) {
        trips.emplace_back(ref.idx, j, ref.sign * val);
      }
    }
  }

  // ── Step 6: column-bound rows
  int col_row = n_lp_rows;
  for (int j = 0; j < num_col; ++j) {
    if (!isInf(cu[j])) {
      // kColUB:  x[j] + s = cu[j],  s ≥ 0  →  x[j] ≤ cu[j]
      trips.emplace_back(col_row, j, +1.0);
      b[col_row] = cu[j];
      cones.push_back(clarabel::NonnegativeConeT<double>(1));
      row_map.push_back({RowKind::kColUB, j});
      ++col_row;
    }
    if (!isInf(cl[j])) {
      // kColLB: -x[j] + s = -cl[j], s ≥ 0  →  x[j] ≥ cl[j]
      trips.emplace_back(col_row, j, -1.0);
      b[col_row] = -cl[j];
      cones.push_back(clarabel::NonnegativeConeT<double>(1));
      row_map.push_back({RowKind::kColLB, j});
      ++col_row;
    }
  }

  assert(col_row == total_rows);
  assert(static_cast<int>(cones.size())   == total_rows);
  assert(static_cast<int>(row_map.size()) == total_rows);

  // ── Step 7: assemble Ã
  A.resize(total_rows, num_col);
  A.setFromTriplets(trips.begin(), trips.end());
  A.makeCompressed();
}

// ── makeSettings ─────────────────────────────────────────────────────────────
//
// Maps HiGHS solver options onto Clarabel DefaultSettings.
// Only non-zero / valid values override the Clarabel defaults.

clarabel::DefaultSettings<double>
HighsClarabelInterface::makeSettings(const HighsOptions& opts) {
  auto s = clarabel::DefaultSettings<double>::default_settings();

  // Time limit: HiGHS uses kHighsInf (== +inf) or 0 for "no limit"
  if (opts.time_limit > 0.0 && !std::isinf(opts.time_limit))
    s.time_limit = opts.time_limit;

  // Iteration limit (0 means "not set" in HiGHS)
  if (opts.ipm_iteration_limit > 0)
    s.max_iter = static_cast<uint32_t>(opts.ipm_iteration_limit);

  // Tolerances: HiGHS optimality tolerance → Clarabel primal+dual gap tolerances
  if (opts.ipm_optimality_tolerance > 0.0) {
    s.tol_gap_abs = opts.ipm_optimality_tolerance;
    s.tol_gap_rel = opts.ipm_optimality_tolerance;
  }
  if (opts.primal_feasibility_tolerance > 0.0)
    s.tol_feas = opts.primal_feasibility_tolerance;

  s.verbose = opts.output_flag;
  return s;
}
