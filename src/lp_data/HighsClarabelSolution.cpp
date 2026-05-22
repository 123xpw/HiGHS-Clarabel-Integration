// src/lp_data/HighsClarabelSolution.cpp
//
// Implementation of HighsClarabelSolution (M3).
// See HighsClarabelSolution.h for dual sign derivation and KKT cross-check.

#include "HighsClarabelSolution.h"

#include <cassert>

using clarabel::SolverStatus;

// ── mapStatus ────────────────────────────────────────────────────────────────

HighsModelStatus HighsClarabelSolution::mapStatus(SolverStatus cs) {
  switch (cs) {
    case SolverStatus::Solved:
    case SolverStatus::AlmostSolved:
      return HighsModelStatus::kOptimal;

    case SolverStatus::PrimalInfeasible:
    case SolverStatus::AlmostPrimalInfeasible:
      return HighsModelStatus::kInfeasible;

    case SolverStatus::DualInfeasible:
    case SolverStatus::AlmostDualInfeasible:
      // For LP, dual infeasibility of the Clarabel problem implies primal
      // unboundedness of the original HiGHS problem.
      return HighsModelStatus::kUnbounded;

    case SolverStatus::MaxIterations:
      return HighsModelStatus::kIterationLimit;

    case SolverStatus::MaxTime:
      return HighsModelStatus::kTimeLimit;

    case SolverStatus::NumericalError:
    case SolverStatus::InsufficientProgress:
    case SolverStatus::CallbackTerminated:
      return HighsModelStatus::kSolveError;

    case SolverStatus::Unsolved:
    default:
      return HighsModelStatus::kSolveError;
  }
}

// ── isFeasible ───────────────────────────────────────────────────────────────
//
// Returns true when Clarabel produced a usable primal + dual solution.
// AlmostSolved is included: the solution may not meet full tolerances but
// is still a valid warm-start and useful dual estimate.

bool HighsClarabelSolution::isFeasible(SolverStatus cs) {
  return cs == SolverStatus::Solved || cs == SolverStatus::AlmostSolved;
}

// ── recoverObjective ─────────────────────────────────────────────────────────
//
// For kMinimize: q = c, so HiGHS obj = Clarabel obj_val + offset.
// For kMaximize: q = −c, so Clarabel obj_val = −c'x.
//     HiGHS max obj = c'x = −obj_val → stored as −obj_val + offset.
// Compactly: highs_obj = es * obj_val + offset, where es = (kMin ? +1 : −1).

double HighsClarabelSolution::recoverObjective(double clarabel_obj,
                                               ObjSense sense,
                                               double   offset) {
  const double es = (sense == ObjSense::kMinimize) ? 1.0 : -1.0;
  return es * clarabel_obj + offset;
}

// ── write ────────────────────────────────────────────────────────────────────

HighsClarabelSolution::Result
HighsClarabelSolution::write(
    const clarabel::DefaultSolution<double>&     csol,
    const HighsClarabelInterface::ClarabelInput& input,
    const HighsLp&                               lp) {

  Result res;
  res.model_status    = mapStatus(csol.status);
  res.primal_feasible = isFeasible(csol.status);
  res.dual_feasible   = isFeasible(csol.status);

  const int nc = input.num_highs_cols;
  const int nr = input.num_highs_rows;

  res.solution.col_value.assign(nc, 0.0);
  res.solution.row_value.assign(nr, 0.0);
  res.solution.col_dual.assign(nc, 0.0);
  res.solution.row_dual.assign(nr, 0.0);

  if (!res.primal_feasible) {
    res.solution.value_valid = false;
    res.solution.dual_valid  = false;
    return res;
  }

  res.solution.value_valid = true;
  res.solution.dual_valid  = true;

  // ── Objective ──────────────────────────────────────────────────────────────
  res.objective_value =
      recoverObjective(csol.obj_val, input.sense, input.obj_offset);

  // ── Primal: col_value = x ─────────────────────────────────────────────────
  for (int j = 0; j < nc; ++j)
    res.solution.col_value[j] = csol.x[j];

  // ── Row activities: row_value = A * x (CSC traversal) ────────────────────
  const auto& mat = lp.a_matrix_;
  for (int j = 0; j < nc; ++j) {
    for (int k = static_cast<int>(mat.start_[j]);
         k < static_cast<int>(mat.start_[j + 1]); ++k) {
      res.solution.row_value[static_cast<int>(mat.index_[k])] +=
          mat.value_[k] * csol.x[j];
    }
  }

  // ── Duals ─────────────────────────────────────────────────────────────────
  //
  // effective_sign (es):
  //   kMinimize → +1  (shadow prices follow the natural sign)
  //   kMaximize → −1  (M2 negated q to convert max→min; duals flip back)
  //
  // Derivation (see header):  row_dual / col_dual formula per RowKind.
  const double es =
      (input.sense == ObjSense::kMinimize) ? 1.0 : -1.0;

  using RowKind = HighsClarabelInterface::RowKind;
  const int total_rows = static_cast<int>(input.row_map.size());
  assert(static_cast<int>(csol.z.size()) >= total_rows);

  for (int k = 0; k < total_rows; ++k) {
    const auto& entry = input.row_map[k];
    const double zk   = csol.z[k];

    switch (entry.kind) {
      case RowKind::kRowEq:
        res.solution.row_dual[entry.highs_idx]  = -zk * es;
        break;
      case RowKind::kRowUB:
        res.solution.row_dual[entry.highs_idx] += -zk * es;
        break;
      case RowKind::kRowLB:
        res.solution.row_dual[entry.highs_idx] += +zk * es;
        break;
      case RowKind::kColUB:
        res.solution.col_dual[entry.highs_idx] += -zk * es;
        break;
      case RowKind::kColLB:
        res.solution.col_dual[entry.highs_idx] += +zk * es;
        break;
    }
  }

  return res;
}
