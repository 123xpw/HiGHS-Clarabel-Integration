// src/lp_data/HighsClarabelSolution.h
//
// M3 — solution write-back layer.
// Reads a Clarabel DefaultSolution<double> and recovers HighsSolution
// (col_value, row_value, col_dual, row_dual) plus HighsModelStatus.
//
// Dual sign convention
// ────────────────────
// Clarabel solves  min ½x'Px + q'x  s.t.  Ãx + s = b̃,  s ∈ K.
// The Lagrangian is  L = q'x + z'(Ãx + s − b̃), so by the envelope theorem:
//     ∂obj*/∂b̃[k]  =  −z[k]
//
// row_dual[i] = ∂obj*/∂rhs_i (shadow price of constraint i, in HiGHS sense).
// Mapping from Clarabel k to HiGHS row/col dual:
//
//   let es = (sense == kMinimize) ? +1.0 : −1.0
//   (sense == kMaximize negates q; the objective shadow price also negates.)
//
//   row_map[k].kind   Formula                         Rationale
//   ───────────────   ──────────────────────────────  ──────────────────────────────────────
//   kRowEq            row_dual[i]   = −z[k] * es      b̃ = lo = hi,  ∂b̃/∂rhs = +1
//   kRowUB            row_dual[i]  += −z[k] * es      b̃ = hi,       ∂b̃/∂hi  = +1
//   kRowLB            row_dual[i]  += +z[k] * es      b̃ = −lo,      ∂b̃/∂lo  = −1  →  sign flips
//   kColUB            col_dual[j]  += −z[k] * es      same analysis as kRowUB
//   kColLB            col_dual[j]  += +z[k] * es      same analysis as kRowLB
//
// Note: M2's header comments document the kMaximize formulas under the
// kMinimize label (UB/LB signs are swapped there).  The table above is correct.
//
// For double-sided rows (kRowUB + kRowLB pair for the same HiGHS row i):
//     row_dual[i] = (−z[k_ub] + z[k_lb]) * es
// which accumulates correctly with the += formulas above.
//
// Objective recovery
// ──────────────────
// For kMinimize: q = c, so  HiGHS obj = Clarabel obj_val + offset.
// For kMaximize: q = −c,  so Clarabel obj_val = −c'x.
//     HiGHS max obj = c'x = −obj_val,  stored value = −obj_val + offset.
// In both cases:  highs_obj = es * obj_val + offset   (es as above).
//
// KKT residual cross-check (matches HiGHS internal check at HighsSolution.cpp:393)
// ────────────────────────────────────────────────────────────────────────────────
// The LP column KKT condition is:
//     A' * row_dual  +  col_dual  =  c          (for kMinimize)
// Equivalently:  |(A'y − c) + col_dual| ≈ 0    (HiGHS residual form)
// This identity holds because col_dual = −(col-bound part of Ã'z)
// and the LP-row part of Ã'z equals −A'*row_dual.

#pragma once

#include <cmath>
#include <limits>
#include <vector>

#include <clarabel.hpp>

#include "lp_data/HConst.h"    // HighsModelStatus, ObjSense, kHighsInf
#include "lp_data/HStruct.h"   // HighsSolution
#include "lp_data/HighsClarabelInterface.h"
#include "lp_data/HighsLp.h"

// ============================================================
// HighsClarabelSolution
// ============================================================
class HighsClarabelSolution {
 public:
  // Return bundle from write().
  struct Result {
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsSolution    solution;
    double           objective_value = std::numeric_limits<double>::quiet_NaN();
    bool             primal_feasible = false;
    bool             dual_feasible   = false;
  };

  // Convert a Clarabel solution to a HiGHS Result.
  // lp.a_matrix_ must be colwise (CSC).
  static Result write(const clarabel::DefaultSolution<double>&        csol,
                      const HighsClarabelInterface::ClarabelInput&    input,
                      const HighsLp&                                  lp);

 private:
  // Map Clarabel SolverStatus to HighsModelStatus.
  static HighsModelStatus mapStatus(clarabel::SolverStatus cs);

  // True when the Clarabel status carries a usable primal + dual solution.
  static bool isFeasible(clarabel::SolverStatus cs);

  // Recover the HiGHS-sense objective from the Clarabel obj_val.
  static double recoverObjective(double clarabel_obj, ObjSense sense,
                                  double offset);
};
