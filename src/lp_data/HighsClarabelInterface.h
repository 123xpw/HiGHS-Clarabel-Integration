// src/lp_data/HighsClarabelInterface.h
//
// M2 — LP/QP data mapping layer.
// Converts a HiGHS HighsLp / HighsModel to the input format required by
// Clarabel's DefaultSolver<double>.
//
// Problem correspondence
// ─────────────────────
//   HiGHS          :  min/max  c'x         s.t.  lo_r ≤ Ax ≤ hi_r,  lo_c ≤ x ≤ hi_c
//   Clarabel       :  min      ½x'Px + q'x  s.t.  Ãx + s = b̃,  s ∈ K
//
// Row expansion rules (spec §5.M2)
// ─────────────────────────────────
//   HiGHS row type   Condition                  Clarabel rows         Cone per row
//   ─────────────    ─────────────────────────  ────────────────────  ───────────────
//   Equality         lo == hi                    Ã[k,:] =  A[i,:],
//                                                b̃[k]   =  lo[i]      ZeroConeT(1)
//   Upper bound      isinf(lo) && !isinf(hi)     Ã[k,:] =  A[i,:],
//                                                b̃[k]   =  hi[i]      NonnegativeConeT(1)
//   Lower bound      !isinf(lo) && isinf(hi)     Ã[k,:] = -A[i,:],
//                                                b̃[k]   = -lo[i]      NonnegativeConeT(1)
//   Double-sided     !isinf(lo) && !isinf(hi)    (UB row) + (LB row)  Two NonnegativeConeT(1)
//                      && lo != hi
//   Free             isinf(lo) && isinf(hi)      skipped               —
//
// Column bounds treated identically, with A = ±e_j (standard basis vector).
// Equalities before inequalities within each group to keep ZeroCone blocks
// contiguous, but individual ZeroConeT(1)/NonnegativeConeT(1) are used so the
// row_map index always lines up with the cone list index 1-to-1.
//
// Dual sign convention (for M3)
// ──────────────────────────────
//   Clarabel solves  min ½x'Px + q'x  and returns z (dual of Ãx + s = b̃).
//   KKT:  Px + q + Ã'z = 0,  z ∈ K*.
//   Envelope theorem:  ∂obj*/∂b̃[k] = −z[k]  (note the minus sign).
//
//   HiGHS row_dual[i] = shadow price = ∂obj*/∂rhs_i.
//   Let  es = (sense==kMinimize ? +1 : −1)  (kMaximize negates q, so duals flip).
//
//   For kMinimize (es=+1):
//     kRowEq   → row_dual[i]  = −z[k]           ∂b̃/∂rhs=+1 → −z[k]·(+1)
//     kRowUB   → row_dual[i] += −z[k]            ∂b̃/∂hi  =+1 → −z[k]·(+1)
//     kRowLB   → row_dual[i] += +z[k]            ∂b̃/∂lo  =−1 → −z[k]·(−1)=+z[k]
//     kColUB   → col_dual[j] += −z[k]            same analysis as kRowUB
//     kColLB   → col_dual[j] += +z[k]            same analysis as kRowLB
//
//   For kMaximize (es=−1): every sign above is negated.
//   Compact form used in M3:
//     kRowEq  : row_dual[i]  = −z[k] * es
//     kRowUB  : row_dual[i] += −z[k] * es
//     kRowLB  : row_dual[i] += +z[k] * es
//     kColUB  : col_dual[j] += −z[k] * es
//     kColLB  : col_dual[j] += +z[k] * es
//
//   KKT residual cross-check: A'·row_dual + col_dual = c  (for kMinimize)
//   Dual feasibility (kMinimize): row/col at LB → dual ≥ 0; at UB → dual ≤ 0.
//
// Lifetime note
// ─────────────
//   DefaultSolver copies all problem data into Rust memory at construction.
//   ClarabelInput therefore does NOT need to outlive the solver.

#pragma once

#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include <clarabel.hpp>

#include "lp_data/HConst.h"    // ObjSense, kHighsInf (= inf)
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "model/HighsHessian.h"
#include "model/HighsModel.h"

// ============================================================
// HighsClarabelInterface
// ============================================================
class HighsClarabelInterface {
 public:
  // ----------------------------------------------------------
  // RowKind — identifies the origin of each Clarabel constraint
  // row so that M3 can recover HiGHS row/column duals.
  // ----------------------------------------------------------
  enum class RowKind : uint8_t {
    kRowEq,   ///< ZeroCone:    Ã[k,:] =  A[i,:],  b̃[k] =  lo[i]   (equality)
    kRowUB,   ///< NonnegCone:  Ã[k,:] =  A[i,:],  b̃[k] =  hi[i]   (≤ bound)
    kRowLB,   ///< NonnegCone:  Ã[k,:] = -A[i,:],  b̃[k] = -lo[i]   (≥ bound, negated)
    kColUB,   ///< NonnegCone:  x[j] + s = col_upper[j]
    kColLB,   ///< NonnegCone: -x[j] + s = -col_lower[j]
  };

  // One entry per Clarabel row, in the same order as the rows of Ã.
  struct RowEntry {
    RowKind kind;
    int     highs_idx;  ///< HighsLp row index (kRowXxx) or column index (kColXxx)
  };

  // ----------------------------------------------------------
  // ClarabelInput — all data needed to construct DefaultSolver.
  // Can be discarded once DefaultSolver is constructed because
  // Clarabel copies the matrices into Rust.
  // ----------------------------------------------------------
  struct ClarabelInput {
    // Problem matrices (Clarabel convention: P upper-triangular CSC)
    Eigen::SparseMatrix<double> P;
    Eigen::VectorXd             q;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd             b;

    // Cone list: one entry per constraint row, parallel to rows of A / b.
    std::vector<clarabel::SupportedConeT<double>> cones;

    // Solver settings derived from HighsOptions.
    clarabel::DefaultSettings<double> settings;

    // Row mapping: row_map[k] describes Clarabel row k.
    // Length == A.rows() == b.size() == cones.size().
    std::vector<RowEntry> row_map;

    // Original problem dimensions (needed by M3 to size output vectors).
    int      num_highs_cols;
    int      num_highs_rows;

    // Objective metadata needed by M3 to recover the HiGHS objective value.
    //   highs_obj = sense_sign * clarabel_obj_val + obj_offset
    // where sense_sign = (sense == kMinimize) ? +1 : -1.
    ObjSense sense;
    double   obj_offset;
  };

  // ----------------------------------------------------------
  // Public entry points
  // ----------------------------------------------------------

  /// Convert a pure LP (no hessian) to ClarabelInput.
  /// lp.a_matrix_ is ensured colwise internally.
  static ClarabelInput buildFromLp(const HighsLp&      lp,
                                   const HighsOptions& opts);

  /// Convert an LP or QP (model.isQp() ? fills P : P = 0).
  static ClarabelInput buildFromModel(const HighsModel&   model,
                                      const HighsOptions& opts);

 private:
  // ----------------------------------------------------------
  // Internal helpers
  // ----------------------------------------------------------

  /// Build the Eigen P matrix from HighsHessian.
  /// HighsHessian::kTriangular is stored as upper triangle (column j has
  /// row indices i ≤ j), which matches Clarabel's expected upper-triangular
  /// CSC format exactly.  kSquare is reduced to upper triangle here.
  static void fillEigenHessian(const HighsHessian&          H,
                                Eigen::SparseMatrix<double>& P);

  /// Build Ã, b̃, cones, and row_map from HighsLp.
  ///
  /// Pass order:
  ///   1. LP rows (equality, then inequality per row in order):
  ///      equality rows → kRowEq / ZeroConeT(1)
  ///      inequality rows → kRowUB and/or kRowLB / NonnegativeConeT(1)
  ///      free rows → skipped
  ///   2. Column bounds (for each column j, finite bounds only):
  ///      finite col_upper → kColUB / NonnegativeConeT(1)
  ///      finite col_lower → kColLB / NonnegativeConeT(1)
  ///
  /// lp.a_matrix_ must be colwise (CSC) before calling.
  static void buildConstraints(
      const HighsLp&                                   lp,
      Eigen::SparseMatrix<double>&                     A,
      Eigen::VectorXd&                                 b,
      std::vector<clarabel::SupportedConeT<double>>&   cones,
      std::vector<RowEntry>&                           row_map);

  /// Populate q from lp.col_cost_, adjusting for sense and offset.
  ///   • Maximisation: negate costs so Clarabel always minimises.
  ///   • offset_ is stored separately in ClarabelInput::obj_offset.
  static Eigen::VectorXd buildCostVector(const HighsLp& lp);

  /// Build Clarabel settings from HiGHS solver options.
  static clarabel::DefaultSettings<double> makeSettings(
      const HighsOptions& opts);

  // ----------------------------------------------------------
  // Infinity helpers (prefer std::isinf over == kHighsInf)
  // ----------------------------------------------------------
  static bool isInf(double v) { return std::isinf(v); }
};
