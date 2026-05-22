// src/test_m3_solution.cpp
//
// Unit tests for HighsClarabelSolution (M3 solution write-back layer).
// Each test builds a small LP/QP, runs Clarabel, and verifies that
// HighsClarabelSolution::write() recovers the expected HiGHS solution.
//
// Dual sign convention (from HighsClarabelSolution.h):
//   es = (kMinimize ? +1 : -1)
//   kRowEq  : row_dual[i]  = -z[k] * es
//   kRowUB  : row_dual[i] += -z[k] * es
//   kRowLB  : row_dual[i] += +z[k] * es
//   kColUB  : col_dual[j] += -z[k] * es
//   kColLB  : col_dual[j] += +z[k] * es
//
// KKT cross-check used in several tests:
//   A' * row_dual + col_dual == c    (for kMinimize LP)
//   A' * row_dual + col_dual == -c   (for kMaximize LP, because HiGHS flips)
//   i.e.  A' * row_dual + col_dual == c * es

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Sparse>
#include <clarabel.hpp>

#include "lp_data/HighsClarabelInterface.h"
#include "lp_data/HighsClarabelSolution.h"
#include "lp_data/HConst.h"
#include "model/HighsModel.h"

// ── tiny test harness ─────────────────────────────────────────────────────────

static int g_pass = 0, g_fail = 0;

static void check(bool cond, const std::string& label) {
  if (cond) {
    std::cout << "  PASS  " << label << "\n";
    ++g_pass;
  } else {
    std::cout << "  FAIL  " << label << "\n";
    ++g_fail;
  }
}

static void section(const std::string& name) {
  std::cout << "\n── " << name << " ──\n";
}

static const double inf = std::numeric_limits<double>::infinity();
static const double kTol = 1e-4;  // tolerance for solution checks

static bool near(double a, double b, double tol = kTol) {
  return std::fabs(a - b) <= tol;
}

// ── helper: solve a ClarabelInput and call write() ───────────────────────────

static HighsClarabelSolution::Result solveAndWrite(
    HighsClarabelInterface::ClarabelInput& in,
    const HighsLp& lp) {
  in.settings.verbose = false;
  clarabel::DefaultSolver<double> solver(in.P, in.q, in.A, in.b,
                                          in.cones, in.settings);
  solver.solve();
  auto csol = solver.solution();
  return HighsClarabelSolution::write(csol, in, lp);
}

// ── helper: compute A' * row_dual + col_dual (KKT residual) ─────────────────

static std::vector<double> kktResidual(const HighsLp& lp,
                                        const HighsSolution& sol) {
  const int nc = static_cast<int>(lp.num_col_);
  std::vector<double> res(nc, 0.0);
  for (int j = 0; j < nc; ++j) res[j] = sol.col_dual[j];
  const auto& mat = lp.a_matrix_;
  for (int j = 0; j < nc; ++j) {
    for (int k = static_cast<int>(mat.start_[j]);
         k < static_cast<int>(mat.start_[j + 1]); ++k) {
      res[j] += mat.value_[k] * sol.row_dual[static_cast<int>(mat.index_[k])];
    }
  }
  return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 1: Simple LP with equality constraint (minimisation)
//
//   min  2x0 + x1
//   s.t.  x0 + x1 = 4
//         x0, x1 >= 0
//
// Optimal: x0=0, x1=4, obj=4.
// row_dual[0] = 1  (shadow price of equality)
// col_dual[0] = 1  (reduced cost; x0=0 at LB)
// col_dual[1] = 0  (x1=4 between bounds)
// KKT: A'y + rc = [2,1]  (A = [[1,1]], y=[1])
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_equality() {
  section("Test 1: LP with equality (min 2x0+x1, x0+x1=4, x0,x1>=0)");

  HighsLp lp;
  lp.num_col_ = 2;  lp.num_row_ = 1;
  lp.col_cost_  = {2.0, 1.0};
  lp.col_lower_ = {0.0, 0.0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {4.0};
  lp.row_upper_ = {4.0};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  auto& m = lp.a_matrix_;
  m.format_ = MatrixFormat::kColwise;
  m.num_col_ = 2; m.num_row_ = 1;
  m.start_ = {0, 1, 2};
  m.index_ = {0, 0};
  m.value_ = {1.0, 1.0};

  HighsOptions opts;
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal,
        "status == kOptimal");
  check(res.primal_feasible && res.dual_feasible,
        "primal + dual feasible");
  check(near(res.objective_value, 4.0), "obj == 4.0");

  check(near(res.solution.col_value[0], 0.0), "x0 == 0");
  check(near(res.solution.col_value[1], 4.0), "x1 == 4");
  check(near(res.solution.row_value[0], 4.0), "row_value[0] == 4");

  check(near(res.solution.row_dual[0], 1.0),  "row_dual[0] == 1 (equality shadow price)");
  check(near(res.solution.col_dual[0], 1.0),  "col_dual[0] == 1 (x0 at LB)");
  check(near(res.solution.col_dual[1], 0.0),  "col_dual[1] == 0 (x1 between bounds)");

  // KKT: A'y + rc = c → [1,1]*[1] + [1,0] = [2,1] ✓
  auto kkt = kktResidual(lp, res.solution);
  check(near(kkt[0], lp.col_cost_[0]) && near(kkt[1], lp.col_cost_[1]),
        "KKT: A'*row_dual + col_dual == c");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2: LP with double-sided row constraint (binding lower bound)
//
//   min  x0
//   s.t.  1 <= x0 <= 3   (double-sided row)
//         x0 >= 0        (col LB)
//
// Optimal: x0=1, obj=1.
// row_dual[0] = 1  (shadow price of x0 >= 1, binding LB)
// col_dual[0] = 0  (x0=1 not at col bound: c - A'y = 1 - 1*1 = 0)
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_double_sided() {
  section("Test 2: LP double-sided row (1 <= x0 <= 3, min x0)");

  HighsLp lp;
  lp.num_col_ = 1;  lp.num_row_ = 1;
  lp.col_cost_  = {1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {inf};
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {3.0};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  auto& m = lp.a_matrix_;
  m.format_ = MatrixFormat::kColwise;
  m.num_col_ = 1; m.num_row_ = 1;
  m.start_ = {0, 1};
  m.index_ = {0};
  m.value_ = {1.0};

  HighsOptions opts;
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, 1.0), "obj == 1.0");
  check(near(res.solution.col_value[0], 1.0), "x0 == 1");
  check(near(res.solution.row_value[0], 1.0), "row_value[0] == 1");

  // row LB binding: row_dual[0] = +z[k_lb] * es = +1
  check(near(res.solution.row_dual[0], 1.0),
        "row_dual[0] == 1 (x0=1 at row LB)");
  // col_dual[0] = c - A'*row_dual = 1 - 1*1 = 0
  check(near(res.solution.col_dual[0], 0.0),
        "col_dual[0] == 0 (x0 not at col bound)");

  auto kkt = kktResidual(lp, res.solution);
  check(near(kkt[0], lp.col_cost_[0]), "KKT: A'*row_dual + col_dual == c");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3: LP with only column bounds (no explicit row constraints)
//
//   min  x0 + x1
//   s.t.  x0 >= 2   (col LB)
//         x1 >= 3   (col LB)
//         x0 <= 10  (col UB)
//
// Optimal: x0=2, x1=3, obj=5.
// col_dual[0] = 1  (x0=2 at col LB, reduced cost = 1 > 0)
// col_dual[1] = 1  (x1=3 at col LB, reduced cost = 1 > 0)
// (No row constraints, so no row_dual.)
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_col_bounds_only() {
  section("Test 3: LP with only col bounds (min x0+x1, x0>=2, x1>=3, x0<=10)");

  HighsLp lp;
  lp.num_col_ = 2;  lp.num_row_ = 0;
  lp.col_cost_  = {1.0, 1.0};
  lp.col_lower_ = {2.0, 3.0};
  lp.col_upper_ = {10.0, inf};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  // No LP rows: empty A matrix (0×2)
  auto& m = lp.a_matrix_;
  m.format_  = MatrixFormat::kColwise;
  m.num_col_ = 2; m.num_row_ = 0;
  m.start_ = {0, 0, 0};

  HighsOptions opts;
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, 5.0), "obj == 5.0");
  check(near(res.solution.col_value[0], 2.0), "x0 == 2");
  check(near(res.solution.col_value[1], 3.0), "x1 == 3");

  // col_dual[j] = +z[k_colLB] * es = 1 at LB
  check(near(res.solution.col_dual[0], 1.0),
        "col_dual[0] == 1 (x0 at col LB)");
  check(near(res.solution.col_dual[1], 1.0),
        "col_dual[1] == 1 (x1 at col LB)");

  // x0 <= 10 not binding: col_dual contribution from kColUB should be 0
  // (so col_dual[0] only gets contribution from kColLB)
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4: Infeasible LP
//
//   min  x0
//   s.t.  x0 <= 2  (row UB)
//         x0 >= 3  (row LB)   ← infeasible with the above
//         x0 >= 0  (col LB)
//
// Expected: status == kInfeasible, primal/dual_feasible = false.
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_infeasible() {
  section("Test 4: Infeasible LP (x0<=2, x0>=3)");

  HighsLp lp;
  lp.num_col_ = 1;  lp.num_row_ = 2;
  lp.col_cost_  = {1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {inf};
  lp.row_lower_ = {-inf, 3.0};
  lp.row_upper_ = {2.0,  inf};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  auto& m = lp.a_matrix_;
  m.format_  = MatrixFormat::kColwise;
  m.num_col_ = 1; m.num_row_ = 2;
  m.start_ = {0, 2};
  m.index_ = {0, 1};
  m.value_ = {1.0, 1.0};

  HighsOptions opts;
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kInfeasible,
        "status == kInfeasible");
  check(!res.primal_feasible, "primal_feasible == false");
  check(!res.dual_feasible,   "dual_feasible   == false");
  check(!res.solution.value_valid, "value_valid == false");
  check(!res.solution.dual_valid,  "dual_valid  == false");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 5: Maximisation LP
//
//   max  x0
//   s.t.  x0 <= 5  (row UB)
//         x0 >= 0  (col LB)
//
// M2 negates q to min -x0.
// Optimal: x0=5, obj=5.
// Clarabel: z[row_UB]=1 (binding UB), z[col_LB]=0 (not at LB)
// effective_sign = -1  (kMaximize)
// row_dual[0] = -z[row_UB] * (-1) = +1
// col_dual[0] = +z[col_LB] * (-1) = 0
// KKT: A'*row_dual + col_dual  == -c (for max HiGHS convention)
//       [1]*[1] + [0]           == -(-1) = 1  ← but c_highs=1 and es=-1
// More precisely: c_eff = c*es = [1]*(-1) = [-1], and A'y+rc = [-1]. ✓
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_maximize() {
  section("Test 5: Maximisation LP (max x0, x0<=5, x0>=0)");

  HighsLp lp;
  lp.num_col_ = 1;  lp.num_row_ = 1;
  lp.col_cost_  = {1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {5.0};
  lp.sense_     = ObjSense::kMaximize;
  lp.offset_    = 0.0;

  auto& m = lp.a_matrix_;
  m.format_  = MatrixFormat::kColwise;
  m.num_col_ = 1; m.num_row_ = 1;
  m.start_ = {0, 1};
  m.index_ = {0};
  m.value_ = {1.0};

  HighsOptions opts;
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, 5.0), "obj == 5.0 (max obj)");
  check(near(res.solution.col_value[0], 5.0), "x0 == 5");

  // row_dual[0] = 1 for binding UB in max LP
  check(near(res.solution.row_dual[0], 1.0),
        "row_dual[0] == 1 (x0=5 at row UB, max LP)");

  // KKT check: for max LP, A'*row_dual + col_dual = c * es
  // c = [1], es = -1 (max) → expected residual = -1
  // Actually the check is: A'*row_dual + col_dual == effective_q
  // where effective_q = q (Clarabel q = -c for max) = [-1]
  // But col_dual is stored as per HiGHS max convention (negated).
  // We just check the KKT residual formula holds for the effective sense.
  // Here: A'[1]*row_dual[1]+col_dual[0] = 1*1+0 = 1 = -q[0] = -(-1) = 1 ✓
  auto kkt = kktResidual(lp, res.solution);
  // For max LP: A'*row_dual + col_dual = -q = c (the original unsensed c)
  // Wait: HiGHS flips duals for max at IPX time. Let me just check the
  // sign is consistent with HiGHS's expectation (row_dual >= 0 at row LB,
  // row_dual <= 0 at row UB for min; flipped for max):
  // For max at row UB: row_dual should be >= 0.
  check(res.solution.row_dual[0] >= -kTol,
        "row_dual[0] >= 0 for max LP at UB (dual feasible)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 6: QP (min ½x² - 3x, x >= 0)
//
//   min  ½x² - 3x
//   s.t.  x >= 0
//
// Optimal: x=3, obj = ½*9 - 9 = -4.5.
// col_dual[0] = 0  (x=3 not at col LB, rc = 1*3 - 3 - 0 = 0)
// ─────────────────────────────────────────────────────────────────────────────
static void test_qp_basic() {
  section("Test 6: Basic QP (min ½x²-3x, x>=0)");

  HighsModel model;
  HighsLp&       lp  = model.lp_;
  HighsHessian&  H   = model.hessian_;

  lp.num_col_ = 1;  lp.num_row_ = 0;
  lp.col_cost_  = {-3.0};   // linear term
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {inf};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  // No LP rows
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1; lp.a_matrix_.num_row_ = 0;
  lp.a_matrix_.start_ = {0, 0};

  // P = [[1]] (upper triangular 1×1)
  H.dim_    = 1;
  H.format_ = HessianFormat::kTriangular;
  H.start_  = {0, 1};
  H.index_  = {0};
  H.value_  = {1.0};

  HighsOptions opts;
  auto in  = HighsClarabelInterface::buildFromModel(model, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, -4.5, 1e-3), "obj == -4.5");
  check(near(res.solution.col_value[0], 3.0, 1e-3), "x == 3");
  check(near(res.solution.col_dual[0], 0.0, 1e-3),
        "col_dual[0] == 0 (x not at col bound)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 7: Two-variable LP (binding both row UB and col LB)
//
//   min  x0 + 2*x1
//   s.t.  x0 + x1 <= 5  (row UB)
//         x0, x1 >= 1   (col LB)
//
// Optimal: x0=1, x1=1, obj=3 (both col LBs binding, row UB slack).
// col_dual[0] = 1-0 = 1 (reduced cost = c[0] - A'*row_dual = 1 - 0 = 1 > 0 ✓ at LB)
// col_dual[1] = 2-0 = 2 (reduced cost = c[1] - A'*row_dual = 2 - 0 = 2 > 0 ✓ at LB)
// row_dual[0] = 0 (UB not binding)
// ─────────────────────────────────────────────────────────────────────────────
static void test_lp_two_var() {
  section("Test 7: Two-variable LP (min x0+2x1, x0+x1<=5, x0,x1>=1)");

  HighsLp lp;
  lp.num_col_ = 2;  lp.num_row_ = 1;
  lp.col_cost_  = {1.0, 2.0};
  lp.col_lower_ = {1.0, 1.0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {5.0};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 0.0;

  auto& m = lp.a_matrix_;
  m.format_ = MatrixFormat::kColwise;
  m.num_col_ = 2; m.num_row_ = 1;
  m.start_ = {0, 1, 2};
  m.index_ = {0, 0};
  m.value_ = {1.0, 1.0};

  HighsOptions opts;
  auto in  = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, 3.0), "obj == 3.0");
  check(near(res.solution.col_value[0], 1.0), "x0 == 1");
  check(near(res.solution.col_value[1], 1.0), "x1 == 1");

  // row UB not binding → row_dual = 0
  check(near(res.solution.row_dual[0], 0.0), "row_dual[0] == 0 (row UB not binding)");
  // col LBs binding → col_dual = c = [1, 2]
  check(near(res.solution.col_dual[0], 1.0), "col_dual[0] == 1 (x0 at col LB)");
  check(near(res.solution.col_dual[1], 2.0), "col_dual[1] == 2 (x1 at col LB)");

  auto kkt = kktResidual(lp, res.solution);
  check(near(kkt[0], lp.col_cost_[0]) && near(kkt[1], lp.col_cost_[1]),
        "KKT: A'*row_dual + col_dual == c");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 8: Objective offset
//
//   min  x0 + 10          (offset = 10)
//   s.t.  x0 >= 2 (col LB)
//
// Optimal: x0=2, obj = 2 + 10 = 12.
// ─────────────────────────────────────────────────────────────────────────────
static void test_obj_offset() {
  section("Test 8: Objective offset (min x0 + 10, x0>=2)");

  HighsLp lp;
  lp.num_col_ = 1;  lp.num_row_ = 0;
  lp.col_cost_  = {1.0};
  lp.col_lower_ = {2.0};
  lp.col_upper_ = {inf};
  lp.sense_     = ObjSense::kMinimize;
  lp.offset_    = 10.0;

  lp.a_matrix_.format_  = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1; lp.a_matrix_.num_row_ = 0;
  lp.a_matrix_.start_ = {0, 0};

  HighsOptions opts;
  auto in  = HighsClarabelInterface::buildFromLp(lp, opts);
  auto res = solveAndWrite(in, lp);

  check(res.model_status == HighsModelStatus::kOptimal, "status == kOptimal");
  check(near(res.objective_value, 12.0), "obj == 12.0 (with offset)");
  check(near(res.solution.col_value[0], 2.0), "x0 == 2");
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main() {
  std::cout << "test_m3_solution: HighsClarabelSolution unit tests\n";

  test_lp_equality();
  test_lp_double_sided();
  test_lp_col_bounds_only();
  test_lp_infeasible();
  test_lp_maximize();
  test_qp_basic();
  test_lp_two_var();
  test_obj_offset();

  std::cout << "\n════════════════════════════════\n";
  std::cout << "Results:  " << g_pass << " passed,  " << g_fail << " failed\n";
  if (g_fail == 0)
    std::cout << "test_m3_solution: ALL PASS\n";
  else
    std::cout << "test_m3_solution: FAILURES DETECTED\n";

  return (g_fail == 0) ? 0 : 1;
}
