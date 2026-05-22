// src/test_m2_mapping.cpp
//
// Unit tests for HighsClarabelInterface (M2 data-mapping layer).
// Verifies constraint expansion, row_map contents, b values, and A matrix
// entries WITHOUT invoking the Clarabel solver.
//
// Test LP (4 rows × 3 cols):
//
//   min  x0 + x1 + x2
//   s.t.
//     Row 0 (equality):     2x0 +  x1        = 4      lo=hi=4
//     Row 1 (UB only):       x0        + x2  ≤ 5      lo=-inf, hi=5
//     Row 2 (LB only):             x1  + x2  ≥ 1      lo=1,   hi=+inf
//     Row 3 (double-sided):  x0  + x1  + x2  ∈ [1,3]  lo=1,   hi=3
//
//   Column bounds:
//     x0: (-inf, 10]   → kColUB only
//     x1: [2,   +inf)  → kColLB only
//     x2: [0,   8]     → kColUB + kColLB
//
// Expected Clarabel rows (9 total):
//   [0]   kRowEq  / ZeroCone        — Row 0 equality
//   [1]   kRowUB  / NonnegCone      — Row 1 UB
//   [2]   kRowLB  / NonnegCone      — Row 2 LB  (negated)
//   [3]   kRowUB  / NonnegCone      — Row 3 double UB
//   [4]   kRowLB  / NonnegCone      — Row 3 double LB (negated)
//   [5]   kColUB  / NonnegCone      — x0 ≤ 10
//   [6]   kColLB  / NonnegCone      — x1 ≥ 2  (negated: -x1 ≥ -2)
//   [7]   kColUB  / NonnegCone      — x2 ≤ 8
//   [8]   kColLB  / NonnegCone      — x2 ≥ 0  (negated: -x2 ≥ 0)

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "lp_data/HighsClarabelInterface.h"

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

// ── helper: build the test HighsLp ───────────────────────────────────────────

static HighsLp buildTestLp() {
  const double inf = std::numeric_limits<double>::infinity();

  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 4;

  lp.col_cost_  = {1.0, 1.0, 1.0};
  lp.col_lower_ = {-inf, 2.0, 0.0};
  lp.col_upper_ = {10.0, inf, 8.0};
  lp.row_lower_ = {4.0, -inf, 1.0, 1.0};
  lp.row_upper_ = {4.0,  5.0,  inf, 3.0};

  lp.sense_  = ObjSense::kMinimize;
  lp.offset_ = 0.0;

  // Constraint matrix (colwise CSC):
  //        x0  x1  x2
  // Row 0:  2   1   0
  // Row 1:  1   0   1
  // Row 2:  0   1   1
  // Row 3:  1   1   1
  //
  // Col 0: (0,2),(1,1),(3,1)
  // Col 1: (0,1),(2,1),(3,1)
  // Col 2: (1,1),(2,1),(3,1)
  auto& m = lp.a_matrix_;
  m.format_  = MatrixFormat::kColwise;
  m.num_col_ = 3;
  m.num_row_ = 4;
  m.start_ = {0, 3, 6, 9};
  m.index_ = {0, 1, 3,   0, 2, 3,   1, 2, 3};
  m.value_ = {2.0, 1.0, 1.0,   1.0, 1.0, 1.0,   1.0, 1.0, 1.0};

  return lp;
}

// ── main ─────────────────────────────────────────────────────────────────────

int main() {
  using RowKind = HighsClarabelInterface::RowKind;
  using Tag     = clarabel::SupportedConeT<double>::Tag;

  std::cout << "test_m2_mapping: HighsClarabelInterface unit tests\n";

  HighsLp lp = buildTestLp();
  HighsOptions opts;  // default settings
  auto in = HighsClarabelInterface::buildFromLp(lp, opts);

  // ── 1. Dimensions ───────────────────────────────────────────────────────────
  // 1 eq + 1 UB + 1 LB + 2 (double-sided) = 5 LP rows
  // x0→1 UB + x1→1 LB + x2→1 UB+1 LB    = 4 col-bound rows
  // total = 9

  section("Dimensions");
  check(in.A.rows()         == 9, "A.rows() == 9");
  check(in.A.cols()         == 3, "A.cols() == 3");
  check(in.b.size()         == 9, "b.size() == 9");
  check((int)in.cones.size()   == 9, "cones.size() == 9");
  check((int)in.row_map.size() == 9, "row_map.size() == 9");
  check(in.num_highs_rows   == 4, "num_highs_rows == 4");
  check(in.num_highs_cols   == 3, "num_highs_cols == 3");

  // ── 2. row_map kinds and HiGHS indices ──────────────────────────────────────
  section("row_map kinds");
  check(in.row_map[0].kind == RowKind::kRowEq && in.row_map[0].highs_idx == 0,
        "row[0]: kRowEq, highs_idx=0");
  check(in.row_map[1].kind == RowKind::kRowUB && in.row_map[1].highs_idx == 1,
        "row[1]: kRowUB, highs_idx=1");
  check(in.row_map[2].kind == RowKind::kRowLB && in.row_map[2].highs_idx == 2,
        "row[2]: kRowLB, highs_idx=2");
  check(in.row_map[3].kind == RowKind::kRowUB && in.row_map[3].highs_idx == 3,
        "row[3]: kRowUB, highs_idx=3 (double-sided UB)");
  check(in.row_map[4].kind == RowKind::kRowLB && in.row_map[4].highs_idx == 3,
        "row[4]: kRowLB, highs_idx=3 (double-sided LB)");
  check(in.row_map[5].kind == RowKind::kColUB && in.row_map[5].highs_idx == 0,
        "row[5]: kColUB, highs_idx=0 (x0 ≤ 10)");
  check(in.row_map[6].kind == RowKind::kColLB && in.row_map[6].highs_idx == 1,
        "row[6]: kColLB, highs_idx=1 (x1 ≥ 2)");
  check(in.row_map[7].kind == RowKind::kColUB && in.row_map[7].highs_idx == 2,
        "row[7]: kColUB, highs_idx=2 (x2 ≤ 8)");
  check(in.row_map[8].kind == RowKind::kColLB && in.row_map[8].highs_idx == 2,
        "row[8]: kColLB, highs_idx=2 (x2 ≥ 0)");

  // ── 3. Cone types ──────────────────────────────────────────────────────────
  section("Cone types");
  check(in.cones[0].tag == Tag::ZeroConeT,        "cones[0]: ZeroCone (equality)");
  check(in.cones[1].tag == Tag::NonnegativeConeT, "cones[1]: NonnegCone (UB)");
  check(in.cones[2].tag == Tag::NonnegativeConeT, "cones[2]: NonnegCone (LB)");
  check(in.cones[3].tag == Tag::NonnegativeConeT, "cones[3]: NonnegCone (double UB)");
  check(in.cones[4].tag == Tag::NonnegativeConeT, "cones[4]: NonnegCone (double LB)");
  check(in.cones[5].tag == Tag::NonnegativeConeT, "cones[5]: NonnegCone (kColUB)");
  check(in.cones[6].tag == Tag::NonnegativeConeT, "cones[6]: NonnegCone (kColLB)");
  check(in.cones[7].tag == Tag::NonnegativeConeT, "cones[7]: NonnegCone (kColUB)");
  check(in.cones[8].tag == Tag::NonnegativeConeT, "cones[8]: NonnegCone (kColLB)");

  // ── 4. b values ─────────────────────────────────────────────────────────────
  // b[0] =  4.0  (equality lo=hi=4)
  // b[1] =  5.0  (UB hi=5)
  // b[2] = -1.0  (LB: -lo = -1, negated because row was flipped)
  // b[3] =  3.0  (double-sided UB hi=3)
  // b[4] = -1.0  (double-sided LB: -lo = -1)
  // b[5] = 10.0  (kColUB: cu[0]=10)
  // b[6] = -2.0  (kColLB: -cl[1]=-2)
  // b[7] =  8.0  (kColUB: cu[2]=8)
  // b[8] =  0.0  (kColLB: -cl[2]=0)

  section("b vector");
  check(in.b[0] ==  4.0, "b[0] ==  4.0  (equality)");
  check(in.b[1] ==  5.0, "b[1] ==  5.0  (UB)");
  check(in.b[2] == -1.0, "b[2] == -1.0  (LB, negated)");
  check(in.b[3] ==  3.0, "b[3] ==  3.0  (double-sided UB)");
  check(in.b[4] == -1.0, "b[4] == -1.0  (double-sided LB, negated)");
  check(in.b[5] == 10.0, "b[5] == 10.0  (kColUB x0)");
  check(in.b[6] == -2.0, "b[6] == -2.0  (kColLB x1, negated)");
  check(in.b[7] ==  8.0, "b[7] ==  8.0  (kColUB x2)");
  check(in.b[8] ==  0.0, "b[8] ==  0.0  (kColLB x2, -0)");

  // ── 5. A matrix entries ─────────────────────────────────────────────────────
  // LP rows (sign expansion):
  //   Clarabel row 0 (kRowEq,  sign=+1): col0→2, col1→1
  //   Clarabel row 1 (kRowUB,  sign=+1): col0→1, col2→1
  //   Clarabel row 2 (kRowLB,  sign=-1): col1→-1, col2→-1
  //   Clarabel row 3 (double UB, sign=+1): col0→1, col1→1, col2→1
  //   Clarabel row 4 (double LB, sign=-1): col0→-1, col1→-1, col2→-1
  // Column bound rows:
  //   row 5 (kColUB x0): A[5,0]=+1
  //   row 6 (kColLB x1): A[6,1]=-1
  //   row 7 (kColUB x2): A[7,2]=+1
  //   row 8 (kColLB x2): A[8,2]=-1

  section("A matrix (LP constraint rows)");
  check(in.A.coeff(0,0) ==  2.0, "A[0,0]==2  (equality, 2*x0)");
  check(in.A.coeff(0,1) ==  1.0, "A[0,1]==1  (equality,   x1)");
  check(in.A.coeff(0,2) ==  0.0, "A[0,2]==0  (equality, no x2)");
  check(in.A.coeff(1,0) ==  1.0, "A[1,0]==1  (UB row, x0)");
  check(in.A.coeff(1,1) ==  0.0, "A[1,1]==0  (UB row, no x1)");
  check(in.A.coeff(1,2) ==  1.0, "A[1,2]==1  (UB row, x2)");
  check(in.A.coeff(2,1) == -1.0, "A[2,1]==-1 (LB row, x1 negated)");
  check(in.A.coeff(2,2) == -1.0, "A[2,2]==-1 (LB row, x2 negated)");
  check(in.A.coeff(3,0) ==  1.0, "A[3,0]==1  (double UB, x0)");
  check(in.A.coeff(3,1) ==  1.0, "A[3,1]==1  (double UB, x1)");
  check(in.A.coeff(3,2) ==  1.0, "A[3,2]==1  (double UB, x2)");
  check(in.A.coeff(4,0) == -1.0, "A[4,0]==-1 (double LB, x0 negated)");
  check(in.A.coeff(4,1) == -1.0, "A[4,1]==-1 (double LB, x1 negated)");
  check(in.A.coeff(4,2) == -1.0, "A[4,2]==-1 (double LB, x2 negated)");

  section("A matrix (column bound rows)");
  check(in.A.coeff(5,0) ==  1.0, "A[5,0]==+1 (kColUB x0)");
  check(in.A.coeff(6,1) == -1.0, "A[6,1]==-1 (kColLB x1)");
  check(in.A.coeff(7,2) ==  1.0, "A[7,2]==+1 (kColUB x2)");
  check(in.A.coeff(8,2) == -1.0, "A[8,2]==-1 (kColLB x2)");
  // Cross-check: col-bound rows touch only one column each
  check(in.A.coeff(5,1) == 0.0 && in.A.coeff(5,2) == 0.0,
        "row 5 is sparse (only col 0)");
  check(in.A.coeff(6,0) == 0.0 && in.A.coeff(6,2) == 0.0,
        "row 6 is sparse (only col 1)");

  // ── 6. Cost vector (minimise: no sign flip) ─────────────────────────────────
  section("Cost vector");
  check(in.q[0] == 1.0 && in.q[1] == 1.0 && in.q[2] == 1.0,
        "q == [1,1,1] for minimise");

  // ── 7. Metadata ────────────────────────────────────────────────────────────
  section("Metadata");
  check(in.sense      == ObjSense::kMinimize, "sense == kMinimize");
  check(in.obj_offset == 0.0,                 "obj_offset == 0");
  check(in.P.rows()   == 3 && in.P.nonZeros() == 0,
        "P is 3×3 zero (LP, no Hessian)");

  // ── Summary ─────────────────────────────────────────────────────────────────
  std::cout << "\n════════════════════════════════\n";
  std::cout << "Results:  " << g_pass << " passed,  " << g_fail << " failed\n";
  if (g_fail == 0)
    std::cout << "test_m2_mapping: ALL PASS\n";
  else
    std::cout << "test_m2_mapping: FAILURES DETECTED\n";

  return (g_fail == 0) ? 0 : 1;
}
