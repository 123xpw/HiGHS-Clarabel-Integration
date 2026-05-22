// src/test_m_mip_integration.cpp
//
// M_MIP integration test: root-node LP relaxation solved by Clarabel
// (mip_lp_solver = "clarabel"), subsequent B&B by HiGHS simplex.
//
// Tests:
//   1. Simple MIP:  min x + y  s.t. x + y >= 3, x,y integer, x,y >= 0
//      Root LP relaxation optimal via Clarabel; B&B closes to integer solution.
//   2. MIP with forced Clarabel at root: verify same optimal as default solver.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "Highs.h"

static bool approxEq(double a, double b, double tol = 0.5) {
  return std::fabs(a - b) <= tol;
}

static double getObj(const Highs& h) {
  double v = 0.0;
  h.getInfoValue("objective_function_value", v);
  return v;
}

// Solve a MIP with the given solver setting and return its optimal value.
static double solveMip(const std::string& mip_lp_solver) {
  Highs h;
  h.setOptionValue("output_flag", false);

  // min x + y  s.t. x + y >= 3, x,y integer >= 0
  // Optimal: x=3, y=0 (or x=0, y=3, etc.), obj = 3
  h.addVar(0.0, 1e30);   // x
  h.addVar(0.0, 1e30);   // y
  h.changeColCost(0, 1.0);
  h.changeColCost(1, 1.0);
  h.changeColIntegrality(0, HighsVarType::kInteger);
  h.changeColIntegrality(1, HighsVarType::kInteger);
  {
    std::vector<int>    idx = {0, 1};
    std::vector<double> val = {1.0, 1.0};
    h.addRow(3.0, 1e30, 2, idx.data(), val.data());  // x + y >= 3
  }
  h.setOptionValue("mip_lp_solver", mip_lp_solver);
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kOptimal);
  return getObj(h);
}

static void test_mip_clarabel_vs_default() {
  const double obj_clarabel = solveMip("clarabel");
  const double obj_default  = solveMip("choose");

  std::printf("  MIP obj (clarabel root): %.4g\n", obj_clarabel);
  std::printf("  MIP obj (default):       %.4g\n", obj_default);

  assert(approxEq(obj_clarabel, 3.0));
  assert(approxEq(obj_default,  3.0));
  assert(approxEq(obj_clarabel, obj_default));
  std::printf("  test_mip_clarabel_vs_default: PASS\n");
}

// Slightly larger MIP to stress-test the root-node flow.
static void test_mip_larger() {
  Highs h;
  h.setOptionValue("output_flag", false);

  // min  sum(c_i * x_i)  s.t. sum(x_i) >= 5, x_i in {0,1}
  // Optimal obj = 5 * 1 = 5  (5 cheapest binary vars set to 1)
  const int n = 8;
  for (int i = 0; i < n; ++i) {
    h.addVar(0.0, 1.0);
    h.changeColCost(i, 1.0 + i * 0.1);
    h.changeColIntegrality(i, HighsVarType::kInteger);
  }
  {
    std::vector<int>    idx(n);
    std::vector<double> val(n, 1.0);
    for (int i = 0; i < n; ++i) idx[i] = i;
    h.addRow(5.0, 1e30, n, idx.data(), val.data());
  }
  h.setOptionValue("mip_lp_solver", "clarabel");
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kOptimal);
  double obj = getObj(h);
  // cheapest 5 vars: costs 1.0, 1.1, 1.2, 1.3, 1.4 → total = 6.0
  assert(approxEq(obj, 6.0, 0.1));
  std::printf("  test_mip_larger: PASS (obj=%.4g)\n", obj);
}

int main() {
  std::printf("=== M_MIP Integration Tests (mip_lp_solver = clarabel) ===\n");
  test_mip_clarabel_vs_default();
  test_mip_larger();
  std::printf("=== All M_MIP tests PASSED ===\n");
  return 0;
}
