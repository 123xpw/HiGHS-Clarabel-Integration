// src/test_m4_integration.cpp
//
// M4 integration test: routes LP through solver = "clarabel" in Highs::run().
// Tests:
//   1. Simple LP (min x, x >= 1) — optimal
//   2. Infeasible LP (x <= -1, x >= 1)
//   3. Dual feasibility via the shadow price of an equality constraint
//   4. Maximization LP

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "Highs.h"

static bool approxEq(double a, double b, double tol = 1e-4) {
  return std::fabs(a - b) <= tol;
}

static double getObj(const Highs& h) {
  double v = 0.0;
  h.getInfoValue("objective_function_value", v);
  return v;
}

// ── Test 1: simple LP  min x  s.t. x >= 1 ──────────────────────────────────
static void test_simple_lp() {
  Highs h;
  h.setOptionValue("output_flag", false);
  h.addVar(1.0, 1e30);  // x with lower bound 1
  h.changeColCost(0, 1.0);
  h.setOptionValue("solver", "clarabel");
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution& sol = h.getSolution();
  assert(approxEq(sol.col_value[0], 1.0));
  double obj = getObj(h);
  assert(approxEq(obj, 1.0));
  std::printf("  test_simple_lp: PASS (x=%.4g obj=%.4g)\n",
              sol.col_value[0], obj);
}

// ── Test 2: infeasible LP ───────────────────────────────────────────────────
static void test_infeasible_lp() {
  Highs h;
  h.setOptionValue("output_flag", false);
  h.addVar(-1e30, 1e30);
  h.changeColCost(0, 1.0);
  // x <= -1
  {
    std::vector<int>    idx = {0};
    std::vector<double> val = {1.0};
    h.addRow(-1e30, -1.0, 1, idx.data(), val.data());
  }
  // x >= 1
  {
    std::vector<int>    idx = {0};
    std::vector<double> val = {1.0};
    h.addRow(1.0, 1e30, 1, idx.data(), val.data());
  }
  h.setOptionValue("solver", "clarabel");
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kInfeasible);
  std::printf("  test_infeasible_lp: PASS (status = Infeasible)\n");
}

// ── Test 3: equality constraint shadow price ────────────────────────────────
static void test_dual() {
  Highs h;
  h.setOptionValue("output_flag", false);
  // min x + 2y  s.t. x + y == 3, x >= 0, y >= 0
  // Optimal: x=3, y=0, obj=3.
  h.addVar(0.0, 1e30);
  h.addVar(0.0, 1e30);
  h.changeColCost(0, 1.0);
  h.changeColCost(1, 2.0);
  {
    std::vector<int>    idx = {0, 1};
    std::vector<double> val = {1.0, 1.0};
    h.addRow(3.0, 3.0, 2, idx.data(), val.data());
  }
  h.setOptionValue("solver", "clarabel");
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution& sol = h.getSolution();
  assert(approxEq(sol.col_value[0], 3.0, 1e-3));
  assert(approxEq(sol.col_value[1], 0.0, 1e-3));
  // shadow price of equality = 1 (cost increase for tighter rhs)
  assert(approxEq(std::fabs(sol.row_dual[0]), 1.0, 0.1));
  double obj = getObj(h);
  assert(approxEq(obj, 3.0, 1e-3));
  std::printf("  test_dual: PASS (obj=%.4g row_dual=%.4g)\n",
              obj, sol.row_dual[0]);
}

// ── Test 4: maximization LP ─────────────────────────────────────────────────
static void test_max_lp() {
  Highs h;
  h.setOptionValue("output_flag", false);
  // max x  s.t. x <= 5
  h.addVar(-1e30, 5.0);
  h.changeColCost(0, 1.0);
  h.changeObjectiveSense(ObjSense::kMaximize);
  h.setOptionValue("solver", "clarabel");
  h.run();

  assert(h.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution& sol = h.getSolution();
  assert(approxEq(sol.col_value[0], 5.0, 1e-3));
  double obj = getObj(h);
  assert(approxEq(obj, 5.0, 1e-3));
  std::printf("  test_max_lp: PASS (x=%.4g obj=%.4g)\n",
              sol.col_value[0], obj);
}

// ── main ────────────────────────────────────────────────────────────────────
int main() {
  std::printf("=== M4 Integration Tests (solver = clarabel) ===\n");
  test_simple_lp();
  test_infeasible_lp();
  test_dual();
  test_max_lp();
  std::printf("=== All M4 tests PASSED ===\n");
  return 0;
}
