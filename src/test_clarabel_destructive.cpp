// src/test_clarabel_destructive.cpp
//
// Destructive end-to-end tests for the public HiGHS Clarabel integration.
// These tests intentionally use only the Highs API so they exercise the real
// solver routing, status propagation, and MIP root-LP path.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "Highs.h"

namespace {

int failures = 0;

bool approx(double a, double b, double tol = 1e-5) {
  return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

void check(bool ok, const char* msg) {
  if (ok) {
    std::printf("  PASS  %s\n", msg);
  } else {
    std::printf("  FAIL  %s\n", msg);
    ++failures;
  }
}

double objective(const Highs& h) {
  double obj = 0.0;
  h.getInfoValue("objective_function_value", obj);
  return obj;
}

void test_negative_column_bounds() {
  std::printf("\n-- negative finite column bounds --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("solver", "clarabel");

  h.addVar(-3.0, 10.0);
  h.addVar(-5.0, 10.0);
  h.changeColCost(0, 1.0);
  h.changeColCost(1, 1.0);

  const HighsStatus run_status = h.run();
  const HighsSolution& sol = h.getSolution();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kOptimal,
        "model status is optimal");
  check(approx(sol.col_value[0], -3.0, 1e-4), "x0 reaches negative lower");
  check(approx(sol.col_value[1], -5.0, 1e-4), "x1 reaches negative lower");
  check(approx(objective(h), -8.0, 1e-4), "objective includes negative bounds");
}

void test_negative_row_bounds() {
  std::printf("\n-- negative row bounds --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("solver", "clarabel");

  h.addVar(-10.0, 10.0);
  h.changeColCost(0, 1.0);
  {
    std::vector<int> idx = {0};
    std::vector<double> val = {1.0};
    h.addRow(-3.0, -3.0, 1, idx.data(), val.data());
  }

  const HighsStatus run_status = h.run();
  const HighsSolution& sol = h.getSolution();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kOptimal,
        "model status is optimal");
  check(approx(sol.col_value[0], -3.0, 1e-4), "x satisfies negative equality");
  check(approx(sol.row_value[0], -3.0, 1e-4),
        "row activity preserves negative rhs");
  check(approx(objective(h), -3.0, 1e-4), "objective is negative rhs");
}

void test_unbounded_status() {
  std::printf("\n-- unbounded LP status --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("solver", "clarabel");

  h.addVar(0.0, kHighsInf);
  h.changeColCost(0, -1.0);

  const HighsStatus run_status = h.run();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kUnbounded,
        "dual-infeasible Clarabel result maps to unbounded");
}

void test_free_rows_do_not_constrain() {
  std::printf("\n-- free rows are ignored by mapping --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("solver", "clarabel");

  h.addVar(2.0, 10.0);
  h.changeColCost(0, 1.0);
  {
    std::vector<int> idx = {0};
    std::vector<double> val = {100.0};
    h.addRow(-kHighsInf, kHighsInf, 1, idx.data(), val.data());
  }

  const HighsStatus run_status = h.run();
  const HighsSolution& sol = h.getSolution();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kOptimal,
        "model status is optimal");
  check(approx(sol.col_value[0], 2.0, 1e-4),
        "free row does not move bound optimum");
}

void test_qp_clarabel_gap_is_detectable() {
  std::printf("\n-- QP solver=clarabel currently falls back to QP ASM --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("solver", "clarabel");

  h.addVar(0.0, kHighsInf);
  h.changeColCost(0, -3.0);
  HighsInt q_start[2] = {0, 1};
  HighsInt q_index[1] = {0};
  double q_value[1] = {1.0};
  check(h.passHessian(1, 1, static_cast<HighsInt>(HessianFormat::kTriangular),
                      q_start, q_index, q_value) == HighsStatus::kOk,
        "pass convex one-variable Hessian");

  const HighsStatus run_status = h.run();
  const HighsInfo& info = h.getInfo();
  const HighsSolution& sol = h.getSolution();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kOptimal,
        "model status is optimal");
  check(info.qp_iteration_count > 0,
        "QP ASM was used; solver=clarabel is not a QP Clarabel route yet");
  check(approx(sol.col_value[0], 3.0, 1e-4), "QP optimum x=3");
  check(approx(objective(h), -4.5, 1e-4), "QP objective includes Hessian");
}

void test_mip_root_clarabel() {
  std::printf("\n-- MIP root LP via mip_lp_solver=clarabel --\n");
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("mip_lp_solver", "clarabel");

  h.addVar(0.0, 1.0);
  h.addVar(0.0, 1.0);
  h.addVar(0.0, 1.0);
  h.changeColCost(0, 1.0);
  h.changeColCost(1, 2.0);
  h.changeColCost(2, 4.0);
  h.changeColIntegrality(0, HighsVarType::kInteger);
  h.changeColIntegrality(1, HighsVarType::kInteger);
  h.changeColIntegrality(2, HighsVarType::kInteger);
  {
    std::vector<int> idx = {0, 1, 2};
    std::vector<double> val = {1.0, 1.0, 1.0};
    h.addRow(2.0, kHighsInf, 3, idx.data(), val.data());
  }

  const HighsStatus run_status = h.run();
  check(run_status == HighsStatus::kOk, "run status is kOk");
  check(h.getModelStatus() == HighsModelStatus::kOptimal,
        "model status is optimal");
  check(approx(objective(h), 3.0, 1e-4), "MIP optimum after root crossover");
}

}  // namespace

int main() {
  std::printf("test_clarabel_destructive: public API end-to-end checks\n");
  test_negative_column_bounds();
  test_negative_row_bounds();
  test_unbounded_status();
  test_free_rows_do_not_constrain();
  test_qp_clarabel_gap_is_detectable();
  test_mip_root_clarabel();

  std::printf("\nResults: %s\n", failures == 0 ? "ALL PASS" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
