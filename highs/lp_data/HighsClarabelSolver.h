/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsClarabelSolver.h
 * @brief M4 — Clarabel LP solver integration for HiGHS.
 *
 * Provides solveLpClarabel(), which converts a HighsLp to Clarabel's input
 * format (M2 logic), solves via clarabel::DefaultSolver<double>, recovers
 * the HiGHS solution (M3 logic), and optionally runs IPX crossover to obtain
 * a simplex basis.
 *
 * Compiled only when HIGHS_USE_CLARABEL is defined (i.e., cmake -DUSE_CLARABEL=ON).
 */
#ifndef LP_DATA_HIGHS_CLARABEL_SOLVER_H_
#define LP_DATA_HIGHS_CLARABEL_SOLVER_H_

#ifdef HIGHS_USE_CLARABEL

#include "lp_data/HighsLpSolverObject.h"
#include "model/HighsHessian.h"

// Solve the LP (or QP) in solver_object.lp_ with Clarabel.
//
// hessian_ptr: when non-null, the quadratic term 1/2 x^T H x is included
//   (QP mode). When null, P = 0 and the problem is treated as LP.
//
// On success (optimal / almost-optimal):
//   solver_object.solution_.value_valid = true
//   solver_object.solution_.dual_valid  = true
//   solver_object.model_status_         = kOptimal
//   solver_object.basis_.valid          = true  (after crossover, if enabled)
//
// On infeasibility / error:
//   model_status_ is set to kInfeasible / kUnbounded / kSolveError / etc.
//   solution_.value_valid = false
//
// run_crossover option is respected: "on" or "choose" → calls callCrossover().
HighsStatus solveLpClarabel(HighsLpSolverObject& solver_object,
                             const HighsHessian* hessian_ptr = nullptr);

#endif  // HIGHS_USE_CLARABEL
#endif  // LP_DATA_HIGHS_CLARABEL_SOLVER_H_
