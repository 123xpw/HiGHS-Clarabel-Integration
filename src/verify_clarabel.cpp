// Verify that the Clarabel -> C++ -> HiGHS build chain is correctly wired.
// Compile with -DUSE_CLARABEL=ON.  No actual solve is required to pass.

#include <clarabel.hpp>

#include <Eigen/Eigen>
#include <iostream>
#include <vector>

int main() {
    using namespace clarabel;
    using namespace Eigen;

    // Minimal LP:  min  x  s.t.  x >= 1   (x in R)
    // Clarabel form:  min  q'x  s.t.  Ax + s = b,  s in NonnegativeCone
    //   P = 0,  q = [1],  A = [-1],  b = [-1],  cone = Nonneg(1)

    SparseMatrix<double> P(1, 1);  // zero (LP)

    VectorXd q(1);
    q << 1.0;

    SparseMatrix<double> A(1, 1);
    A.insert(0, 0) = -1.0;
    A.makeCompressed();

    VectorXd b(1);
    b << -1.0;

    std::vector<SupportedConeT<double>> cones{NonnegativeConeT<double>(1)};

    DefaultSettings<double> settings = DefaultSettings<double>::default_settings();
    settings.verbose = false;

    DefaultSolver<double> solver(P, q, A, b, cones, settings);
    solver.solve();

    auto sol = solver.solution();
    std::cout << "verify_clarabel: status = " << static_cast<int>(sol.status)
              << "  obj = " << sol.obj_val
              << "  x[0] = " << sol.x[0] << "\n";

    // SolverStatus::Solved == 1
    bool ok = (sol.status == SolverStatus::Solved) &&
              (std::abs(sol.x[0] - 1.0) < 1e-4);
    if (ok) {
        std::cout << "verify_clarabel: PASS - Rust->C++->HiGHS chain functional\n";
        return 0;
    } else {
        std::cout << "verify_clarabel: FAIL - unexpected solver result\n";
        return 1;
    }
}
