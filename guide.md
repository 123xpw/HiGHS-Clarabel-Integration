# HiGHS + Clarabel Integration Guide

## Overview

This project integrates [Clarabel.cpp](https://github.com/oxfordcontrol/Clarabel.cpp) (a Rust/C++ interior-point conic solver) into HiGHS as an optional LP/QP solver. Two new solver options are added:

```
solver         = clarabel   # for standalone LP or QP
mip_lp_solver  = clarabel   # for root-node LP relaxation inside MIP
```

---

## Prerequisites

| Dependency | Version | Notes |
|------------|---------|-------|
| CMake | ≥ 3.18 | |
| C++ compiler | GCC ≥ 10 / Clang ≥ 12 | C++17 required |
| Rust toolchain | stable | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |
| Eigen | 3.x | bundled in `extern/eigen` — no separate install needed |

The Clarabel.cpp source is included as a git submodule under `extern/Clarabel.cpp`.

---

## Build

```bash
# Make sure Cargo is on PATH
export PATH="$HOME/.cargo/bin:$PATH"

# Configure (Clarabel off by default; turn it on explicitly)
cmake -S . -B build_clarabel \
      -DUSE_CLARABEL=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_TESTING=OFF

# Build
cmake --build build_clarabel -- -j$(nproc)
```

The binary is at `build_clarabel/bin/highs`.

> **Offline build**: if the machine has no internet access, set `CARGO_NET_OFFLINE=true` before running cmake.

---

## Running the Tests

```bash
# Build with tests enabled
cmake -S . -B build_clarabel -DUSE_CLARABEL=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_clarabel -- -j$(nproc)

./build_clarabel/bin/test_m4_integration    # 5/5 — LP, infeasible, dual, maximization, QP
./build_clarabel/bin/test_m_mip_integration # 2/2 — MIP root-node integration
./build_clarabel/bin/test_clarabel_destructive  # Netlib regression
```

---

## Usage

### C++ API

```cpp
#include "Highs.h"

Highs h;
h.readModel("problem.mps");

// Use Clarabel for LP or QP
h.setOptionValue("solver", "clarabel");
h.run();

// Use Clarabel for MIP root-node LP relaxation
h.setOptionValue("mip_lp_solver", "clarabel");
h.run();
```

### Options file (`.opt`)

```
solver clarabel
```
or
```
mip_lp_solver clarabel
```

---

## Source Files

### New files

| File | Purpose |
|------|---------|
| `src/lp_data/HighsClarabelInterface.h/.cpp` | Maps HiGHS LP/QP to Clarabel's problem format (objective, constraint matrix, cones) |
| `src/lp_data/HighsClarabelSolution.h/.cpp` | Writes Clarabel's solution back to HiGHS (primal/dual variables, status) |
| `highs/lp_data/HighsClarabelSolver.h/.cpp` | Orchestrates a single LP/QP solve: variable shift, Clarabel call, KKT quality gate |
| `cmake/clarabel-integration.cmake` | Builds the `clarabel_c` static library via Cargo and sets include paths |

### Modified files

| File | Change |
|------|--------|
| `CMakeLists.txt` / `highs/CMakeLists.txt` | `USE_CLARABEL` option; links Clarabel and Eigen |
| `highs/lp_data/Highs.cpp` | QP branch routes through `solveLpClarabel(hessian_ptr)` |
| `highs/lp_data/HighsSolve.h/.cpp` | `solveLpClarabel()` entry point; LP routing logic |
| `highs/lp_data/HighsOptions.h/.cpp` | Registers `clarabel` as a valid value for `solver` and `mip_lp_solver` |
| `highs/mip/HighsLpRelaxation.cpp` | MIP root-node: calls Clarabel for first LP relaxation, falls back to simplex on failure or missing basis |
| `highs/presolve/ICrashX.cpp` | Bug fix: resizes `x[]` to cover slack columns added for doubly-bounded rows |

---

## Key Implementation Details

### Variable shift (LP only)
Before passing a problem to Clarabel, lower bounds are shifted to zero (`x = x' + lb`). This prevents a crash in Clarabel's `NonnegCone` when the RHS is negative. The shift is reversed when writing the solution back. QP problems skip this step.

### KKT quality gate
After Clarabel returns a solution, the integration checks:
1. **Numerical safety**: rejects solutions containing non-finite values or constraint violations beyond 100× tolerance.
2. **KKT residuals** (LP only): rejects solutions where primal/dual residuals or the objective gap exceed `kkt_tolerance`.

A rejected solution falls through to HiGHS simplex for cleanup.

### MIP root-node fallback
When `mip_lp_solver = clarabel`, only the first LP relaxation at the root node uses Clarabel. Two-level fallback:
1. Clarabel returns an error → clear and re-solve with simplex.
2. Clarabel succeeds but produces no valid basis → clear and re-solve with simplex.

This guarantees that the B&B engine always receives a valid simplex basis for warm-starting subsequent nodes.

---

## Benchmark Results

### Test setup
- **Problem**: Security-Constrained Unit Commitment (SCUC) MIP
- **Instances**: 12 monthly instances per grid (random seed 42, one date per month, year 2017)
- **Settings**: `time_limit = 3600 s`, `mip_rel_gap = 1%`
- **Solvers**: HiPO (HiGHS simplex), IPX (HiGHS interior-point), Clarabel (this integration)

### case118 — 36 periods, ~40k variables

| Month | IPX (s) | HiPO (s) | Clarabel (s) |
|-------|---------|---------|-------------|
| Jan | 7.16 | 10.45 | **5.46** |
| Feb | 7.26 | 7.87 | **5.53** |
| Mar | 7.60 | 8.66 | **5.74** |
| Apr | 9.50 | 11.50 | **6.47** |
| May | 6.69 | 7.33 | **5.65** |
| Jun | 8.36 | 9.88 | **5.56** |
| Jul | 10.57 | 12.26 | **6.37** |
| Aug | 7.55 | 8.35 | **6.08** |
| Sep† | 31.17 | 38.16 | **31.87** |
| Oct | 6.48 | 8.25 | **5.79** |
| Nov | 7.22 | 9.20 | **6.09** |
| Dec | 6.60 | 9.11 | **6.04** |
| **Mean (excl. Sep)** | **7.73** | **9.35** | **5.89** |
| **Optimal / 12** | 12/12 | 12/12 | **12/12** |

†Sep is a hard instance (initial LP gap ≈ 2.38%); all three solvers take ~31–38 s.  
Excluding it, Clarabel is **24% faster than IPX** and **37% faster than HiPO**.

### case2383wp — 72 periods, ~767k variables

| Month | IPX (s) | HiPO (s) | Clarabel (s) |
|-------|---------|---------|-------------|
| Jan | 1030.87 | 902.11 | **590.07** |
| Feb | 679.50 | 518.81 | **410.65** |
| Mar | 592.68 | 818.68 | **441.35** |
| Apr | **359.11** | 642.53 | 453.73 |
| May | 362.97 | 297.53 | **267.19** |
| Jun | 1088.70 | 1308.41 | **712.49** |
| Jul | 1203.11 | 1044.50 | **612.92** |
| Aug | ❌ 3600 (2.5% gap) | ❌ 3600 (2.6% gap) | ✅ **1894.25** |
| Sep | 651.22 | 605.91 | **414.85** |
| Oct | ❌ 7173‡ (45.7% gap) | ❌ 3602 (1.1% gap) | ✅ **2500.35** |
| Nov | 1562.17 | ❌ 8169‡ (51.4% gap) | ✅ **783.95** |
| Dec | 1388.26 | ❌ 3600 (1.3% gap) | ✅ **736.33** |
| **Mean (s)** | **1641** | **2093** | **818** |
| **Optimal / 12** | **10/12** | **8/12** | **12/12** |

‡Root LP relaxation itself exceeded the time limit before any B&B node was processed.

**Summary**: Clarabel is the **only solver to achieve 100% solve rate** on case2383wp, with a mean solve time ~2× faster than IPX and ~2.6× faster than HiPO. On the two hardest instances (Aug, Oct) where both IPX and HiPO failed, Clarabel found the optimal solution within the time limit.
