#!/usr/bin/env bash
# run_netlib_regression.sh
#
# Batch regression: compare Clarabel vs default solver on LP instances under
# check/instances/.  Skips MIP (.mps files containing "INTEGERS" section),
# QP (files containing QUADRATIC section), and files > 100 MB.
#
# Usage:
#   ./src/run_netlib_regression.sh [build_dir] [instances_dir]
#
# Defaults: build_dir=./build_clarabel, instances_dir=./check/instances

set -euo pipefail

BUILD_DIR="${1:-./build_clarabel}"
INST_DIR="${2:-./check/instances}"
HIGHS_BIN="$BUILD_DIR/bin/highs"
TOL=1e-4

if [[ ! -x "$HIGHS_BIN" ]]; then
  echo "ERROR: HiGHS binary not found at $HIGHS_BIN" >&2
  exit 1
fi

pass=0; fail=0; skip=0
failed_cases=()

for mps in "$INST_DIR"/*.mps; do
  [[ -f "$mps" ]] || continue

  # Skip large files (>100 MB)
  size=$(stat -c%s "$mps" 2>/dev/null || stat -f%z "$mps")
  if (( size > 100000000 )); then
    echo "SKIP (large) $mps"
    (( skip++ )) || true
    continue
  fi

  # Skip MIP instances
  if grep -Eqi "INTORG|INTEND|MARKER|BINARY|INTEGER" "$mps" 2>/dev/null; then
    (( skip++ )) || true
    continue
  fi

  # Skip QP instances
  if grep -Eqi "QUADOBJ|QUADRATIC|QMATRIX|QSECTION" "$mps" 2>/dev/null; then
    (( skip++ )) || true
    continue
  fi

  name=$(basename "$mps" .mps)

  # Solve with default solver
  default_out=$("$HIGHS_BIN" --solver choose "$mps" 2>&1 || true)
  default_status=$(echo "$default_out" | grep -oP "Model\s+status\s*:\s*\K[^\n]+" | head -1 || true)
  default_obj=$(echo    "$default_out" | grep -oP "Objective value\s*:\s*\K[-0-9.eE+]+" | head -1 || true)

  # Solve with Clarabel
  cla_out=$("$HIGHS_BIN" --solver clarabel "$mps" 2>&1 || true)
  cla_status=$(echo "$cla_out" | grep -oP "Model\s+status\s*:\s*\K[^\n]+" | head -1 || true)
  cla_obj=$(echo    "$cla_out" | grep -oP "Objective value\s*:\s*\K[-0-9.eE+]+" | head -1 || true)

  # Compare
  status_match=0
  obj_match=0
  [[ -n "$default_status" && -n "$cla_status" && "$default_status" == "$cla_status" ]] && status_match=1

  if [[ "$default_status" == *"Optimal"* && "$cla_status" == *"Optimal"* ]]; then
    if [[ -n "$default_obj" && -n "$cla_obj" ]]; then
      rel_err=$(python3 -c "
a,b = float('$default_obj'), float('$cla_obj')
denom = max(abs(a), abs(b), 1.0)
print(abs(a-b)/denom)
" 2>/dev/null || echo "999")
      obj_ok=$(python3 -c "print(1 if float('$rel_err') < $TOL else 0)" 2>/dev/null || echo "0")
      [[ "$obj_ok" == "1" ]] && obj_match=1
    fi
  else
    obj_match=$status_match  # for non-optimal, status agreement suffices
  fi

  if [[ $status_match -eq 1 && $obj_match -eq 1 ]]; then
    echo "PASS  $name  status=$cla_status  obj=$cla_obj"
    (( pass++ )) || true
  else
    echo "FAIL  $name  default=[$default_status/$default_obj]  clarabel=[$cla_status/$cla_obj]"
    failed_cases+=("$name")
    (( fail++ )) || true
  fi
done

echo ""
echo "=== Netlib Regression Summary ==="
echo "  PASS: $pass   FAIL: $fail   SKIP: $skip"
if (( ${#failed_cases[@]} > 0 )); then
  echo "  Failed instances:"
  for f in "${failed_cases[@]}"; do echo "    - $f"; done
fi
(( fail == 0 )) && exit 0 || exit 1
