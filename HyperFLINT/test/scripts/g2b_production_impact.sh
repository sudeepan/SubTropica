#!/bin/bash
# Gate G2b (spec 2026-06-10-carry-option-design.md 4c): default-request
# (carry_discharge field ABSENT) diff, new binary (absent=>false) vs
# deployed v1.2.3 (absent=>true), modulo schema_version / hf_version /
# timing fields.  Every residual diff must be carry-attributable:
# NOLR restored, a different best_order selected, or carry profile
# fields absent from the new strict response.  Output: one line per
# fixture for notes/carry_option/G2B_FINDINGS.md.
#
# Usage:
#   g2b_production_impact.sh <new-hyperflint> [<deployed-hyperflint>]
#
# Fixture inventory:
#   multilinear-2var    : db1 LR face, no algebraic -- control (IDENTICAL expected)
#   deg2-conic-2var     : db2 algebraic, trivial deg-2, clean carry -- control (IDENTICAL expected)
#   shadow-3poly        : db2 algebraic, shadow face (Task-3 fixture stripped) -- DIFF class(ii)+(iii) expected
#   gauged-uq5          : db2 algebraic, real NOLR<->LR flip (gate-1 fixture stripped) -- DIFF class(i) expected
#   synthetic-flip      : db2 algebraic, small 2-var synthetic flip -- DIFF class(i) expected
#   lr-db1-control      : db1, no algebraic_letters key at all -- IDENTICAL expected
#   lr-db2-clean        : db2 algebraic, 3 clean linear polys -- DIFF class(ii)+(iii) expected (carry DFS order flip, carried_sqrts==0)
set -u
NEW=${1:?new cli}
DEP=${2:-$HOME/Library/Wolfram/Paclets/Repository/HyperFLINT-1.2.3/dist/macos-arm64/hyperflint}
unset HF_LR_MAX_DEG
strip() { sed -E 's/"timing_compute_s":[0-9.e+\-]+//g; s/"schema_version":[0-9]+//g; s/"hf_version":"[^"]*"//g'; }

declare -a NAMES=() BODIES=()

# Control fixtures (algebraic_letters absent or trivial -- IDENTICAL expected)
NAMES+=("lr-db1-no-algebraic-control")
BODIES+=('{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x*y + s*x + y"]}')

NAMES+=("multilinear-2var")
BODIES+=('{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x*y + s*x + y"],"algebraic_letters":true}')

NAMES+=("deg2-conic-2var")
BODIES+=('{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x^2*y + x*(1 + y) + y"],"algebraic_letters":true}')

# Carry-flip fixtures: deployed absent=>true (carry), new absent=>false (strict)

# Fixture (a): gauged-uq5 with carry_discharge removed -- known NOLR<->LR flip
# Source: test_find_lr_orders_carry_discharge.cpp kGaugedUq5{Polys,Xvars,Coeffs}
NAMES+=("gauged-uq5-absent")
BODIES+=('{"op":"find_lr_orders","xvars":["x1","x2","x4","x5"],"coeff_vars":["qq1","qq2","wb1","wb2","yb"],"polys":["x1+x2+1","-qq1*x1*x2-qq2*x1+2*wb1*x4-x4^2+2*wb2*x2*x5-x5^2+2*yb*x4*x5"],"algebraic_letters":true}')

# Fixture (b): synthetic flip fixture with carry_discharge removed
# Source: test_find_lr_orders_carry_discharge.cpp kSynth{Polys,Xvars,Coeffs}
NAMES+=("synthetic-flip-absent")
BODIES+=('{"op":"find_lr_orders","xvars":["x1","x2"],"coeff_vars":["s"],"polys":["x1+x2+1","x1^2+x2^2+x1*x2+s"],"algebraic_letters":true}')

# Fixture (c): FindRoots=False body (algebraic_letters absent) -- pure control
NAMES+=("nolr-db1-control")
BODIES+=('{"op":"find_lr_orders","xvars":["x1","x2"],"polys":["1+x1^2+x2^2"]}')

# Fixture (d1): strict-LR db2 face, 3 clean linear polys
# Source: test_find_lr_orders_strategy_roundtrip.cpp lr_polys + algebraic_letters:true (row3 body)
# NB: deployed runs carry DFS (absent=>true) and may return a different order at equal score.
NAMES+=("lr-db2-clean-polys")
BODIES+=('{"op":"find_lr_orders","xvars":["x1","x2"],"polys":["x1","x2","x1+x2"],"algebraic_letters":true}')

# Fixture (d2): shadow-3poly absent -- Task-3 fixture stripped of carry_discharge:false
# Source: g2a_strict_byte_identity.sh fixture 2, carry_discharge field removed.
# Deployed (absent=>true) returns carry-cheapest order (x,y) score~43 with carried_sqrts:1;
# new (absent=>false) returns strict-cheapest (y,x) score~100 (no carry path taken).
NAMES+=("shadow-3poly-absent")
BODIES+=('{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x+y+1","x^2 + x*y + y + 1","y + (s^4+s^3+s^2+s+1)"],"algebraic_letters":true}')

for i in "${!BODIES[@]}"; do
  a=$(echo "${BODIES[$i]}" | "$NEW" eval-json 2>/dev/null | grep '^{' | strip)
  d=$(echo "${BODIES[$i]}" | "$DEP" eval-json 2>/dev/null | grep '^{' | strip)
  if [ -z "$a" ] || [ -z "$d" ]; then
    echo "${NAMES[$i]}: ERROR (empty response: new=${#a}B dep=${#d}B)"
    continue
  fi
  if [ "$a" = "$d" ]; then
    echo "${NAMES[$i]}: IDENTICAL"
  else
    echo "${NAMES[$i]}: DIFF"
    echo "  new(strict-default): $a"
    echo "  deployed(carry-default): $d"
  fi
done
