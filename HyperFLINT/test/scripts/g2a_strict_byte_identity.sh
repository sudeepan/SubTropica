#!/bin/bash
# Gate G2a (spec 2026-06-10-carry-option-design.md 4c): with
# carry_discharge:false the new binary's responses are byte-identical
# to the pre-WIP reference (29175ac8e), modulo timing_compute_s,
# schema_version, hf_version.  DEFAULT env: HF_LR_MAX_DEG must be
# unset (b48a4f154 neutrality holds only there).
#
# Usage:
#   g2a_strict_byte_identity.sh <new-hyperflint> <ref-hyperflint>
#
# All 3 fixtures were verified to parse cleanly on the reference binary
# (29175ac8e) before inclusion here — no fixture swap was required.
#
# Fixture inventory:
#   0: algebraic_letters:false, single poly, simplest clean LR case
#   1: algebraic_letters:true, single deg-2 poly (root_polys emitted)
#   2: algebraic_letters:true, 3-poly group incl. fat kinematic coeff
#      (the Task-3 shadow-face body with carry_discharge:false == Strict)
set -u
NEW=${1:?new hyperflint-cli}; REF=${2:?reference hyperflint-cli}
unset HF_LR_MAX_DEG
strip() {
  sed -E \
    's/"timing_compute_s":[0-9.eE+\-]+//g; s/"schema_version":[0-9]+//g; s/"hf_version":"[^"]*"//g'
}
declare -a BODIES=(
  '{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x*y + s*x + y"],"algebraic_letters":false,"carry_discharge":false}'
  '{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x^2*y + x*(1 + y) + y"],"algebraic_letters":true,"carry_discharge":false}'
  '{"op":"find_lr_orders","xvars":["x","y"],"coeff_vars":["s"],"polys":["x+y+1","x^2 + x*y + y + 1","y + (s^4+s^3+s^2+s+1)"],"algebraic_letters":true,"carry_discharge":false}'
)
fail=0
for i in "${!BODIES[@]}"; do
  a=$(echo "${BODIES[$i]}" | "$NEW" eval-json 2>/dev/null | grep '^{' | strip)
  b=$(echo "${BODIES[$i]}" | "$REF" eval-json 2>/dev/null | grep '^{' | strip)
  if [ -z "$a" ] || [ -z "$b" ]; then
    echo "G2A ERROR fixture $i: empty response (new=${#a}B ref=${#b}B)"
    fail=1
    continue
  fi
  if [ "$a" != "$b" ]; then
    echo "G2A FAIL fixture $i"
    echo "new: $a"
    echo "ref: $b"
    fail=1
  fi
done
[ $fail -eq 0 ] && echo "G2A PASS (${#BODIES[@]} fixtures)"
exit $fail
