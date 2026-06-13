#!/usr/bin/env python3
"""validate_results_split.py — alignment validator + de-migration guard
(results-split spec docs/superpowers/specs/2026-06-11-results-split-design.md
§8.1/§8.2).

Exit 1 if any library-bundled config violates the split invariants:

  DISARMED (default, pre-migration CI mode): split-pair consistency only —
    wherever a split exists (a record carries resultDataId, or a sibling
    results.json is present), stub ids <-> sibling keys must be in bijection,
    stubs must carry no heavy fields, and ids must be unique. Purely
    monolithic entries (heavy fields inline, no sibling) PASS, so the
    still-monolithic library is green before Task 10 migrates it.

  ARMED (--armed, post-migration CI mode; Task 10 flips the flag in the
    migration commit): additionally, ANY record carrying a heavy field
    without resultDataId fails ANYWHERE — the de-migration guard, so a
    re-run legacy script writing monolithic records fails CI here.

Heavy-field-free records (vacuum-period, singularities, proposed-alphabet)
pass in both modes by construction: the heavy-field scan finds no key.

Usage:
  python3 scripts/library_audit/validate_results_split.py           # disarmed
  python3 scripts/library_audit/validate_results_split.py --armed   # guard on
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _results_split_common as rs

armed = "--armed" in sys.argv
problems = []
n_configs = 0
for entry_path in sorted((ROOT / "library-bundled").glob("*/*/entry.json")):
    cfg = entry_path.parent
    n_configs += 1
    problems += rs.validate_split_pair(cfg)   # bijection + stubs-carry-no-heavy
    if armed:
        try:
            entry = json.loads(entry_path.read_text(encoding="utf-8"))
        except Exception:
            continue   # already reported as unparseable by validate_split_pair
        for i, rec in enumerate(entry.get("Results", [])):
            if "resultDataId" not in rec:
                for k in rs.HEAVY_RESULT_FIELDS:
                    if k in rec:
                        problems.append(
                            f"{entry_path}: Results[{i}] monolithic heavy field {k} (de-migration?)")
for p in problems:
    print("FAIL", p)
print(f"checked {n_configs} configs ({'armed' if armed else 'disarmed'}); "
      f"{len(problems)} problem(s)")
sys.exit(1 if problems else 0)
