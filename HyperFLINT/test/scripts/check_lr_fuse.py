#!/usr/bin/env python3
"""Issue #52 round 6 regression: the predicted-cost fuse (HF_LR_MAX_STEP_COST).

The round-6 pole face {1+x3, x3+x6, Q1 (13 terms), A' (84 terms)} in
{x1,x2,x3,x4,x6} wedged both verify_order (25 min, no verdict) and the
exhaustive find_lr_orders (30 min, no verdict) on v1.2.14, inside single
FLINT discriminants/resultants that the operand-size fuse cannot see.

  mode "verify" -> verify_order of the order that integrates this face
                   ({x1,x4,x2,x6,x3}) must RETURN (exit 0, JSON) within the
                   ctest timeout, with either order_is_lr true or
                   verify_inconclusive true; never a bare NOT-LR, never a hang.
  mode "search" -> the exhaustive search must RETURN a verdict within the
                   timeout: an order, or nolr with search_complete false (the
                   fuse skipped paths), or a structured budget abort; never a
                   hang.
  mode "small"  -> a face far below the cap must be byte-identical in
                   verdict: an order of length 2, search_complete true,
                   search_skipped_paths 0.
  mode "nofuse" -> the small face with HF_LR_MAX_STEP_COST=0 and NO time
                   budget (the fuse is inert): verdict identical to "small".
  mode "off"    -> the pole face with an EXPLICIT cap (HF_LR_MAX_STEP_COST=1e12,
                   a hundred times the default; the operations it admits
                   complete in about a second on this face) must still return
                   the same order: the knob is honoured and a looser cap does
                   not change the verdict.
"""
import json, os, subprocess, sys, time

binary, mode = sys.argv[1], sys.argv[2]
fixture = sys.argv[3] if len(sys.argv) > 3 else None

env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "OMP_NUM_THREADS": "4",
       "HF_LR_TIME_BUDGET_S": "900"}
if mode in ("verify", "search", "off"):
    req = json.load(open(fixture))
    if mode == "verify":
        req["verify_order"] = ["x1", "x4", "x2", "x6", "x3"]
else:
    req = {"op": "find_lr_orders", "schema_version_min": 1,
           "groups": [["x1", "x2", "x1+x2", "1+x1+x2"]], "xvars": ["x1", "x2"],
           "coeff_vars": []}
if mode == "off":
    env["HF_LR_MAX_STEP_COST"] = "1e12"
if mode == "nofuse":
    env["HF_LR_MAX_STEP_COST"] = "0"
    env.pop("HF_LR_TIME_BUDGET_S", None)

t0 = time.time()
p = subprocess.run([binary, "eval-json"], input=json.dumps(req), capture_output=True,
                   text=True, timeout=1500, env=env)
wall = time.time() - t0
if p.returncode != 0:
    print(f"FAIL: exit {p.returncode}; stderr tail: {p.stderr.strip().splitlines()[-3:]}")
    sys.exit(1)
try:
    resp = json.loads(p.stdout.strip())
except Exception as e:
    print(f"FAIL: non-JSON stdout ({e}): {p.stdout[:300]}")
    sys.exit(1)
print(f"mode={mode} wall={wall:.1f}s keys={sorted(resp.keys())[:12]}")

if mode == "verify":
    ok = resp.get("order_is_lr") is True or resp.get("verify_inconclusive") is True
    print(f"order_is_lr={resp.get('order_is_lr')} inconclusive={resp.get('verify_inconclusive')} "
          f"reason={resp.get('verify_reason','')[:120]!r}")
    if not ok:
        print("FAIL: verify returned neither LR nor inconclusive (a bare NOT-LR on the face that integrates in this order)")
        sys.exit(1)
elif mode in ("search", "off"):
    if resp.get("budget_exceeded") is True:
        print(f"FAIL: structured budget abort instead of a verdict: {resp.get('reason','')[:160]!r}")
        sys.exit(1)
    elif resp.get("nolr") is True:
        print(f"FAIL: nolr (search_complete={resp.get('search_complete')} skipped_paths={resp.get('search_skipped_paths')}) on the face whose order {{x1,x4,x2,x6,x3}} is known to integrate")
        sys.exit(1)
    else:
        print(f"order found: {resp.get('best_order')} score={resp.get('score')} "
              f"complete={resp.get('search_complete')} skipped={resp.get('search_skipped_paths')}")
        # The order that integrates this face (issue #52 round 6, verified end
        # to end) is the calibration target: the default cap must keep finding
        # exactly it, with the incomplete flag set (paths were skipped).
        if resp.get("best_order") != ["x1", "x4", "x2", "x6", "x3"]:
            print("FAIL: the default cap no longer returns the order that integrates this face")
            sys.exit(1)
        if resp.get("search_complete") is not False or not resp.get("search_skipped_paths"):
            print("FAIL: a search that skipped paths must report search_complete false and a positive skipped count")
            sys.exit(1)
else:
    if len(resp.get("best_order", [])) != 2 or resp.get("search_complete") is not True \
       or resp.get("search_skipped_paths") != 0 or resp.get("nolr") is True:
        print(f"FAIL: small face verdict changed: {resp}")
        sys.exit(1)
    print(f"small face: order={resp['best_order']} complete={resp['search_complete']} skipped={resp['search_skipped_paths']}")
print("PASS")
