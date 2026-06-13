#!/usr/bin/env python3
"""Submit-time TeX gate remediation (results-split spec §8.6).

Called by .github/workflows/submit.yml when the TeX parseability gate
(scripts/tests/test_tex_parseability.mjs) fails after a submission write.

Contract:
  argv[1] = path to the gate's --json report. Failure objects look like
      {file, CNI, result, item?, field, category, error, rawError, snippet}
  with `file` repo-relative (e.g. "library-bundled/<topo>/<cfg>/entry.json")
  and, for entry-level fields, `result` = integer index into Results[].

  - report missing, or report has an empty failure list: nothing to do,
    exit 0 (the workflow re-runs the gate, which then reports its own
    verdict, including setup/version-pin failures that produce no report);
  - EVERY failure is field "resultTeXPreview" inside the entry.json of the
    config directory this workflow run just wrote (recomputed from
    payload.json's cnickelIndex with the same sanitize() as the write step):
    blank those previews, set resultTeXTruncated true, write atomically and
    indent-preserving via scripts/_results_split_common.py, exit 0;
  - ANY other failure: print the failure list, exit 1 (submission rejected).

Degrading the preview is safe because the full resultTeX lives in the
sibling results.json; the UI treats an empty preview + resultTeXTruncated
as "open the detail view for the rendered result" (spec §3.4).
"""
import json
import pathlib
import sys

# Import the split primitives from the repo's scripts/ dir, anchored to this
# file's location so the helper works regardless of CWD; file operations
# below (payload.json, library-bundled/...) intentionally use the CWD, which
# the workflow sets to the repo root.
_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))
import _results_split_common as _rs  # noqa: E402


def sanitize(s):
    """MUST match sanitize() in the submit.yml write step."""
    return s.replace("|", "_").replace(":", "_")


def _print_failures(failures):
    for f in failures:
        item = f.get("item")
        loc = f"result {f.get('result')}" + (f", item {item}" if item is not None else "")
        print(f"  {f.get('category')}  {f.get('field')}  {f.get('file')} "
              f"[{loc}] CNI={f.get('CNI')}\n      error: {f.get('error')}")


def main(argv):
    if len(argv) != 2:
        print("usage: texgate_remediate.py <gate-report.json>")
        return 1
    report_path = pathlib.Path(argv[1])
    if not report_path.exists():
        print(f"texgate_remediate: no report at {report_path}; nothing to do "
              "(gate failed before producing one; its re-run will report).")
        return 0
    report = json.loads(report_path.read_text(encoding="utf-8"))
    failures = report.get("failures") or []
    if not failures:
        print("texgate_remediate: report lists no failures; nothing to do.")
        return 0

    # The config directory this run just wrote, recomputed exactly as the
    # write step does it (cnickelIndex -> topo/mass -> sanitized dirs).
    try:
        payload = json.load(open("payload.json", encoding="utf-8"))
        cnickel = payload["cnickelIndex"]
    except Exception as e:
        print(f"texgate_remediate: cannot determine the submitted config "
              f"(payload.json unreadable: {e}); refusing to remediate.")
        _print_failures(failures)
        return 1
    parts = cnickel.split(":", 1)
    topo_nickel = parts[0]
    mass_nickel = (parts[1] if len(parts) > 1
                   else "".join("0" if c != "|" else "|" for c in topo_nickel))
    entry_rel = (f"library-bundled/{sanitize(topo_nickel)}/"
                 f"{sanitize(mass_nickel)}/entry.json")

    def remediable(f):
        return (f.get("field") == "resultTeXPreview"
                and str(pathlib.PurePosixPath(f.get("file", ""))) == entry_rel)

    hard = [f for f in failures if not remediable(f)]
    if hard:
        print(f"texgate_remediate: {len(hard)} failure(s) are NOT a "
              f"resultTeXPreview in the submitted config ({entry_rel}); "
              "only that one degradation is permitted (spec §8.6). Failures:")
        _print_failures(hard)
        return 1

    entry_path = pathlib.Path(entry_rel)
    text = entry_path.read_text(encoding="utf-8")
    entry = json.loads(text)
    results = entry.get("Results", [])
    touched = set()
    for f in failures:
        idx = f.get("result")
        if not isinstance(idx, int) or not (0 <= idx < len(results)):
            print(f"texgate_remediate: failure references Results index "
                  f"{idx!r}, absent from {entry_rel}; refusing to remediate.")
            _print_failures(failures)
            return 1
        results[idx]["resultTeXPreview"] = ""
        results[idx]["resultTeXTruncated"] = True
        touched.add(idx)
    _rs._atomic_write(entry_path, entry, _rs.detect_indent(text))
    print(f"texgate_remediate: blanked resultTeXPreview and set "
          f"resultTeXTruncated on Results{sorted(touched)} of {entry_rel} "
          "(degrade-not-reject, spec §8.6); full resultTeX remains in the "
          "sibling results.json.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
