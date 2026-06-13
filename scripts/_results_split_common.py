"""Shared primitives for the entry.json / results.json split
(results-split spec docs/superpowers/specs/2026-06-11-results-split-design.md §3).

HEAVY_RESULT_FIELDS is defined in THREE languages; keep in sync:
  - here (Python)
  - SubTropica.wl  ($STHeavyResultFields)
  - ui/app.js      (HEAVY_RESULT_FIELDS)
"""
import hashlib
import json
import os
import pathlib
import tempfile

HEAVY_RESULT_FIELDS = (
    "resultCompressed", "resultTeX", "symbolTeX", "normalizedSymbolTeX",
    "wDefinitions", "algebraicLetters", "resultInputForm",
)
SPLIT_MINTED_FIELDS = ("resultDataId", "resultTeXPreview", "resultTeXTruncated")
PREVIEW_MAX_BYTES = 16384


class SplitError(RuntimeError):
    pass


def mint_result_data_id(heavy_record, taken):
    """Opaque per-writer id: sha256 of this writer's canonical JSON of the
    heavy record, first 16 hex; '<sha16>.<n>' on within-entry collision.
    Cross-writer byte agreement is NOT required (spec §3.2): readers resolve
    by the stored id, never by re-hashing."""
    base = hashlib.sha256(
        json.dumps(heavy_record, sort_keys=True, separators=(",", ":"),
                   ensure_ascii=False).encode("utf-8")).hexdigest()[:16]
    if base not in taken:
        return base
    n = 1
    while f"{base}.{n}" in taken:
        n += 1
    return f"{base}.{n}"


def split_record(record, taken_ids, max_preview_bytes=PREVIEW_MAX_BYTES):
    """Return (stub, heavy_or_None). Records with no heavy fields (period /
    singularities / proposed-alphabet / already-split stubs) come back
    unchanged with heavy=None. resultTeX <= cap is copied into
    resultTeXPreview; oversized -> empty preview + resultTeXTruncated flag
    (kernel preview generation is a separate, optional upgrade pass)."""
    heavy = {k: record[k] for k in HEAVY_RESULT_FIELDS if k in record}
    if not heavy or "resultDataId" in record:
        return dict(record), None
    stub = {k: v for k, v in record.items() if k not in HEAVY_RESULT_FIELDS}
    rid = mint_result_data_id(heavy, taken_ids)
    stub["resultDataId"] = rid
    if "resultTeXPreview" not in stub:
        tex = record.get("resultTeX", "")
        if len(tex.encode("utf-8")) <= max_preview_bytes:
            stub["resultTeXPreview"] = tex
        else:
            stub["resultTeXPreview"] = ""
            stub["resultTeXTruncated"] = True
    return stub, heavy


def merge_record(stub, sibling_results):
    """Inverse of split_record for one record. Reader convention (§3.3):
    no resultDataId -> self-contained, returned as-is."""
    rid = stub.get("resultDataId")
    if rid is None:
        return dict(stub)
    if rid not in sibling_results:
        raise SplitError(f"resultDataId {rid} missing from sibling results.json")
    merged = dict(stub)
    merged.update(sibling_results[rid])
    return merged


def detect_indent(text):
    """Indent style of an existing JSON file: first indented line wins.
    Defaults to tab (the submit.yml / json.dump(indent='\\t') convention)."""
    for line in text.splitlines():
        stripped = line.lstrip(" \t")
        if stripped and stripped != line:
            ws = line[: len(line) - len(stripped)]
            return "\t" if "\t" in ws else ws
    return "\t"


def _atomic_write(path, obj, indent):
    path = pathlib.Path(path)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".tmp")
    try:
        os.fchmod(fd, 0o644)   # mkstemp yields 0600; library convention is 0644
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=indent, ensure_ascii=False)
            f.write("\n")
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def write_split_entry(entry_path, entry):
    """Write entry.json (stubs) + sibling results.json (heavy map), both
    atomically (temp + rename), preserving each file's existing indent.
    Idempotent on already-split input. Merges into an existing sibling
    (never drops ids it does not know about)."""
    entry_path = pathlib.Path(entry_path)
    sib_path = entry_path.parent / "results.json"
    sibling = {"SchemaVersion": 1,
               "CNickelIndex": entry.get("CNickelIndex", ""),
               "Results": {}}
    sib_indent = "\t"
    if sib_path.exists():
        old = sib_path.read_text(encoding="utf-8")
        sib_indent = detect_indent(old)
        sibling = json.loads(old)
        sibling.setdefault("Results", {})
    taken = set(sibling["Results"])
    new_results = []
    changed_sibling = False
    for rec in entry.get("Results", []):
        stub, heavy = split_record(rec, taken)
        new_results.append(stub)
        if heavy is not None:
            sibling["Results"][stub["resultDataId"]] = heavy
            taken.add(stub["resultDataId"])
            changed_sibling = True
    out = dict(entry)
    out["Results"] = new_results
    entry_indent = "\t"
    if entry_path.exists():
        entry_indent = detect_indent(entry_path.read_text(encoding="utf-8"))
    _atomic_write(entry_path, out, entry_indent)
    if changed_sibling or (sib_path.exists() and not sibling["Results"]):
        _atomic_write(sib_path, sibling, sib_indent)
    return out


def read_sibling(config_dir):
    sib = pathlib.Path(config_dir) / "results.json"
    if not sib.exists():
        return {}
    return json.loads(sib.read_text(encoding="utf-8")).get("Results", {})


def read_merged_entry(entry_path):
    """entry.json with heavy fields merged back from the sibling (in memory
    only; files untouched)."""
    entry_path = pathlib.Path(entry_path)
    entry = json.loads(entry_path.read_text(encoding="utf-8"))
    sib = read_sibling(entry_path.parent)
    entry["Results"] = [merge_record(r, sib) for r in entry.get("Results", [])]
    return entry


def validate_split_pair(config_dir):
    """Alignment validator (spec §8.1). Returns a list of problem strings
    (empty = clean): stub ids and sibling keys must be in bijection; stubs
    must carry no heavy fields; id-less records need no sibling entry."""
    config_dir = pathlib.Path(config_dir)
    problems = []
    entry_path = config_dir / "entry.json"
    if not entry_path.exists():
        return problems
    try:
        entry = json.loads(entry_path.read_text(encoding="utf-8"))
    except Exception as e:
        return [f"{entry_path}: unparseable ({e})"]
    sib = read_sibling(config_dir)
    stub_ids = []
    for i, rec in enumerate(entry.get("Results", [])):
        rid = rec.get("resultDataId")
        if rid is not None:
            stub_ids.append(rid)
            for k in HEAVY_RESULT_FIELDS:
                if k in rec:
                    problems.append(f"{entry_path}: Results[{i}] is a stub but carries heavy field {k}")
            if rid not in sib:
                problems.append(f"{entry_path}: Results[{i}] id {rid} missing from results.json")
    for rid in sib:
        if rid not in stub_ids:
            problems.append(f"{config_dir / 'results.json'}: orphan id {rid} (no stub references it)")
    if len(set(stub_ids)) != len(stub_ids):
        problems.append(f"{entry_path}: duplicate resultDataId among stubs")
    return problems


def move_config(src_dir, dst_dir):
    """Sanctioned config relocation: moves entry.json AND results.json
    together (spec §7). Merge semantics when dst exists, strict by design:
      - src record identical to an existing dst record (and, for stubs,
        identical heavy under the same id): dedup, skip it;
      - same resultDataId but differing heavy OR differing stub: SplitError
        (manual resolution; never silently overwrite, never duplicate ids);
      - otherwise append.
    The merged destination is validated logically before any write; source
    files are removed only after both destination writes succeed, so a
    failure never loses data (worst case: src intact, dst untouched)."""
    src_dir, dst_dir = pathlib.Path(src_dir), pathlib.Path(dst_dir)
    src_entry_path = src_dir / "entry.json"
    dst_entry_path = dst_dir / "entry.json"
    dst_dir.mkdir(parents=True, exist_ok=True)
    if not dst_entry_path.exists():
        src_entry_path.rename(dst_entry_path)
        if (src_dir / "results.json").exists():
            (src_dir / "results.json").rename(dst_dir / "results.json")
        return dst_dir
    src_entry = json.loads(src_entry_path.read_text(encoding="utf-8"))
    src_sib = read_sibling(src_dir)
    dst_entry = json.loads(dst_entry_path.read_text(encoding="utf-8"))
    dst_results = dst_entry.setdefault("Results", [])
    # full sibling dict preserved (extra top-level keys, indent)
    dst_sib_path = dst_dir / "results.json"
    if dst_sib_path.exists():
        dst_sib_text = dst_sib_path.read_text(encoding="utf-8")
        dst_sib_full = json.loads(dst_sib_text)
        dst_sib_indent = detect_indent(dst_sib_text)
    else:
        dst_sib_full = {"SchemaVersion": 1,
                        "CNickelIndex": dst_entry.get("CNickelIndex", ""),
                        "Results": {}}
        dst_sib_indent = "\t"
    dst_sib = dst_sib_full.setdefault("Results", {})
    changed_sib = False
    for rec in src_entry.get("Results", []):
        rid = rec.get("resultDataId")
        if rid is None:
            if rec not in dst_results:
                dst_results.append(rec)
            continue
        if rid in dst_sib:
            if dst_sib[rid] != src_sib.get(rid):
                raise SplitError(f"id collision with differing heavy data: {rid}")
            twins = [d for d in dst_results if d.get("resultDataId") == rid]
            if rec in twins:
                continue                       # exact dup: skip
            raise SplitError(f"id {rid} present at dst with a differing stub")
        if rid not in src_sib:
            raise SplitError(f"src stub id {rid} has no heavy data in src sibling")
        dst_results.append(rec)
        dst_sib[rid] = src_sib[rid]
        changed_sib = True
    ids = [r["resultDataId"] for r in dst_results if "resultDataId" in r]
    if len(set(ids)) != len(ids):
        raise SplitError("merge would duplicate resultDataId among dst stubs")
    _atomic_write(dst_entry_path, dst_entry,
                  detect_indent(dst_entry_path.read_text(encoding="utf-8")))
    if changed_sib or dst_sib:
        _atomic_write(dst_sib_path, dst_sib_full, dst_sib_indent)
    src_entry_path.unlink()
    if (src_dir / "results.json").exists():
        (src_dir / "results.json").unlink()
    try:
        src_dir.rmdir()
    except OSError:
        pass   # other files remain; leave the directory
    return dst_dir
