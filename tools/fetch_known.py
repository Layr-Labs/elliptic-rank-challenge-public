#!/usr/bin/env python3
"""Maintainer tool: fetch high-rank curves from ICARM, Dujella and LMFDB into a raw mirror.

This is step 1 of the four-step record-store pipeline (docs/records.md):
  1. tools/fetch_known.py → tools/known_raw.json (raw mirror v2: ainvs, points, source torsion)
  2. tools/certify_known.py → tools/known_raw.json (verification block per entry)
  3. tools/build_known.py → ec_verify/known_curves.json (ledger v2: fingerprints, torsion, H, status)
  4. tools/build_records.py → ec_verify/records.json (per-(torsion, rank) records and baselines)

Sources:
  * ICARM Elliptic Curve Rank Leaderboard: https://elliptic-rank.icarm.cloud/database.json
  * Andrej Dujella's record pages in tools/cache: the rank pages rk*.html / rkeq*.html and the
    per-torsion pages z*.html with their "old records" pages (https://web.math.pmf.unizg.hr/~duje/tors/),
    every curve on every page, parsed by tools/fetch_records.py:parse_curves
  * the LMFDB: the classic anchor curves (tools/anchor_curves.json) and the smallest curves of every
    (torsion, rank) from Cremona's ecdata dump with their generators (tools/lmfdb_curves.json,
    written by tools/fetch_lmfdb.py)

The fingerprint step is in tools/build_known.py; this tool never runs gp.

    python tools/fetch_known.py [--raw tools/known_raw.json]
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from ec_verify.cells import MAZUR_GROUPS  # noqa: E402  (pure Python: fetch_known never runs gp)
from fetch_records import parse_curves  # noqa: E402

ICARM_DB = "https://elliptic-rank.icarm.cloud/database.json"
TOOLS = ROOT / "tools"
RAW = TOOLS / "known_raw.json"
CACHE = TOOLS / "cache"
_PAGE = re.compile(r"(rk(?:eq)?\d+|z\d+(?:z\d+)?(?:old\d+)?)\.html")   # record pages, not tors/rankhist/generic

_SOURCE_ORDER = {"icarm": 0, "dujella": 1, "lmfdb": 2}   # also the tie-break order in build_known / build_records

# Curves a Dujella page misprints so badly that no listed point fits any reading of the equation,
# keyed by the a-invariants of the most plausible reading (what parse_curves reports). Each must be
# accounted for by hand (where the real curve is) or the fetch aborts: a dropped published curve is
# exactly a curve the gate would then accept.
DUJELLA_MISPRINTS = {
    # z1.html curve 13, Ritik Jain (2026), rank 31: the equation as printed; the curve is icarm-3856
    (0, 1, 0, -396137825789612441729444338863306159122947430920364402299145,
     3512038234266144362533916241901734617982932157407905278925467439474016818118945979916226837),
}


def _curl(url: str, timeout: int = 30) -> str | None:
    try:
        p = subprocess.run(["curl", "-fsSL", "-m", str(timeout), url], capture_output=True, text=True)
    except OSError:
        return None
    return p.stdout if p.returncode == 0 and p.stdout.strip() else None


def _ints_from_bracketed(s: str):
    """'[1, 1, 1, -12...]' -> [1, 1, 1, -12...]; entries may be 'a/b'."""
    body = s.strip().strip("[]")
    out = []
    for tok in body.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "/" in tok:
            n, d = tok.split("/")
            out.append((int(n), int(d)))
        else:
            out.append(int(tok))
    return out


class BadRecord(ValueError):
    """An ICARM record we cannot turn into a curve. Never skipped silently: a dropped published
    curve is exactly a curve the gate would then accept."""


def _points(d, cid):
    pts = d.get("points") or []
    out = []
    for p in pts:
        if not (isinstance(p, (list, tuple)) and len(p) == 2):
            raise BadRecord(f"id {cid}: point {p!r} is not a pair")
        out.append([str(p[0]).replace(" ", ""), str(p[1]).replace(" ", "")])
    return out


def _exact_int(a):
    """An a-invariant from a JSON list: an int or a decimal string, exactly. Anything else (a
    float, a bool, a fraction) is returned unchanged so the five-integers check rejects it:
    int(-2.5) == -2 would mirror a different curve than the published one."""
    if type(a) is int:
        return a
    if isinstance(a, str):
        try:
            return int(a.strip())
        except ValueError:
            return a
    return a


def icarm_row(d: dict) -> dict:
    if "id" not in d:
        raise BadRecord(f"record without id: {str(d)[:120]}")
    cid = d["id"]
    ainvs = d.get("ainvs")
    if isinstance(ainvs, str):
        try:
            ainvs = _ints_from_bracketed(ainvs)
        except ValueError as e:
            raise BadRecord(f"id {cid}: ainvs {ainvs!r}: {e}") from None
    if isinstance(ainvs, list):
        ainvs = [_exact_int(a) for a in ainvs]
    if not isinstance(ainvs, list) or len(ainvs) != 5 or not all(type(a) is int for a in ainvs):
        raise BadRecord(f"id {cid}: ainvs {d.get('ainvs')!r} is not five integers (non-integral models must be minimalised by hand)")
    tors = d.get("torsion")
    if tors is not None:
        if not isinstance(tors, list) or tuple(int(x) for x in tors) not in MAZUR_GROUPS:
            raise BadRecord(f"id {cid}: torsion {tors!r} is not a Mazur group")
        tors = [int(x) for x in tors]
    nh = d.get("naive_height")
    return {"id": f"icarm-{cid}", "source": "icarm", "source_id": str(cid), "ainvs": ainvs, "points": _points(d, cid),
            "rank_claimed": int(d.get("rank_lower_bound") or 0), "torsion_source": tors,
            "naive_height_source": None if nh is None else float(nh),
            "source_created_at": d.get("created_at"), "source_updated_at": d.get("updated_at")}


def icarm_records(db) -> list:
    """The list of records in database.json: either the document itself or its single list value."""
    if isinstance(db, list):
        return db
    lists = [v for v in db.values() if isinstance(v, list)] if isinstance(db, dict) else []
    if len(lists) != 1:
        keys = sorted(db) if isinstance(db, dict) else type(db).__name__
        raise SystemExit(f"!! unexpected database shape (expected one list): {keys}")
    return lists[0]


def icarm_rows(records: list) -> list:
    rows, bad = [], []
    for d in records:
        try:
            rows.append(icarm_row(d))
        except BadRecord as e:
            bad.append(str(e))
    for b in bad:
        print(f"!! icarm: {b}", file=sys.stderr)
    if bad:
        msg = f"!! icarm: {len(bad)} record(s) could not be parsed; nothing written"
        print(msg, file=sys.stderr)
        raise SystemExit(msg)
    rows.sort(key=lambda r: int(r["source_id"]))
    return rows


def fetch_icarm():
    raw = _curl(ICARM_DB, timeout=300)
    if raw is None:
        raise SystemExit(f"!! could not download {ICARM_DB}")
    rows = icarm_rows(icarm_records(json.loads(raw)))
    print(f"icarm: {len(rows)} curves (ids up to {rows[-1]['id'] if rows else 'none'})")
    return rows


def dujella_hash_id(stem: str, ainvs) -> str:
    """`dujella-<page>-<8 hex of sha1(a-invariants)>`: stable when the page is re-fetched with new
    curves prepended (Dujella lists the newest first), unlike a position on the page."""
    h = hashlib.sha1(",".join(str(int(a)) for a in ainvs).encode()).hexdigest()[:8]
    return f"dujella-{stem}-{h}"


def _dujella_row(cid: str, ainvs, points, rank: int) -> dict:
    return {"id": cid, "source": "dujella", "source_id": cid[len("dujella-"):], "ainvs": list(ainvs),
            "points": [[str(x), str(y)] for x, y in points], "rank_claimed": rank, "torsion_source": None,
            "naive_height_source": None, "source_created_at": None, "source_updated_at": None}


def dujella_rows():
    """Every curve on every cached Dujella record page. The first curve of each rank page keeps its
    historical id `dujella-rk20` (rank from the file name); every other curve, on the rank pages and
    on the per-torsion pages, gets a hash id and its announced point count as the claimed rank. A
    curve printed on several pages is mirrored once (its historical row if it has one, else the
    first page in name order). A page that does not parse aborts the run; a misprinted point is
    reported and left out; a misprinted equation aborts unless listed in DUJELLA_MISPRINTS."""
    rows, seen, bad, problems = [], set(), [], []
    pages = sorted(p for p in CACHE.glob("*.html") if _PAGE.fullmatch(p.name))
    parsed = {}
    for p in pages:
        try:
            parsed[p] = parse_curves(p.read_text(errors="replace"), problems)
        except Exception as e:
            bad.append(f"{p.name}: {e}")
            continue
        for pr in problems:
            print(f"!! dujella: {p.name}: {pr['message']}", file=sys.stderr)
            if pr["kind"] == "omitted" and tuple(pr["ainvs"]) not in DUJELLA_MISPRINTS:
                bad.append(f"{p.name}: {pr['message']} (not in DUJELLA_MISPRINTS)")
        problems.clear()
    legacy = {}                                        # a-invariants -> (page, rank) of every historical row
    for p, curves in parsed.items():
        m = re.match(r"rk(eq)?(\d+)\.html", p.name)
        if m and curves:
            legacy.setdefault(tuple(curves[0]["ainvs"]), (p, int(m.group(2))))
    for p, curves in parsed.items():
        for c in curves:
            key = tuple(c["ainvs"])
            if key in seen:
                continue
            if key in legacy:
                lp, rank = legacy[key]
                if lp != p:
                    continue                           # written under its historical id when its own page comes up
                rows.append(_dujella_row(f"dujella-{lp.stem}", c["ainvs"], c["points"][:rank], rank))
            else:
                rows.append(_dujella_row(dujella_hash_id(p.stem, c["ainvs"]), c["ainvs"], c["points"], c["listed"]))
            seen.add(key)
    for b in bad:
        print(f"!! dujella: {b}", file=sys.stderr)
    if bad:
        msg = f"!! dujella: {len(bad)} problem(s); nothing written"
        print(msg, file=sys.stderr)
        raise SystemExit(msg)
    print(f"dujella: {len(rows)} curves from {len(pages)} pages in tools/cache")
    return rows


def anchor_rows():
    """The LMFDB source: the classic anchors (no points) and, when tools/lmfdb_curves.json exists,
    the smallest curves per (torsion, rank) from Cremona's ecdata with their generators and
    published torsion. The same label in both takes the ecdata row (it has the points)."""
    rows = {}
    for a in json.loads((TOOLS / "anchor_curves.json").read_text()):
        if a["source"] == "lmfdb":
            rows[a["id"]] = {"id": f"lmfdb-{a['id']}", "source": "lmfdb", "source_id": a["id"], "ainvs": list(a["ainvs"]),
                             "points": [], "rank_claimed": a["rank"], "torsion_source": None, "naive_height_source": None,
                             "source_created_at": None, "source_updated_at": None}
    p = TOOLS / "lmfdb_curves.json"
    if p.exists():
        for c in json.loads(p.read_text())["curves"]:
            rows[c["label"]] = {"id": f"lmfdb-{c['label']}", "source": "lmfdb", "source_id": c["label"],
                                "ainvs": [int(a) for a in c["ainvs"]], "points": [[str(x), str(y)] for x, y in c["gens"]],
                                "rank_claimed": int(c["rank"]), "torsion_source": [int(t) for t in c["torsion"]],
                                "naive_height_source": None, "source_created_at": None, "source_updated_at": None}
    print(f"lmfdb: {len(rows)} curves (anchors + tools/lmfdb_curves.json)")
    return list(rows.values())


def _sort_key(c):
    sid = c["source_id"]
    num = int(sid) if sid.isdigit() else 10 ** 9
    return (_SOURCE_ORDER.get(c["source"], 9), num, sid)


CONTENT_FIELDS = ("ainvs", "points", "rank_claimed", "torsion_source")


def _content(entry: dict) -> tuple:
    """What certification depends on. Dujella and LMFDB rows carry no updated_at, so a refreshed
    page is detected by content, not by timestamp."""
    return tuple(entry.get(k) for k in CONTENT_FIELDS)


def merge(existing: dict, rows: list, now: str) -> dict:
    """Append-only merge: keep all existing entries, update/add from rows. An entry is kept as is
    (verification included) only when its source timestamp and its content are both unchanged."""
    old = {c["id"]: c for c in existing.get("curves", [])}
    seen = set()
    out = []
    for r in rows:
        seen.add(r["id"])
        prev = old.get(r["id"])
        if prev is not None and prev.get("source_updated_at") == r.get("source_updated_at") and _content(prev) == _content(r):
            out.append(prev)                   # unchanged: keep its verification block
            continue
        entry = dict(r)
        entry["first_seen"] = prev.get("first_seen", now) if prev is not None else now
        out.append(entry)
    # Append-only: keep entries that vanished from upstream
    for cid, entry in old.items():
        if cid not in seen:
            out.append(entry)
    out.sort(key=_sort_key)
    return {"version": 2, "fetched_at": now, "curves": out}


def load_raw(path: pathlib.Path = RAW) -> dict:
    return json.loads(path.read_text()) if path.exists() else {"version": 2, "curves": []}


def dump_per_line(doc: dict, list_key: str) -> str:
    """JSON with one compact record per line inside doc[list_key] and the other top-level keys
    sorted, one per line. Still plain JSON for json.loads; a refresh that adds a curve changes one
    line, so the ledger (4k records) and the raw mirror stay reviewable in a pull request. A
    json.dumps with indent put every fingerprint entry on its own line: 816k lines."""
    head = [f"{json.dumps(k)}: {json.dumps(doc[k], sort_keys=True)}," for k in sorted(doc) if k != list_key]
    items = [json.dumps(item, sort_keys=True, separators=(",", ":")) for item in doc[list_key]]
    body = ",\n".join(items)
    return "{\n" + "\n".join(head) + "\n" + f"{json.dumps(list_key)}: [\n" + body + ("\n" if items else "") + "]\n}\n"


def save_raw(path: pathlib.Path, raw: dict) -> None:
    path.write_text(dump_per_line(raw, "curves"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=pathlib.Path, default=RAW)
    args = ap.parse_args()
    rows = fetch_icarm() + dujella_rows() + anchor_rows()
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    raw = merge(load_raw(args.raw), rows, now)
    save_raw(args.raw, raw)
    print(f"wrote {args.raw}: {len(raw['curves'])} curves")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
