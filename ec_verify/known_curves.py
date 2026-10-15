"""The frozen known-curve ledger (schema v2, ec_verify/known_curves.json): the match lookup and
the grouping into isogeny classes per rank."""

from __future__ import annotations

import json
import pathlib

from . import fingerprint as fp_mod
from .curve import naive_height

_SET = None


def load(path: pathlib.Path | None = None) -> dict:
    global _SET
    if path is None and _SET is not None:
        return _SET
    p = path or pathlib.Path(__file__).with_name("known_curves.json")
    t = json.loads(p.read_text())
    if path is None:
        _SET = t
    return t


def match(fp: list, table: dict | None = None):
    """Return the known-curve record matching the fingerprint, or None."""
    t = table or load()
    for c in t["curves"]:
        if fp_mod.matches(fp, c["fingerprint"]):
            return c
    return None


def classes(table: dict | None = None) -> dict:
    """Group the known set by (fingerprint, rank_claimed): one key per isogeny class per claimed rank.

    Isogenous curves share conductor and L-function, so their |a_p| fingerprints coincide and
    the gate already treats them as one known curve; a quadratic twist usually has another
    rank or other bad primes and stays separate. Anything that counts "known curves of rank r"
    must use these classes, or a class uploaded whole counts once per member (798 of the 4048
    records in the 2026-10-05 freeze are such repeats). Members keep the table's order."""
    t = table or load()
    out: dict = {}
    for c in t["curves"]:
        out.setdefault((tuple(c["fingerprint"]), int(c["rank_claimed"])), []).append(c)
    return out


def _id_key(cid: str):
    head, _, tail = cid.rpartition("-")
    return (head, int(tail)) if tail.isdigit() else (cid, -1)


def representative(members: list) -> dict:
    """The class member a participant could have submitted: smallest naive height, then the
    numerically smallest id (icarm-9 before icarm-10)."""
    return min(members, key=lambda c: (naive_height([int(a) for a in c["ainvs"]]), _id_key(c["id"])))


_STATUSES = {"certified", "failed", "unverified"}
_REQUIRED = ("id", "source", "source_id", "ainvs", "fingerprint", "rank_claimed", "torsion", "torsion_key",
             "naive_height", "verification", "duplicate_of_earlier")


def validate_ledger(table: dict) -> None:
    if table.get("version") != 2:
        raise ValueError("ledger version must be 2")
    for c in table["curves"]:
        for f in _REQUIRED:
            if f not in c:
                raise ValueError(f"{c.get('id')}: missing {f}")
        v = c["verification"]
        if not isinstance(v, dict) or "status" not in v or "certified_rank" not in v:
            raise ValueError(f"{c['id']}: verification block lacks status or certified_rank")
        if v["status"] not in _STATUSES:
            raise ValueError(f"{c['id']}: bad verification status {v['status']!r}")
        if v["status"] == "certified" and not isinstance(v["certified_rank"], int):
            raise ValueError(f"{c['id']}: certified without certified_rank")
        if v["status"] != "certified" and v["certified_rank"] is not None:
            raise ValueError(f"{c['id']}: certified_rank set without status certified")
        if len(c["ainvs"]) != 5 or len(c["fingerprint"]) != len(fp_mod.PRIMES):
            raise ValueError(f"{c['id']}: ainvs/fingerprint length")
