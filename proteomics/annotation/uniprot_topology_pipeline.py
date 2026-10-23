#!/usr/bin/env python3
"""
Membrane topology from UniProt feature annotations (Fig 4E, S5C)
================================================================

Classifies each protein by its curated UniProt TRANSMEM / SIGNAL / TOPO_DOM
features and flags MitoCarta 3.0 proteins. The paper used UniProt release
2026_01 (the release and run parameters are written to <output>.meta.txt).

Pipeline:
  1. Parse UniProt-style accessions from a text file.
  2. Batch-retrieve feature JSON from UniProt's REST search endpoint.
  3. Parse TRANSMEM / SIGNAL / TOPO_DOM features.
  4. Cross-reference accessions against MitoCarta 3.0.
  5. Classify topology (TA, signal-anchored / Type II, Type I,
     multipass +/- SP, beta-barrel, soluble, secreted/lumenal).
  6. Write annotated TSV.

REQUIREMENTS
------------
    pip install pandas requests openpyxl xlrd

USAGE
-----
    python export_accessions.py --outdir <dir>     # silac_uniprotIDs.csv, total_uniprotIDS.csv
    python uniprot_topology_pipeline.py \\
        --input <dir>/total_uniprotIDS.csv \\
        --mitocarta Human.MitoCarta3.0.xls \\
        --output total_topology.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100  # accessions per UniProt request
FIELDS = "accession,length,ft_transmem,ft_signal,ft_topo_dom,ft_intramem"

ACCESSION_RE = re.compile(
    r"(?:^|[|>\s])"
    r"("
    r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4}"
    r")"
    r"(?:-\d+)?"           # optional isoform suffix, e.g. -2, -13
    r"(?:[|\s]|$)"
)

# Detects an isoform suffix in the original input token, for reporting only.
ISOFORM_SUFFIX_RE = re.compile(
    r"("
    r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4}"
    r")-(\d+)"
)


# ============================================================
# 1. Parse IDs
# ============================================================

def parse_accession(text: str) -> str | None:
    """Return the canonical 6/10-char UniProt accession, stripping any
    -N isoform suffix and any surrounding sp|...|... wrapping."""
    if not text:
        return None
    m = ACCESSION_RE.search(text)
    if m:
        return m.group(1)
    token = text.strip().split()[0] if text.strip() else ""
    # Strip optional -N isoform suffix from bare token, then fullmatch.
    bare_token = re.sub(r"-\d+$", "", token)
    bare = re.fullmatch(
        r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]"
        r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"
        r"|[A-NR-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4}",
        bare_token,
    )
    return bare_token if bare else None


def _isoform_of(token: str) -> str | None:
    """If token contains an isoform suffix (P12345-3), return '3'. Else None."""
    m = ISOFORM_SUFFIX_RE.search(token)
    return m.group(2) if m else None


def load_accessions(path: str) -> list[tuple[str, str, str | None]]:
    """
    Returns a list of (input_token, canonical_accession, isoform_or_None)
    in input order, preserving duplicates. The canonical accession is what
    gets queried against UniProt; the input token is preserved for the
    output TSV so each input line maps to exactly one output row.
    """
    out: list[tuple[str, str, str | None]] = []
    n_iso = 0
    canonicals_with_iso: set[str] = set()
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            acc = parse_accession(line)
            if acc is None:
                print(f"[warn] could not parse: {line!r}", file=sys.stderr)
                continue
            iso = _isoform_of(line)
            if iso is not None:
                n_iso += 1
                canonicals_with_iso.add(acc)
            out.append((line, acc, iso))
    if n_iso:
        print(
            f"[info] {n_iso} isoform-suffixed inputs collapsed to "
            f"{len(canonicals_with_iso)} canonical accessions for feature "
            f"lookup (UniProt features are mapped to canonical sequences)",
            file=sys.stderr,
        )
    return out


# ============================================================
# 2. Fetch features from UniProt (JSON, cached)
# ============================================================

def chunked(seq: list, n: int) -> Iterable[list]:
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def fetch_features(
    accessions: list[str],
    cache_path: Path,
    batch_size: int = BATCH_SIZE,
    retries: int = 3,
    sleep: float = 0.5,
) -> tuple[dict[str, dict], dict[str, str]]:
    """
    Returns ({accession: uniprot_entry_dict}, release_info_dict).
    Cached on disk as a JSON file. release_info_dict carries the UniProt
    release tag and release date scraped from response headers.
    """
    if cache_path.exists():
        with open(cache_path) as fh:
            cache = json.load(fh)
    else:
        cache = {}

    # release of the cached entries, recorded when they were fetched
    release_path = cache_path.with_name("uniprot_release.json")
    release_info: dict[str, str] = {}
    if release_path.exists():
        release_info = json.loads(release_path.read_text())
    missing = [a for a in accessions if a not in cache]
    if not missing:
        print(f"[info] all {len(accessions)} entries cached", file=sys.stderr)
        return cache, release_info

    print(f"[info] fetching {len(missing)} entries from UniProt", file=sys.stderr)
    n_batches = (len(missing) + batch_size - 1) // batch_size

    for i, batch in enumerate(chunked(missing, batch_size), 1):
        query = " OR ".join(f"accession:{a}" for a in batch)
        params = {
            "query": query,
            "format": "json",
            "fields": FIELDS,
            "size": batch_size,
        }
        r = None
        for attempt in range(1, retries + 1):
            try:
                r = requests.get(UNIPROT_SEARCH, params=params, timeout=120)
                r.raise_for_status()
                data = r.json()
                break
            except requests.RequestException as e:
                print(f"[warn] batch {i} attempt {attempt}/{retries}: {e}",
                      file=sys.stderr)
                if attempt == retries:
                    raise
                time.sleep(sleep * attempt)

        # Pull release info from this batch's headers (overwrite each time;
        # all batches in one run should agree).
        if r is not None:
            hdr = _release_from_headers(r.headers)
            if hdr:
                release_info = hdr
                release_path.write_text(json.dumps(release_info))

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            if acc:
                cache[acc] = entry

        # Persist incrementally
        with open(cache_path, "w") as fh:
            json.dump(cache, fh)

        time.sleep(sleep)
        print(f"[info] UniProt batch {i}/{n_batches} done", file=sys.stderr)

    still_missing = [a for a in accessions if a not in cache]
    if still_missing:
        sample = still_missing[:10]
        ellipsis = "..." if len(still_missing) > 10 else ""
        print(f"[warn] {len(still_missing)} accessions not returned: "
              f"{sample}{ellipsis}", file=sys.stderr)

    return cache, release_info


def _release_from_headers(headers) -> dict[str, str]:
    """
    Extract UniProt release info from REST response headers. Headers used
    historically: 'X-UniProt-Release', 'X-UniProt-Release-Date',
    'X-Total-Results'. Names are case-insensitive.
    """
    info: dict[str, str] = {}
    for src_key, out_key in [
        ("X-UniProt-Release",       "uniprot_release"),
        ("X-UniProt-Release-Date",  "uniprot_release_date"),
        ("X-API-Deployment-Date",   "uniprot_api_deployment_date"),
    ]:
        val = headers.get(src_key)
        if val:
            info[out_key] = val.strip()
    return info


# ============================================================
# 3. Feature extraction
# ============================================================

def _loc(feat: dict) -> tuple[int, int]:
    """Return (start, end) as 1-based inclusive positions; (-1,-1) on failure."""
    loc = feat.get("location", {})
    try:
        return (int(loc["start"]["value"]), int(loc["end"]["value"]))
    except (KeyError, TypeError, ValueError):
        return (-1, -1)


def extract_features(entry: dict) -> dict:
    """
    Pull TM, SP, intramembrane, and topological-domain features
    out of a UniProt entry dict.
    """
    tms: list[dict] = []        # alpha-helical or beta-strand TM
    sps: list[dict] = []
    intramems: list[dict] = []
    topo_doms: list[dict] = []  # for N-term orientation

    for feat in entry.get("features", []) or []:
        ftype = feat.get("type", "")
        start, end = _loc(feat)
        desc = feat.get("description", "") or ""
        if ftype == "Transmembrane":
            kind = "beta" if "Beta" in desc else "alpha"  # "Helical" -> alpha
            tms.append({"start": start, "end": end, "kind": kind, "desc": desc})
        elif ftype == "Signal":
            sps.append({"start": start, "end": end})
        elif ftype == "Intramembrane":
            intramems.append({"start": start, "end": end, "desc": desc})
        elif ftype == "Topological domain":
            topo_doms.append({"start": start, "end": end, "desc": desc})

    return {
        "length": entry.get("sequence", {}).get("length")
                  or int(entry.get("length") or 0)
                  or _max_end(tms + sps + topo_doms),
        "tms": sorted(tms, key=lambda x: x["start"]),
        "sps": sps,
        "intramems": intramems,
        "topo_doms": sorted(topo_doms, key=lambda x: x["start"]),
    }


def _max_end(feats: list[dict]) -> int:
    return max((f["end"] for f in feats if f["end"] > 0), default=0)


# ============================================================
# 4. MitoCarta 3.0
# ============================================================

def load_mitocarta_accessions(path: Path) -> set[str]:
    if not path.exists():
        sys.exit(f"MitoCarta file not found: {path}")
    try:
        df = pd.read_excel(path, sheet_name="A Human MitoCarta3.0")
    except Exception:
        df = pd.read_excel(path, sheet_name=0)

    upcol = next((c for c in df.columns if c.strip().lower() == "uniprot"), None)
    if upcol is None:
        upcol = next((c for c in df.columns if "uniprot" in c.lower()), None)
    if upcol is None:
        raise RuntimeError(
            "No UniProt column found in MitoCarta sheet. Columns: "
            + ", ".join(map(str, df.columns))
        )

    accs: set[str] = set()
    for val in df[upcol].dropna():
        for token in re.split(r"[;,|/\s]+", str(val)):
            acc = parse_accession(token)
            if acc:
                accs.add(acc)
    print(f"[info] {len(accs)} UniProt accessions loaded from MitoCarta",
          file=sys.stderr)
    return accs


# ============================================================
# 5. Classifier
# ============================================================

# Map UniProt TOPO_DOM descriptions to a canonical side label.
# Cytoplasmic / matrix / mitochondrial intermembrane in/out — we collapse to "in" (cyto-side)
# vs "out" (lumen/extracellular/periplasmic).
_IN_TERMS = {"Cytoplasmic", "Mitochondrial matrix", "Nuclear",
             "Peroxisomal matrix", "Chloroplast stroma", "Stromal"}
_OUT_TERMS = {"Extracellular", "Lumenal", "Lumen",
              "Mitochondrial intermembrane", "Periplasmic",
              "Chloroplast thylakoid lumen", "Vesicular"}


def _topo_side(desc: str) -> str | None:
    if not desc:
        return None
    d = desc.strip()
    for t in _IN_TERMS:
        if t.lower() in d.lower():
            return "in"
    for t in _OUT_TERMS:
        if t.lower() in d.lower():
            return "out"
    return None


def classify(feat: dict, ta_c_window: int = 40) -> dict:
    L = feat["length"] or 0
    tms = feat["tms"]
    sps = feat["sps"]
    topo_doms = feat["topo_doms"]

    alpha_tms = [t for t in tms if t["kind"] == "alpha"]
    beta_tms = [t for t in tms if t["kind"] == "beta"]
    n_tm = len(alpha_tms)
    n_beta = len(beta_tms)
    has_sp = len(sps) > 0
    notes: list[str] = []

    if n_beta > 0:
        cls = "beta_barrel"
    elif n_tm == 0:
        cls = "secreted_or_lumenal" if has_sp else "soluble"
    elif has_sp:
        cls = "Type_I_or_SP_single_TM" if n_tm == 1 else "SP_plus_multipass"
    else:
        if n_tm == 1:
            tm = alpha_tms[0]
            c_term_after = (L - tm["end"]) if L else -1
            n_term_side = None
            # First topological domain strictly before the TM (closest preceding)
            preceding = [td for td in topo_doms if td["end"] < tm["start"]]
            if preceding:
                n_term_side = _topo_side(preceding[-1]["desc"])

            if 0 <= c_term_after <= ta_c_window:
                cls = "tail_anchored"
                notes.append(f"TM ends {c_term_after} aa from C-term")
            elif n_term_side == "in":
                cls = "signal_anchored_TypeII"
                notes.append("N-term cytoplasmic")
            elif n_term_side == "out":
                cls = "TypeIII_or_unusual"
                notes.append("N-term lumenal but no SP annotated")
            else:
                cls = "single_TM_unclear"
                notes.append("no flanking topological-domain annotation")
        else:
            cls = "multipass_no_SP"

    return {
        "class": cls,
        "n_tm_alpha": n_tm,
        "n_tm_beta": n_beta,
        "has_sp": has_sp,
        "length": L,
        "tm_starts": ",".join(str(t["start"]) for t in alpha_tms),
        "tm_ends":   ",".join(str(t["end"])   for t in alpha_tms),
        "notes": "; ".join(notes),
    }


# ============================================================
# 6. Driver
# ============================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--mitocarta", required=True)
    ap.add_argument("--output", default="topology.tsv")
    ap.add_argument("--workdir", default="topology_work")
    ap.add_argument("--ta-window", type=int, default=40)
    args = ap.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    parsed = load_accessions(args.input)
    print(f"[info] {len(parsed)} input rows parsed", file=sys.stderr)
    if not parsed:
        sys.exit("No accessions parsed.")

    # Unique canonicals for the actual UniProt fetch (preserves order)
    unique_canonicals = list(dict.fromkeys(c for _, c, _ in parsed))
    print(f"[info] {len(unique_canonicals)} unique canonical accessions to query",
          file=sys.stderr)

    cache_path = workdir / "uniprot_features.json"
    entries, release_info = fetch_features(unique_canonicals, cache_path)

    mito = load_mitocarta_accessions(Path(args.mitocarta))

    rows = []
    for input_token, acc, iso in parsed:
        entry = entries.get(acc)
        if entry is None:
            rows.append({
                "input_id": input_token,
                "accession": acc,
                "isoform": iso or "",
                "class": "no_uniprot_entry",
                "n_tm_alpha": pd.NA, "n_tm_beta": pd.NA, "has_sp": pd.NA,
                "length": pd.NA, "tm_starts": "", "tm_ends": "",
                "in_mitocarta3": acc in mito, "notes": "missing from UniProt",
            })
            continue
        feat = extract_features(entry)
        cls = classify(feat, ta_c_window=args.ta_window)
        rows.append({
            "input_id": input_token,
            "accession": acc,
            "isoform": iso or "",
            **cls,
            "in_mitocarta3": acc in mito,
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[info] wrote {len(df)} rows to {args.output}", file=sys.stderr)

    # Metadata file alongside the TSV
    out_path = Path(args.output)
    meta_path = out_path.with_suffix(out_path.suffix + ".meta.txt")
    n_iso = sum(1 for _, _, iso in parsed if iso)
    write_metadata(
        meta_path,
        release_info=release_info,
        n_input=len(parsed),
        n_unique_canonical=len(unique_canonicals),
        n_isoform_inputs=n_iso,
        n_returned=int(df["class"].ne("no_uniprot_entry").sum()),
        n_mito_hits=int(df["in_mitocarta3"].sum()),
        args=args,
        mito_count=len(mito),
    )
    print(f"[info] wrote run metadata to {meta_path}", file=sys.stderr)

    print("\n--- Class summary ---", file=sys.stderr)
    print(df["class"].value_counts().to_string(), file=sys.stderr)
    print(f"\nIn MitoCarta 3.0: {int(df['in_mitocarta3'].sum())}/{len(df)}",
          file=sys.stderr)


def write_metadata(path: Path, release_info: dict, n_input: int,
                   n_unique_canonical: int, n_isoform_inputs: int,
                   n_returned: int, n_mito_hits: int, args, mito_count: int) -> None:
    """
    Write a human- and machine-readable metadata file capturing the
    UniProt release version and run parameters for reproducibility.
    """
    import datetime
    import hashlib
    import platform

    script_path = Path(__file__).resolve()
    try:
        with open(script_path, "rb") as fh:
            script_sha = hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        script_sha = "unknown"

    lines = [
        "# uniprot_topology_pipeline run metadata",
        f"run_timestamp_utc        = {datetime.datetime.utcnow().isoformat()}Z",
        f"uniprot_release          = {release_info.get('uniprot_release', 'unknown')}",
        f"uniprot_release_date     = {release_info.get('uniprot_release_date', 'unknown')}",
        f"uniprot_api_deployment   = {release_info.get('uniprot_api_deployment_date', 'unknown')}",
        f"uniprot_endpoint         = {UNIPROT_SEARCH}",
        f"uniprot_fields           = {FIELDS}",
        f"n_input_rows             = {n_input}",
        f"n_unique_canonical       = {n_unique_canonical}",
        f"n_isoform_suffixed_input = {n_isoform_inputs}",
        f"n_canonicals_returned    = {n_returned}",
        f"n_in_mitocarta3          = {n_mito_hits}",
        f"mitocarta_path           = {Path(args.mitocarta).resolve()}",
        f"mitocarta_n_proteins     = {mito_count}",
        f"ta_window_aa             = {args.ta_window}",
        f"input_path               = {Path(args.input).resolve()}",
        f"output_path              = {Path(args.output).resolve()}",
        f"workdir                  = {Path(args.workdir).resolve()}",
        f"script                   = {script_path}",
        f"script_sha256            = {script_sha}",
        f"python_version           = {platform.python_version()}",
        f"platform                 = {platform.platform()}",
        "",
        "# Note: UniProt TRANSMEM / SIGNAL / TOPO_DOM features are mapped onto",
        "# canonical isoform sequences. Inputs with isoform suffixes (e.g.",
        "# P12345-3) were collapsed to their canonical accession for feature",
        "# lookup. The 'isoform' column in the output TSV records the original",
        "# suffix when present; topology classification reflects canonical features.",
    ]
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
