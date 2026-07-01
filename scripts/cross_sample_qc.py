#!/usr/bin/env python3
"""Single-organism QC per ITS cluster and cross-barcode comparison."""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import pandas as pd

DATA_ROOT = os.environ.get("FUNGIFLOW_DATA", "/data")
BLAST_DIR = os.path.join(DATA_ROOT, "blast_results")
OUT_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase2")

NUMERIC_COLS = (
    "Percent_Identity",
    "Alignment_Length",
    "Query_Coverage(%)",
    "Query_Length",
    "Top2_Pident_Gap",
)


def parse_taxonomy(tax_string) -> dict[str, str]:
    out: dict[str, str] = {}
    if tax_string is None or (isinstance(tax_string, float) and pd.isna(tax_string)):
        return out
    raw = str(tax_string).strip()
    if not raw or raw.upper() == "NA":
        return out
    parts = raw.split("|")
    tax_segment = parts[-1] if len(parts) > 1 else raw
    for rank in tax_segment.split(";"):
        if "__" not in rank:
            continue
        lvl, val = rank.split("__", 1)
        val = val.strip()
        if not val:
            continue
        if lvl == "g":
            out["Genus"] = val
        elif lvl == "s":
            out["Species"] = val.replace("_", " ")
    if "Genus" not in out and parts and parts[0] not in ("", "NA"):
        token = parts[0].strip()
        if "_" in token:
            out["Genus"] = token.split("_")[0]
            rest = token.split("_", 1)[1]
            if rest and "Species" not in out:
                out["Species"] = rest.replace("_", " ")
        elif token:
            out["Genus"] = token
    return out


def sample_from_path(path: str) -> str:
    name = Path(path).name
    return name.replace("_blast_summary_strict.tsv", "").replace("_blast_summary.tsv", "")


def load_all_blast_tables(blast_dir: str, strict: bool) -> pd.DataFrame:
    pattern = (
        "*_blast_summary_strict.tsv" if strict else "*_blast_summary.tsv"
    )
    paths = sorted(glob.glob(os.path.join(blast_dir, pattern)))
    if not paths:
        raise FileNotFoundError(f"No BLAST tables matching {pattern} in {blast_dir}")

    frames: list[pd.DataFrame] = []
    for path in paths:
        df = pd.read_csv(path, sep="\t")
        df["Barcode"] = sample_from_path(path)
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    for col in NUMERIC_COLS:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    parsed = merged["Species_Name"].apply(parse_taxonomy).apply(pd.Series)
    for col in ("Genus", "Species"):
        if col in parsed.columns:
            merged[col] = parsed[col]
        else:
            merged[col] = pd.NA

    merged["Genus"] = merged["Genus"].fillna("Unknown")
    merged["Species"] = merged["Species"].fillna("")
    merged["Taxon_Key"] = merged.apply(
        lambda r: (
            f"{r['Genus']}|{r['Species']}"
            if str(r["Species"]).strip()
            else f"{r['Genus']}|"
        ),
        axis=1,
    )
    return merged


def classify_single_organism(row: pd.Series) -> tuple[str, str]:
    """Return (single_organism_call, notes)."""
    conf = str(row.get("Confidence", "")).strip()
    level = str(row.get("Assigned_Level", "")).strip()
    gap = row.get("Top2_Pident_Gap")
    gap_ok = pd.isna(gap) or float(gap) >= 0.5

    if conf == "review_long":
        return "manual_review", "Query >900 bp — check manually"
    if conf == "fail":
        return "no_call", "Below identity/coverage thresholds"
    if conf == "ambiguous":
        return "not_single", "Top-2 BLAST hits too similar (gap <0.5%)"
    if conf in ("high", "medium", "low_species") and level == "species":
        if gap_ok:
            return "single_species", f"Species-level ({conf})"
        return "not_single", "Species call but top-2 hits too close"
    if conf == "genus_only" or level == "genus":
        return "single_genus", "Genus-level only — species uncertain"
    return "unclear", f"Confidence={conf}, level={level}"


def add_qc_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    classified = out.apply(classify_single_organism, axis=1, result_type="expand")
    out["Single_Organism"] = classified[0]
    out["QC_Notes"] = classified[1]
    out["Is_Single_Organism"] = out["Single_Organism"].isin(
        ("single_species", "single_genus")
    )
    return out


def build_per_barcode_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for barcode, sub in df.groupby("Barcode", sort=True):
        row = {
            "Barcode": barcode,
            "Clusters": len(sub),
            "Single_species": int((sub["Single_Organism"] == "single_species").sum()),
            "Single_genus": int((sub["Single_Organism"] == "single_genus").sum()),
            "Not_single": int((sub["Single_Organism"] == "not_single").sum()),
            "No_call": int((sub["Single_Organism"] == "no_call").sum()),
            "Manual_review": int((sub["Single_Organism"] == "manual_review").sum()),
            "Unique_genera": int(sub["Genus"].nunique()),
            "Unique_species_keys": int(
                sub.loc[sub["Species"].astype(str).str.strip() != "", "Taxon_Key"].nunique()
            ),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def build_genus_cross_sample(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for genus, sub in df.groupby("Genus", sort=True):
        barcodes = sorted(sub["Barcode"].unique().tolist())
        species_in_genus = sorted(
            {s for s in sub["Species"].astype(str).tolist() if s.strip()}
        )
        rows.append(
            {
                "Genus": genus,
                "Clusters_total": len(sub),
                "Barcodes_n": len(barcodes),
                "Barcodes": ", ".join(barcodes),
                "Species_assigned_n": len(species_in_genus),
                "Species_list": "; ".join(species_in_genus[:8])
                + (" …" if len(species_in_genus) > 8 else ""),
                "Single_species_clusters": int(
                    (sub["Single_Organism"] == "single_species").sum()
                ),
                "In_all_barcodes": len(barcodes) == df["Barcode"].nunique(),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["Barcodes_n", "Clusters_total"], ascending=[False, False]
    ).reset_index(drop=True)


def build_species_conflicts(df: pd.DataFrame) -> pd.DataFrame:
    """Genera where different species names appear across barcodes."""
    rows = []
    species_df = df[df["Species"].astype(str).str.strip() != ""].copy()
    for genus, sub in species_df.groupby("Genus"):
        species_by_barcode: dict[str, set[str]] = {}
        for barcode, bsub in sub.groupby("Barcode"):
            species_by_barcode[barcode] = set(bsub["Species"].astype(str).str.strip())
        all_species = set.union(*species_by_barcode.values()) if species_by_barcode else set()
        if len(all_species) <= 1:
            continue
        rows.append(
            {
                "Genus": genus,
                "Distinct_species": len(all_species),
                "Species": "; ".join(sorted(all_species)),
                "Barcodes": ", ".join(sorted(species_by_barcode.keys())),
                "Detail": " | ".join(
                    f"{b}: {', '.join(sorted(s))}"
                    for b, s in sorted(species_by_barcode.items())
                ),
            }
        )
    return pd.DataFrame(rows)


def write_report(
    df: pd.DataFrame,
    out_dir: str,
    strict: bool,
) -> None:
    os.makedirs(out_dir, exist_ok=True)
    tag = "strict" if strict else "full"

    merged_path = os.path.join(out_dir, f"merged_blast_{tag}.tsv")
    qc_path = os.path.join(out_dir, f"cluster_qc_{tag}.tsv")
    barcode_path = os.path.join(out_dir, f"per_barcode_qc_{tag}.tsv")
    genus_path = os.path.join(out_dir, f"genus_cross_sample_{tag}.tsv")
    conflict_path = os.path.join(out_dir, f"species_conflicts_{tag}.tsv")
    summary_path = os.path.join(out_dir, f"phase2_summary_{tag}.txt")

    qc_cols = [
        "Barcode",
        "Cluster_Name",
        "Query_Length",
        "Percent_Identity",
        "Query_Coverage(%)",
        "Top2_Pident_Gap",
        "Genus",
        "Species",
        "Confidence",
        "Assigned_Level",
        "Single_Organism",
        "Is_Single_Organism",
        "QC_Notes",
        "Species_Name",
    ]
    qc_cols = [c for c in qc_cols if c in df.columns]

    df.to_csv(merged_path, sep="\t", index=False)
    df[qc_cols].to_csv(qc_path, sep="\t", index=False)

    per_barcode = build_per_barcode_summary(df)
    genus_cross = build_genus_cross_sample(df)
    conflicts = build_species_conflicts(df)

    per_barcode.to_csv(barcode_path, sep="\t", index=False)
    genus_cross.to_csv(genus_path, sep="\t", index=False)
    conflicts.to_csv(conflict_path, sep="\t", index=False)

    n_barcodes = df["Barcode"].nunique()
    n_clusters = len(df)
    union_genera = df["Genus"].nunique()
    shared_genera = genus_cross.loc[genus_cross["In_all_barcodes"], "Genus"].tolist()
    single_species = int((df["Single_Organism"] == "single_species").sum())
    single_genus = int((df["Single_Organism"] == "single_genus").sum())
    not_single = int((df["Single_Organism"] == "not_single").sum())
    no_call = int((df["Single_Organism"] == "no_call").sum())

    lines = [
        "FungiFlow Phase 2 — single-organism QC",
        f"BLAST set: {tag}",
        f"Barcodes: {n_barcodes}",
        f"Total clusters: {n_clusters}",
        "",
        "Per-cluster single-organism calls:",
        f"  single_species (confident species): {single_species}",
        f"  single_genus (genus only):          {single_genus}",
        f"  not_single (ambiguous top-2):         {not_single}",
        f"  no_call (fail):                       {no_call}",
        f"  manual_review (long reads):           {(df['Single_Organism'] == 'manual_review').sum()}",
        "",
        f"Unique genera across all barcodes: {union_genera}",
        f"Genera present in all {n_barcodes} barcodes: {len(shared_genera)}",
    ]
    if shared_genera:
        lines.append("  " + ", ".join(shared_genera[:20]))
        if len(shared_genera) > 20:
            lines.append(f"  … +{len(shared_genera) - 20} more")
    lines.extend(
        [
            "",
            f"Genera with conflicting species across barcodes: {len(conflicts)}",
            "",
            "Output files:",
            f"  {merged_path}",
            f"  {qc_path}",
            f"  {barcode_path}",
            f"  {genus_path}",
            f"  {conflict_path}",
        ]
    )

    summary_text = "\n".join(lines) + "\n"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(summary_text)
    print(summary_text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--blast-dir",
        default=BLAST_DIR,
        help="Directory with *_blast_summary.tsv files",
    )
    parser.add_argument(
        "--out-dir",
        default=OUT_DIR,
        help="Output directory for Phase 2 tables",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Use *_blast_summary_strict.tsv instead of full tables",
    )
    args = parser.parse_args()

    df = load_all_blast_tables(args.blast_dir, strict=args.strict)
    df = add_qc_columns(df)
    write_report(df, args.out_dir, strict=args.strict)


if __name__ == "__main__":
    main()
