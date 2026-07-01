#!/usr/bin/env python3
"""Phase 3: compile all BLAST runs, count unique organisms, select top 10."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

DATA_ROOT = os.environ.get("FUNGIFLOW_DATA", "/data")
PHASE2_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase2")
OUT_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase3")
TOP_N = 10


def load_compiled(tag: str) -> pd.DataFrame:
    path = os.path.join(PHASE2_DIR, f"cluster_qc_{tag}.tsv")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Missing {path}. Run cross_sample_qc.py first (Phase 2)."
        )
    return pd.read_csv(path, sep="\t")


def organism_label(genus: str, species: str) -> str:
    g = str(genus).strip() if pd.notna(genus) else "Unknown"
    s = str(species).strip() if pd.notna(species) else ""
    if s and s.lower() != "nan":
        return f"{g} — {s}"
    return g


def count_unique_organisms(df: pd.DataFrame) -> dict[str, int]:
    all_genera = df["Genus"].fillna("Unknown").astype(str).nunique()
    species_rows = df[df["Species"].fillna("").astype(str).str.strip() != ""]
    all_species_labels = species_rows.apply(
        lambda r: organism_label(r["Genus"], r["Species"]), axis=1
    ).nunique()

    confident = df[df["Single_Organism"].isin(["single_species", "single_genus"])]
    conf_genera = confident["Genus"].fillna("Unknown").astype(str).nunique()

    single_sp = df[df["Single_Organism"] == "single_species"]
    conf_species = single_sp.apply(
        lambda r: organism_label(r["Genus"], r["Species"]), axis=1
    ).nunique()

    return {
        "total_clusters": len(df),
        "barcodes": df["Barcode"].nunique(),
        "unique_genera_all_calls": all_genera,
        "unique_species_labels_all_calls": all_species_labels,
        "unique_genera_single_organism": conf_genera,
        "unique_species_single_species_calls": conf_species,
    }


def build_top_organisms(df: pd.DataFrame, level: str, n: int) -> pd.DataFrame:
    work = df.copy()
    if level == "genus":
        work["Organism"] = work["Genus"].fillna("Unknown").astype(str)
    else:
        work["Organism"] = work.apply(
            lambda r: organism_label(r["Genus"], r["Species"]), axis=1
        )
        work = work[work["Species"].fillna("").astype(str).str.strip() != ""]

    rows = []
    for organism, sub in work.groupby("Organism", sort=False):
        barcodes = sorted(sub["Barcode"].astype(str).unique().tolist())
        rows.append(
            {
                "Rank": 0,
                "Organism": organism,
                "Level": level,
                "Clusters": len(sub),
                "Barcodes_n": len(barcodes),
                "Barcodes": ", ".join(barcodes),
                "Mean_identity_pct": round(float(sub["Percent_Identity"].mean()), 2),
                "Mean_query_cov_pct": round(float(sub["Query_Coverage(%)"].mean()), 2),
                "Single_species_clusters": int(
                    (sub["Single_Organism"] == "single_species").sum()
                ),
                "Single_genus_clusters": int(
                    (sub["Single_Organism"] == "single_genus").sum()
                ),
            }
        )

    out = (
        pd.DataFrame(rows)
        .sort_values(["Clusters", "Barcodes_n"], ascending=[False, False])
        .head(n)
        .reset_index(drop=True)
    )
    out["Rank"] = range(1, len(out) + 1)
    return out


def build_top_single_species(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Top organisms using only confident species-level assignments."""
    sub = df[df["Single_Organism"] == "single_species"].copy()
    return build_top_organisms(sub, "species", n)


def write_phase3(df: pd.DataFrame, tag: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    counts = count_unique_organisms(df)

    top_genus = build_top_organisms(df, "genus", TOP_N)
    top_species = build_top_organisms(df, "species", TOP_N)
    top_confident = build_top_single_species(df, TOP_N)

    counts_path = os.path.join(out_dir, f"unique_organisms_{tag}.tsv")
    top_genus_path = os.path.join(out_dir, f"top10_genera_{tag}.tsv")
    top_species_path = os.path.join(out_dir, f"top10_species_{tag}.tsv")
    top_conf_path = os.path.join(out_dir, f"top10_single_species_{tag}.tsv")
    summary_path = os.path.join(out_dir, f"phase3_summary_{tag}.txt")

    pd.DataFrame([counts]).to_csv(counts_path, sep="\t", index=False)
    top_genus.to_csv(top_genus_path, sep="\t", index=False)
    top_species.to_csv(top_species_path, sep="\t", index=False)
    top_confident.to_csv(top_conf_path, sep="\t", index=False)

    lines = [
        "FungiFlow Phase 3 — compiled organisms across all barcodes",
        f"Source: cluster_qc_{tag}.tsv (Phase 2)",
        "",
        "Unique organism counts:",
        f"  Barcodes:                          {counts['barcodes']}",
        f"  Total clusters:                    {counts['total_clusters']}",
        f"  Unique genera (all BLAST calls):     {counts['unique_genera_all_calls']}",
        f"  Unique species labels (all calls): {counts['unique_species_labels_all_calls']}",
        f"  Unique genera (single-organism QC): {counts['unique_genera_single_organism']}",
        f"  Unique species (single_species QC): {counts['unique_species_single_species_calls']}",
        "",
        f"Top {TOP_N} genera by cluster count:",
    ]
    for _, row in top_genus.iterrows():
        lines.append(
            f"  {int(row['Rank'])}. {row['Organism']} — {int(row['Clusters'])} clusters "
            f"in {int(row['Barcodes_n'])} barcodes"
        )

    lines.extend(["", f"Top {TOP_N} species (all calls with species name):"])
    for _, row in top_species.iterrows():
        lines.append(
            f"  {int(row['Rank'])}. {row['Organism']} — {int(row['Clusters'])} clusters"
        )

    lines.extend(["", f"Top {TOP_N} (single_species QC only — best for downstream work):"])
    for _, row in top_confident.iterrows():
        lines.append(
            f"  {int(row['Rank'])}. {row['Organism']} — {int(row['Clusters'])} clusters "
            f"(mean id {row['Mean_identity_pct']}%, cov {row['Mean_query_cov_pct']}%)"
        )

    lines.extend(
        [
            "",
            "Output files:",
            f"  {counts_path}",
            f"  {top_genus_path}",
            f"  {top_species_path}",
            f"  {top_conf_path}",
        ]
    )

    summary = "\n".join(lines) + "\n"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(summary)
    print(summary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tag",
        choices=("full", "strict"),
        default="full",
        help="Use Phase 2 tables from full or strict BLAST set",
    )
    parser.add_argument("--out-dir", default=OUT_DIR)
    args = parser.parse_args()

    df = load_compiled(args.tag)
    write_phase3(df, args.tag, args.out_dir)


if __name__ == "__main__":
    main()
