#!/usr/bin/env python3
"""Select one ITS representative per taxon, align with MAFFT, optional IQ-TREE."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path

import pandas as pd

DATA_ROOT = os.environ.get("FUNGIFLOW_DATA", "/data")
PHASE2_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase2")
PHASE3_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase3")
CONSENSUS_DIR = os.path.join(DATA_ROOT, "consensus_results")
OUT_DIR = os.path.join(DATA_ROOT, "intermediate_data", "phase5")
THREADS = os.environ.get("MAFFT_THREADS", os.environ.get("BLAST_THREADS", "4"))
IQTREE_THREADS = os.environ.get("IQTREE_THREADS", "AUTO")
RUN_TRIMAL = os.environ.get("RUN_TRIMAL", "1") == "1"
MIN_SEQS_IQTREE = 3


def slugify(text: str) -> str:
    s = re.sub(r"[^\w]+", "_", str(text).strip())
    return re.sub(r"_+", "_", s).strip("_")[:80]


def display_tip_label(genus: str, species: str) -> str:
    """Human-readable tip label (species / binomial, no barcode)."""
    g = str(genus or "").strip().replace("_", " ")
    s = str(species or "").strip().replace("_", " ")
    if s and s.lower() not in ("nan", ""):
        g_norm = g.lower()
        s_norm = s.lower()
        if g_norm and (s_norm.startswith(g_norm) or g_norm in s_norm):
            return s
        if g:
            return f"{g} — {s}"
        return s
    return g if g else "Unknown"


def tip_seq_id(genus: str, species: str, rank: int, used: set[str]) -> str:
    """Unique FASTA / IQ-TREE tip name (no spaces)."""
    g = slugify(genus) or "Unknown"
    s = slugify(species) if species and str(species).strip().lower() != "nan" else ""
    base = f"{g}_{s}" if s and s.lower() != g.lower() else g
    candidate = base[:80]
    if candidate in used:
        candidate = f"{candidate}_{rank}"
    used.add(candidate)
    return candidate


def parse_organism(label: str) -> tuple[str, str]:
    if " — " in label:
        genus, species = label.split(" — ", 1)
        return genus.strip(), species.strip()
    return label.strip(), ""


def organism_label(genus: str, species: str) -> str:
    g = str(genus).strip() if pd.notna(genus) else "Unknown"
    s = str(species).strip() if pd.notna(species) else ""
    if s and s.lower() != "nan":
        return f"{g} — {s}"
    return g


def organism_rows(qc: pd.DataFrame, genus: str, species: str) -> pd.DataFrame:
    g = qc["Genus"].fillna("").astype(str) == genus
    s = qc["Species"].fillna("").astype(str) == species
    return qc[g & s].copy()


def load_qc(tag: str) -> pd.DataFrame:
    path = os.path.join(PHASE2_DIR, f"cluster_qc_{tag}.tsv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing {path}. Run cross_sample_qc.py first.")
    df = pd.read_csv(path, sep="\t")
    for col in ("Percent_Identity", "Query_Coverage(%)", "Top2_Pident_Gap", "Query_Length"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_organism_list(scope: str, qc: pd.DataFrame, top_n: int) -> list[str]:
    if scope == "top10":
        path = os.path.join(PHASE3_DIR, "top10_species_full.tsv")
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {path}. Run compile_phase3.py first.")
        top = pd.read_csv(path, sep="\t").head(top_n)
        return top["Organism"].astype(str).tolist()

    if scope == "top10_single_species":
        path = os.path.join(PHASE3_DIR, "top10_single_species_full.tsv")
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {path}. Run compile_phase3.py first.")
        top = pd.read_csv(path, sep="\t").head(top_n)
        return top["Organism"].astype(str).tolist()

    if scope == "all_single_species":
        sub = qc[qc["Single_Organism"].astype(str) == "single_species"].copy()
        sub["Organism"] = sub.apply(
            lambda r: organism_label(r["Genus"], r["Species"]), axis=1
        )
        order = (
            sub.groupby("Organism", sort=False)
            .size()
            .sort_values(ascending=False)
            .head(top_n if top_n > 0 else None)
        )
        return order.index.astype(str).tolist()

    raise ValueError(f"Unknown scope: {scope}")


def pick_best_per_organism(qc: pd.DataFrame, organisms: list[str]) -> pd.DataFrame:
    rows = []
    used_ids: set[str] = set()
    for rank, organism in enumerate(organisms, start=1):
        genus, species = parse_organism(organism)
        sub = organism_rows(qc, genus, species)
        if sub.empty:
            continue
        sub = sub.sort_values(
            ["Percent_Identity", "Query_Coverage(%)", "Top2_Pident_Gap", "Query_Length"],
            ascending=[False, False, False, False],
            na_position="last",
        )
        best = sub.iloc[0].copy()
        best["Organism"] = organism
        best["Representative_Rank"] = rank
        best["Source_Cluster"] = f"{best['Barcode']}_{best['Cluster_Name']}"
        g = str(best.get("Genus", genus) or genus).strip()
        s = str(best.get("Species", species) or species).strip()
        best["Tip_Label"] = display_tip_label(g, s)
        best["Seq_ID"] = tip_seq_id(g, s, rank, used_ids)
        rows.append(best)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def write_representative_fasta(reps: pd.DataFrame, out_path: str) -> int:
    n = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for _, row in reps.iterrows():
            barcode = str(row["Barcode"])
            cluster = str(row["Cluster_Name"])
            seq_id = str(row["Seq_ID"])
            fasta_path = os.path.join(CONSENSUS_DIR, barcode, f"{cluster}.fasta")
            if not os.path.isfile(fasta_path):
                continue
            with open(fasta_path, encoding="utf-8") as fh:
                lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith(">")]
            if not lines:
                continue
            out.write(f">{seq_id}\n")
            out.write("".join(lines) + "\n")
            n += 1
    return n


def run_mafft(n_seqs: int, input_path: str, output_path: str) -> None:
    if n_seqs <= 200:
        cmd = ["mafft", "--thread", str(THREADS), "--auto", input_path]
    elif n_seqs <= 1000:
        cmd = [
            "mafft",
            "--thread",
            str(THREADS),
            "--retree",
            "1",
            "--maxiterate",
            "0",
            input_path,
        ]
    else:
        cmd = [
            "mafft",
            "--thread",
            str(THREADS),
            "--parttree",
            "--retree",
            "1",
            "--maxiterate",
            "0",
            input_path,
        ]
    with open(output_path, "w", encoding="utf-8") as out:
        subprocess.run(cmd, check=True, stdout=out)


def run_trimal(input_path: str, output_path: str) -> bool:
    if not RUN_TRIMAL:
        return False
    try:
        subprocess.run(
            ["trimal", "-in", input_path, "-out", output_path, "-gappyout"],
            check=True,
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def write_tip_meta(reps: pd.DataFrame, csv_path: str) -> None:
    rows = []
    for _, row in reps.iterrows():
        seq_id = str(row["Seq_ID"])
        genus = str(row.get("Genus", "") or "Unknown").strip()
        species = str(row.get("Species", "") or "").strip()
        tip_label = str(row.get("Tip_Label", "") or display_tip_label(genus, species))
        rows.append(
            {
                "cluster": seq_id,
                "genus": genus or "Unknown",
                "species": species,
                "tip_label": tip_label,
                "source_cluster": str(row.get("Source_Cluster", "")),
            }
        )
    pd.DataFrame(rows).to_csv(csv_path, index=False)


def run_iqtree(alignment: str, prefix: str, bootstrap: int, alrt: int) -> None:
    iqtree = "iqtree2" if subprocess.run(["which", "iqtree2"], capture_output=True).returncode == 0 else "iqtree"
    cmd = [
        iqtree,
        "-s",
        alignment,
        "-T",
        str(IQTREE_THREADS),
        "-pre",
        prefix,
        "-redo",
        "-m",
        "MFP",
    ]
    if bootstrap > 0:
        cmd.extend(["-bb", str(bootstrap), "-bnni"])
    if alrt > 0:
        cmd.extend(["-alrt", str(alrt)])
    log_path = f"{prefix}.log"
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scope",
        choices=("top10", "top10_single_species", "all_single_species"),
        default="top10",
        help="Which taxa to include (default: Phase 3 top-10 species, all calls)",
    )
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument(
        "--qc-tag",
        choices=("full", "strict"),
        default="full",
        help="Phase 2 cluster QC table to rank representatives from",
    )
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--run-iqtree", action="store_true", help="Run IQ-TREE after MAFFT")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--alrt", type=int, default=1000)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    qc = load_qc(args.qc_tag)
    organisms = load_organism_list(args.scope, qc, args.top_n)
    reps = pick_best_per_organism(qc, organisms)

    if reps.empty:
        raise SystemExit("No representatives selected. Check Phase 2/3 inputs.")

    reps_path = os.path.join(args.out_dir, "representatives.tsv")
    fasta_path = os.path.join(args.out_dir, "representatives.fasta")
    mafft_path = os.path.join(args.out_dir, "representatives_mafft.fasta")
    trim_path = os.path.join(args.out_dir, "representatives_mafft_trimmed.fasta")
    summary_path = os.path.join(args.out_dir, "phase5_summary.txt")

    reps.to_csv(reps_path, sep="\t", index=False)
    n_fasta = write_representative_fasta(reps, fasta_path)
    if n_fasta < 2:
        raise SystemExit(f"Need ≥2 representative FASTAs; wrote {n_fasta}.")

    print(f"Selected {n_fasta} representatives → {fasta_path}")
    run_mafft(n_fasta, fasta_path, mafft_path)
    print(f"MAFFT → {mafft_path}")

    align_for_tree = mafft_path
    if run_trimal(mafft_path, trim_path):
        align_for_tree = trim_path
        print(f"trimAl → {trim_path}")

    tip_meta_path = f"{align_for_tree}.tip_meta.csv"
    write_tip_meta(reps, tip_meta_path)
    print(f"Tip metadata → {tip_meta_path}")

    treefile = ""
    if args.run_iqtree:
        if n_fasta < MIN_SEQS_IQTREE:
            print(f"Skip IQ-TREE: need ≥{MIN_SEQS_IQTREE} sequences (have {n_fasta}).")
        else:
            iq_prefix = align_for_tree
            print("Running IQ-TREE (MFP + UFBoot + SH-aLRT)…")
            run_iqtree(align_for_tree, iq_prefix, args.bootstrap, args.alrt)
            treefile = f"{iq_prefix}.treefile"
            print(f"IQ-TREE → {treefile}")

    lines = [
        "FungiFlow — representative sequences for phylogeny",
        f"Scope: {args.scope} (n={len(reps)})",
        f"QC source: cluster_qc_{args.qc_tag}.tsv",
        "",
        "Selection rule per taxon: highest Percent_Identity, then Query_Coverage(%),",
        "then Top2_Pident_Gap, then Query_Length.",
        "",
        "Representatives:",
    ]
    show_cols = [
        "Representative_Rank",
        "Organism",
        "Tip_Label",
        "Seq_ID",
        "Source_Cluster",
        "Percent_Identity",
        "Query_Coverage(%)",
        "Confidence",
        "Single_Organism",
    ]
    show_cols = [c for c in show_cols if c in reps.columns]
    for _, row in reps.iterrows():
        pid = row.get("Percent_Identity", "")
        qcov = row.get("Query_Coverage(%)", "")
        lines.append(
            f"  {int(row['Representative_Rank'])}. {row.get('Tip_Label', row['Organism'])} "
            f"← {row.get('Source_Cluster', row.get('Seq_ID', ''))} "
            f"(id {pid}%, cov {qcov}%)"
        )
    lines.extend(
        [
            "",
            "Output files:",
            f"  {reps_path}",
            f"  {fasta_path}",
            f"  {mafft_path}",
        ]
    )
    if os.path.isfile(trim_path):
        lines.append(f"  {trim_path}")
    lines.append(f"  {tip_meta_path}")
    if treefile:
        lines.append(f"  {treefile}")
        lines.append("")
        lines.append("View tree: dashboard → Build IQ-TREE → Tree explorer")

    summary = "\n".join(lines) + "\n"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(summary)
    print(summary)


if __name__ == "__main__":
    main()
