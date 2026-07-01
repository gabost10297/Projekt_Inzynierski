#!/usr/bin/env python3
"""Phase 4: MAFFT alignment per top-N species (all BLAST calls)."""

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
OUT_DIR = os.path.join(DATA_ROOT, "intermediate_data", "mafft", "by_organism")
THREADS = os.environ.get("MAFFT_THREADS", os.environ.get("BLAST_THREADS", "4"))
RUN_TRIMAL = os.environ.get("RUN_TRIMAL", "1") == "1"
MIN_SEQS = 2


def slugify(text: str) -> str:
    s = re.sub(r"[^\w]+", "_", text.strip())
    return re.sub(r"_+", "_", s).strip("_")[:120]


def parse_organism(label: str) -> tuple[str, str]:
    if " — " in label:
        genus, species = label.split(" — ", 1)
        return genus.strip(), species.strip()
    return label.strip(), ""


def load_top_organisms(path: str, n: int) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return df.head(n).copy()


def organism_rows(qc: pd.DataFrame, genus: str, species: str) -> pd.DataFrame:
    g = qc["Genus"].fillna("").astype(str) == genus
    s = qc["Species"].fillna("").astype(str) == species
    return qc[g & s].copy()


def write_input_fasta(rows: pd.DataFrame, out_path: str) -> tuple[int, int]:
    n_written = 0
    n_missing = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for _, row in rows.iterrows():
            barcode = str(row["Barcode"])
            cluster = str(row["Cluster_Name"])
            fasta_path = os.path.join(CONSENSUS_DIR, barcode, f"{cluster}.fasta")
            seq_id = f"{barcode}_{cluster}"
            if not os.path.isfile(fasta_path):
                n_missing += 1
                continue
            with open(fasta_path, encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            seq_lines = [ln for ln in lines if ln and not ln.startswith(">")]
            if not seq_lines:
                n_missing += 1
                continue
            out.write(f">{seq_id}\n")
            out.write("".join(seq_lines) + "\n")
            n_written += 1
    return n_written, n_missing


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


def process_organism(
    rank: int,
    organism: str,
    qc: pd.DataFrame,
    out_dir: str,
) -> dict:
    genus, species = parse_organism(organism)
    slug = slugify(f"{genus}_{species}") if species else slugify(genus)
    rows = organism_rows(qc, genus, species)

    input_fasta = os.path.join(out_dir, f".{slug}.input.fasta")
    raw_fasta = os.path.join(out_dir, f"{slug}_mafft.fasta")
    trim_fasta = os.path.join(out_dir, f"{slug}_mafft_trimmed.fasta")
    cluster_list = os.path.join(out_dir, f"{slug}_clusters.tsv")

    result = {
        "Rank": rank,
        "Organism": organism,
        "Genus": genus,
        "Species": species,
        "Clusters_expected": len(rows),
        "Sequences_aligned": 0,
        "Sequences_missing": 0,
        "Status": "pending",
        "Raw_alignment": "",
        "Trimmed_alignment": "",
        "Cluster_list": cluster_list,
    }

    if rows.empty:
        result["Status"] = "no_clusters"
        return result

    rows.to_csv(cluster_list, sep="\t", index=False)
    n_seqs, n_missing = write_input_fasta(rows, input_fasta)
    result["Sequences_aligned"] = n_seqs
    result["Sequences_missing"] = n_missing

    if n_seqs < MIN_SEQS:
        result["Status"] = f"skipped_lt_{MIN_SEQS}_seqs"
        if os.path.isfile(input_fasta):
            os.remove(input_fasta)
        return result

    print(f">>> [{rank}] {organism} — {n_seqs} sequences <<<")
    try:
        run_mafft(n_seqs, input_fasta, raw_fasta)
        result["Raw_alignment"] = raw_fasta
        if run_trimal(raw_fasta, trim_fasta):
            result["Trimmed_alignment"] = trim_fasta
        result["Status"] = "ok"
        print(f"  Saved: {raw_fasta}")
        if result["Trimmed_alignment"]:
            print(f"  Saved: {trim_fasta}")
    except subprocess.CalledProcessError as exc:
        result["Status"] = f"mafft_failed:{exc.returncode}"
    finally:
        if os.path.isfile(input_fasta):
            os.remove(input_fasta)

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--top10",
        default=os.path.join(PHASE3_DIR, "top10_species_full.tsv"),
        help="Top organisms table from Phase 3",
    )
    parser.add_argument(
        "--qc",
        default=os.path.join(PHASE2_DIR, "cluster_qc_full.tsv"),
        help="Cluster QC table from Phase 2",
    )
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    if not os.path.isfile(args.top10):
        raise FileNotFoundError(f"Missing {args.top10}. Run compile_phase3.py first.")
    if not os.path.isfile(args.qc):
        raise FileNotFoundError(f"Missing {args.qc}. Run cross_sample_qc.py first.")

    os.makedirs(args.out_dir, exist_ok=True)
    top = load_top_organisms(args.top10, args.top_n)
    qc = pd.read_csv(args.qc, sep="\t")

    print("Phase 4 — MAFFT per top organism (all calls)")
    print("=" * 60)
    print(f"Threads: {THREADS} | trimAl: {RUN_TRIMAL}")
    print(f"Output: {args.out_dir}")

    manifest_rows = []
    for _, row in top.iterrows():
        manifest_rows.append(
            process_organism(int(row["Rank"]), str(row["Organism"]), qc, args.out_dir)
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(args.out_dir, "phase4_manifest.tsv")
    summary_path = os.path.join(args.out_dir, "phase4_summary.txt")
    manifest.to_csv(manifest_path, sep="\t", index=False)

    lines = [
        "FungiFlow Phase 4 — MAFFT per top-10 species (all calls)",
        f"Source top list: {args.top10}",
        f"Source clusters: {args.qc}",
        "",
    ]
    for _, m in manifest.iterrows():
        lines.append(
            f"  {int(m['Rank'])}. {m['Organism']}: "
            f"{int(m['Sequences_aligned'])} seqs aligned — {m['Status']}"
        )
    lines.extend(["", f"Manifest: {manifest_path}"])
    summary = "\n".join(lines) + "\n"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(summary)

    print("=" * 60)
    print(summary)


if __name__ == "__main__":
    main()
