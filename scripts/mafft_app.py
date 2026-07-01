import glob
import html
import os
import re
import subprocess
from collections import Counter
from io import BytesIO, StringIO

import pandas as pd
import plotly.express as px
import streamlit as st

MAFFT_DIR = "/data/intermediate_data/mafft"
ORGANISM_MAFFT_DIR = os.path.join(MAFFT_DIR, "by_organism")
BLAST_DIR = "/data/blast_results"
CONSENSUS_DIR = "/data/consensus_results"
MAFFT_SCRIPT = "/data/scripts/run_mafft.sh"
ORGANISM_MAFFT_SCRIPT = "/data/scripts/run_phase4_mafft.py"
MAFFT_LOG = f"{MAFFT_DIR}/mafft_run.log"
MAFFT_RUNNING_FLAG = f"{MAFFT_DIR}/mafft_running.flag"
ORGANISM_MAFFT_LOG = f"{ORGANISM_MAFFT_DIR}/phase4_run.log"
ORGANISM_MAFFT_RUNNING_FLAG = f"{ORGANISM_MAFFT_DIR}/phase4_running.flag"
ORGANISM_MANIFEST = f"{ORGANISM_MAFFT_DIR}/phase4_manifest.tsv"
MANIFEST = f"{MAFFT_DIR}/manifest.tsv"

PLOTLY_LAYOUT = dict(
    template="plotly_white",
    font=dict(color="#0a160a"),
    colorway=["#1b4332", "#2d6a4f", "#40916c", "#52b788", "#74c69d", "#95d5b2"],
)

LABEL_EXTRA_FIELDS = [
    "Genus",
    "Species",
    "Confidence",
    "Percent_Identity",
    "Assigned_Level",
]

LABEL_FIELD_LABELS = {
    "Cluster_Name": "Cluster",
    "Barcode": "Barcode",
    "Genus": "Genus",
    "Species": "Species",
    "Confidence": "Confidence",
    "Percent_Identity": "Identity %",
    "Assigned_Level": "Assigned level",
}

SETS_DISPLAY_COLS = [
    "Cluster_Name",
    "Genus",
    "Species",
    "Percent_Identity",
    "Confidence",
    "Assigned_Level",
    "Query_Length",
    "Length_Tier",
]


def page_hero(title: str, subtitle: str) -> None:
    st.markdown(
        f'<div class="ff-page-hero"><h1>{title}</h1><p>{subtitle}</p></div>',
        unsafe_allow_html=True,
    )


def list_blast_samples() -> list[str]:
    paths = glob.glob(os.path.join(BLAST_DIR, "*_blast_summary.tsv"))
    samples = []
    for p in paths:
        name = os.path.basename(p)
        if name.endswith("_blast_summary.tsv") and "_strict" not in name:
            samples.append(name.replace("_blast_summary.tsv", ""))
    return sorted(set(samples))


def blast_tsv_path(sample: str, mode: str) -> str:
    if mode == "strict":
        return os.path.join(BLAST_DIR, f"{sample}_blast_summary_strict.tsv")
    return os.path.join(BLAST_DIR, f"{sample}_blast_summary.tsv")


def alignment_paths(sample: str, mode: str) -> dict[str, str]:
    return {
        "raw": os.path.join(MAFFT_DIR, f"{sample}_{mode}_mafft.fasta"),
        "trimmed": os.path.join(MAFFT_DIR, f"{sample}_{mode}_mafft_trimmed.fasta"),
    }


def selection_key(sample: str, mode: str) -> str:
    return f"mafft_sel_{sample}_{mode}"


def parse_taxonomy_from_species_name(tax_string) -> dict[str, str]:
    """Extract genus/species from BLAST Species_Name (same logic as blast_app)."""
    out: dict[str, str] = {}
    if tax_string is None or (isinstance(tax_string, float) and pd.isna(tax_string)):
        return out
    try:
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
    except Exception:
        pass
    return out


def enrich_blast_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Add Genus/Species from Species_Name (raw BLAST TSV has no separate columns)."""
    if df.empty or "Species_Name" not in df.columns:
        return df

    parsed_rows = df["Species_Name"].apply(parse_taxonomy_from_species_name)
    parsed = pd.DataFrame(list(parsed_rows))

    for col in ("Genus", "Species"):
        if col not in parsed.columns:
            continue
        if col not in df.columns:
            df[col] = parsed[col]
        else:
            empty = df[col].isna() | (df[col].astype(str).str.strip() == "")
            df.loc[empty, col] = parsed.loc[empty, col]

    return df


@st.cache_data(show_spinner=False)
def load_blast_table(tsv_path: str, mtime: float) -> pd.DataFrame:
    del mtime  # cache bust when file changes
    if not os.path.isfile(tsv_path):
        return pd.DataFrame()
    df = pd.read_csv(tsv_path, sep="\t")
    return enrich_blast_dataframe(df)


@st.cache_data(show_spinner=False)
def load_alignment(fasta_path: str, mtime: float) -> dict[str, str]:
    del mtime
    sequences: dict[str, str] = {}
    current = ""
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].strip()
                sequences[current] = ""
            elif current:
                sequences[current] += line.upper()
    return sequences


def default_cluster_selection(df: pd.DataFrame, mode: str) -> list[str]:
    if df.empty or "Cluster_Name" not in df.columns:
        return []
    if "Confidence" in df.columns:
        conf = df["Confidence"].astype(str)
        if mode == "strict":
            pick = df[conf.isin(["high", "medium"])]["Cluster_Name"]
        else:
            pick = df[~conf.isin(["fail"])]["Cluster_Name"]
        ids = pick.astype(str).tolist()
        if ids:
            return ids
    return df["Cluster_Name"].astype(str).head(min(20, len(df))).tolist()


def preflight_stats(sample: str, mode: str) -> dict:
    tsv = blast_tsv_path(sample, mode)
    consensus_dir = os.path.join(CONSENSUS_DIR, sample)
    n_tsv = 0
    n_fasta = 0
    n_missing = 0

    if os.path.isfile(tsv):
        df = load_blast_table(tsv, os.path.getmtime(tsv))
        if not df.empty and "Cluster_Name" in df.columns:
            clusters = df["Cluster_Name"].astype(str).tolist()
            n_tsv = len(clusters)
            for cluster in clusters:
                if os.path.isfile(os.path.join(consensus_dir, f"{cluster}.fasta")):
                    n_fasta += 1
                else:
                    n_missing += 1

    paths = alignment_paths(sample, mode)
    return {
        "n_tsv": n_tsv,
        "n_fasta": n_fasta,
        "n_missing": n_missing,
        "has_raw": os.path.isfile(paths["raw"]),
        "has_trimmed": os.path.isfile(paths["trimmed"]),
    }


def mafft_is_running() -> bool:
    return os.path.isfile(MAFFT_RUNNING_FLAG)


def read_log_tail(n: int = 20) -> str:
    if not os.path.isfile(MAFFT_LOG):
        return ""
    try:
        with open(MAFFT_LOG, encoding="utf-8", errors="replace") as f:
            return "".join(f.readlines()[-n:])
    except OSError:
        return ""


def launch_mafft(args: list[str]) -> None:
    os.makedirs(MAFFT_DIR, exist_ok=True)
    with open(MAFFT_RUNNING_FLAG, "w", encoding="utf-8") as f:
        f.write("running\n")
    shell_cmd = (
        f"bash {MAFFT_SCRIPT} {' '.join(args)} >> {MAFFT_LOG} 2>&1; "
        f"rm -f {MAFFT_RUNNING_FLAG}"
    )
    with open(MAFFT_LOG, "w", encoding="utf-8") as log:
        log.write(f"$ {shell_cmd}\n\n")
    subprocess.Popen(
        ["bash", "-c", shell_cmd],
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    load_blast_table.clear()
    load_alignment.clear()


def load_manifest() -> pd.DataFrame:
    if not os.path.isfile(MANIFEST):
        return pd.DataFrame()
    try:
        return pd.read_csv(MANIFEST, sep="\t")
    except (pd.errors.ParserError, OSError, ValueError):
        return pd.DataFrame()


def alignment_fasta_bytes(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()


def subset_fasta_bytes(sequences: dict[str, str], selected: list[str]) -> bytes:
    buf = StringIO()
    for seq_id in selected:
        if seq_id not in sequences:
            continue
        buf.write(f">{seq_id}\n")
        buf.write(sequences[seq_id])
        buf.write("\n")
    return buf.getvalue().encode("utf-8")


def render_run_tab(sample: str, mode: str) -> None:
    stats = preflight_stats(sample, mode)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Clusters in TSV", stats["n_tsv"])
    c2.metric("FASTAs found", stats["n_fasta"])
    c3.metric("Missing FASTA", stats["n_missing"])
    c4.metric("Ready for MAFFT", "Yes" if stats["n_fasta"] >= 2 else "No")

    if stats["n_fasta"] < 2:
        st.warning("Need at least 2 consensus FASTAs matching the BLAST table.")

    if mafft_is_running():
        st.warning("MAFFT is running in the background.")
        if st.button("Refresh status", key="mafft_refresh_running"):
            st.rerun()
        st.markdown("**Log (last lines)**")
        st.code(read_log_tail(25) or "(empty)", language="log")
        return

    st.caption(
        "Aligns clusters listed in the selected BLAST table. "
        "Uses faster MAFFT settings for large sets (>200 / >1000 sequences)."
    )

    b1, b2, b3 = st.columns(3)
    if b1.button("This sample only", type="primary", key="mafft_run_one"):
        launch_mafft([sample, mode])
        st.rerun()
    if b2.button("All samples (strict + full)", key="mafft_run_all"):
        launch_mafft(["--all", "both"])
        st.rerun()
    if b3.button("All samples (current set only)", key="mafft_run_all_mode"):
        launch_mafft(["--all", mode])
        st.rerun()

    manifest = load_manifest()
    if not manifest.empty:
        st.markdown("**Recent MAFFT runs**")
        st.dataframe(manifest, hide_index=True, width="stretch")

    with st.expander("Full log", expanded=False):
        st.code(read_log_tail(80) or "No log yet.", language="log")


def render_sets_tab(sample: str, mode: str, df: pd.DataFrame) -> list[str]:
    sel_key = selection_key(sample, mode)

    if df.empty:
        st.warning("BLAST table is empty or missing.")
        return []

    if sel_key not in st.session_state:
        st.session_state[sel_key] = default_cluster_selection(df, mode)


    with st.container(border=True):
        st.markdown("**Filter clusters**")
        f1, f2, f3 = st.columns(3)
        view = df.copy()
        if "Genus" in view.columns:
            genera = ["(all)"] + sorted(view["Genus"].dropna().astype(str).unique().tolist())
            genus_pick = f1.selectbox("Genus", genera, key=f"mafft_sets_genus_{sample}_{mode}")
            if genus_pick != "(all)":
                view = view[view["Genus"].astype(str) == genus_pick]
        if "Confidence" in view.columns:
            conf_opts = sorted(view["Confidence"].dropna().astype(str).unique().tolist())
            conf_pick = f2.multiselect(
                "Confidence",
                conf_opts,
                default=conf_opts,
                key=f"mafft_sets_conf_{sample}_{mode}",
            )
            if conf_pick:
                view = view[view["Confidence"].astype(str).isin(conf_pick)]
        if "Percent_Identity" in view.columns:
            min_pid = f3.slider(
                "Min. identity %",
                0.0,
                100.0,
                0.0,
                0.5,
                key=f"mafft_sets_pid_{sample}_{mode}",
            )
            view = view[view["Percent_Identity"] >= min_pid]

        show_cols = [c for c in SETS_DISPLAY_COLS if c in view.columns]
        st.caption(f"{len(view)} / {len(df)} clusters after filters")
        st.dataframe(view[show_cols] if show_cols else view, hide_index=True, width="stretch")

    options = df["Cluster_Name"].astype(str).tolist()
    valid_defaults = [x for x in st.session_state.get(sel_key, []) if x in options]

    c1, c2 = st.columns(2)
    if c1.button("Select high + medium", key=f"mafft_pick_hm_{sample}_{mode}"):
        st.session_state[sel_key] = default_cluster_selection(df, "strict")
        st.rerun()
    if c2.button("Select all in table", key=f"mafft_pick_all_{sample}_{mode}"):
        st.session_state[sel_key] = options
        st.rerun()

    selected = st.multiselect(
        "Clusters to show in alignment viewer",
        options=options,
        default=valid_defaults,
        key=f"mafft_multiselect_{sample}_{mode}",
    )
    st.session_state[sel_key] = selected
    st.info(f"**{len(selected)}** sequences selected for the View tab.")
    return selected




def cluster_metadata_lookup(sample: str, mode: str) -> dict[str, dict]:
    """Map Cluster_Name -> BLAST row for alignment row labels."""
    tsv = blast_tsv_path(sample, mode)
    if not os.path.isfile(tsv):
        return {}
    df = load_blast_table(tsv, os.path.getmtime(tsv))
    if df.empty or "Cluster_Name" not in df.columns:
        return {}
    lookup: dict[str, dict] = {}
    for _, row in df.iterrows():
        cid = str(row["Cluster_Name"]).strip()
        meta = row.to_dict()
        for key, val in parse_taxonomy_from_species_name(meta.get("Species_Name")).items():
            if key not in meta or pd.isna(meta.get(key)) or not str(meta.get(key)).strip():
                meta[key] = val
        lookup[cid] = meta
        if cid.isdigit():
            lookup[str(int(cid))] = meta
    return lookup


def resolve_cluster_meta(seq_id: str, meta_by_cluster: dict[str, dict]) -> dict | None:
    """Match alignment FASTA headers to BLAST Cluster_Name."""
    sid = str(seq_id).strip()
    if sid in meta_by_cluster:
        return meta_by_cluster[sid]
    for key, row in meta_by_cluster.items():
        if str(key).strip() == sid:
            return row
        if sid.endswith(str(key)) or str(key).endswith(sid):
            return row
    m = re.match(r"^(barcode\d+)_(.+)$", sid)
    if m and sid in meta_by_cluster:
        return meta_by_cluster[sid]
    if m:
        cluster = m.group(2)
        if cluster in meta_by_cluster:
            return meta_by_cluster[cluster]
    if sid.isdigit():
        alt = str(int(sid))
        if alt in meta_by_cluster:
            return meta_by_cluster[alt]
    return None


def format_row_label(
    cluster_id: str,
    meta: dict | None,
    fields: list[str],
) -> tuple[str, str]:
    """HTML display (may include <br>) and tooltip title."""
    lines: list[str] = []
    for field in fields:
        if field == "Cluster_Name":
            lines.append(html.escape(str(cluster_id)))
            continue
        if field == "Barcode":
            barcode_val = ""
            if meta and meta.get("Barcode") and not pd.isna(meta.get("Barcode")):
                barcode_val = str(meta["Barcode"]).strip()
            if not barcode_val:
                m_bc = re.match(r"^(barcode\d+)_", str(cluster_id))
                if m_bc:
                    barcode_val = m_bc.group(1)
            if barcode_val:
                lines.append(html.escape(barcode_val))
            continue
        if not meta or field not in meta or pd.isna(meta[field]):
            continue
        val = meta[field]
        if field == "Percent_Identity":
            try:
                lines.append(html.escape(f"{float(val):.1f}%"))
            except (TypeError, ValueError):
                lines.append(html.escape(str(val)))
        else:
            s = str(val).strip()
            if s:
                lines.append(html.escape(s))

    if not lines:
        lines.append(html.escape(str(cluster_id)))

    display = "<br>".join(lines)
    tip_parts = [f"Cluster: {cluster_id}"]
    if meta:
        for key, label in LABEL_FIELD_LABELS.items():
            if key == "Cluster_Name" or key not in meta or pd.isna(meta[key]):
                continue
            val = str(meta[key]).strip()
            if val:
                tip_parts.append(f"{label}: {val}")
        if "Species_Name" in meta and pd.notna(meta["Species_Name"]):
            raw = str(meta["Species_Name"]).strip()
            if raw and raw.upper() != "NA":
                tip_parts.append(f"Hit: {raw[:120]}{'…' if len(raw) > 120 else ''}")
    return display, html.escape(" · ".join(tip_parts))


def render_view_tab(
    fasta_path: str,
    selected: list[str],
    sample: str,
    mode: str,
) -> None:
    if not os.path.isfile(fasta_path):
        st.info(
            f"No alignment at `{fasta_path}`. Run MAFFT on the **Run** tab first (needs ≥2 sequences)."
        )
        legacy = "/data/intermediate_data/mafft_alignment.fasta"
        if os.path.exists(legacy):
            with st.expander("Legacy global alignment"):
                sequences = load_alignment(legacy, os.path.getmtime(legacy))
                visualize_mafft_alignment(
                    sequences,
                    list(sequences.keys())[:15],
                    None,
                    ["Cluster_Name"],
                    key_prefix="mafft_legacy",
                )
        return

    if not selected:
        st.warning("No sequences selected. Pick clusters on the **Sequence sets** tab.")
        return

    mtime = os.path.getmtime(fasta_path)
    sequences = load_alignment(fasta_path, mtime)
    in_aln = [s for s in selected if s in sequences]
    missing = [s for s in selected if s not in sequences]

    if missing:
        st.caption(f"{len(missing)} selected IDs not in this alignment (run MAFFT again if the BLAST set changed).")

    st.caption(f"`{os.path.basename(fasta_path)}` · showing **{len(in_aln)}** sequences")

    dl1, dl2 = st.columns(2)
    with dl1:
        st.download_button(
            "Download full alignment",
            alignment_fasta_bytes(fasta_path),
            file_name=os.path.basename(fasta_path),
            mime="text/plain",
            width="stretch",
        )
    with dl2:
        if in_aln:
            st.download_button(
                "Download selected subset (FASTA)",
                subset_fasta_bytes(sequences, in_aln),
                file_name=f"{sample}_{mode}_subset.fasta",
                mime="text/plain",
                width="stretch",
            )

    meta_lookup = cluster_metadata_lookup(sample, mode)
    extra_labels = st.multiselect(
        "Also show in row labels (with cluster ID)",
        options=LABEL_EXTRA_FIELDS,
        default=["Genus", "Species"],
        format_func=lambda f: LABEL_FIELD_LABELS.get(f, f),
        key=f"mafft_view_label_fields_{sample}_{mode}",
        help="Shown in the fixed left column next to each sequence.",
    )
    label_fields = ["Cluster_Name"] + [f for f in extra_labels if f in LABEL_EXTRA_FIELDS]

    visualize_mafft_alignment(
        sequences,
        in_aln,
        meta_lookup,
        label_fields,
        key_prefix=f"mafft_view_{sample}_{mode}",
    )


def visualize_mafft_alignment(
    sequences: dict[str, str],
    selected_seqs: list[str],
    meta_by_cluster: dict[str, dict] | None = None,
    label_fields: list[str] | None = None,
    *,
    key_prefix: str = "mafft_view",
) -> None:
    colors = {
        "A": "#9F1B12",
        "T": "#13781b",
        "C": "#3c82d8",
        "G": "#dea023",
        "N": "#7A7A7A",
    }

    if not sequences:
        st.error("The alignment is empty.")
        return

    selected_seqs = [s for s in selected_seqs if s in sequences]
    if not selected_seqs:
        st.warning("No sequences to display.")
        return

    aln_length = max(len(sequences[s]) for s in selected_seqs)

    ctrl1, ctrl2 = st.columns([1, 1])
    with ctrl1:
        default_window = min(500, aln_length)
        start_pos, end_pos = st.slider(
            "Alignment region (bp)",
            1,
            aln_length,
            (1, default_window),
            key=f"{key_prefix}_window",
        )
    with ctrl2:
        highlight_mutations = st.checkbox(
            "Highlight mutations only",
            value=True,
            key=f"{key_prefix}_highlight_mut",
        )
        st.caption(f"Alignment length: **{aln_length}** bp")

    window_seqs = [sequences[s][start_pos - 1 : end_pos] for s in selected_seqs]
    consensus_slice = ""
    for i in range(end_pos - start_pos + 1):
        col_bases = [s[i] for s in window_seqs if i < len(s)]
        if not col_bases:
            consensus_slice += "-"
        else:
            consensus_slice += Counter(col_bases).most_common(1)[0][0]

    html_out = "<div class='aln-container'>"
    if label_fields is None:
        label_fields = ["Cluster_Name"]
    if meta_by_cluster is None:
        meta_by_cluster = {}

    label_header = "Labels"
    if len(label_fields) > 1:
        label_header = "Cluster · " + " · ".join(
            LABEL_FIELD_LABELS.get(f, f) for f in label_fields if f != "Cluster_Name"
        )

    html_out += (
        f"<div class='aln-pos'><strong class='sticky-label'>"
        f"{html.escape(label_header)}</strong> [{start_pos} ... {end_pos}]</div>"
    )
    html_out += (
        "<div class='aln-consensus'><strong class='sticky-label'>CONSENSUS</strong> "
    )
    for char in consensus_slice:
        cls = "base gap" if char == "-" else "base"
        style = (
            f"style='background-color: {colors.get(char, '#ffffff')}; font-weight: bold;'"
            if char != "-"
            else ""
        )
        html_out += f"<span class='{cls}' {style}>{char}</span>"
    html_out += "</div>"

    for seq_id in selected_seqs:
        seq_slice = sequences[seq_id][start_pos - 1 : end_pos]
        meta = resolve_cluster_meta(seq_id, meta_by_cluster)
        display_label, tip = format_row_label(seq_id, meta, label_fields)
        html_out += (
            f"<div class='aln-seq'><span class='copyable-seq sticky-label mafft-row-label' "
            f"contenteditable='true' spellcheck='false' title='{tip}'>"
            f"{display_label}</span> "
        )
        for i, char in enumerate(seq_slice):
            if char == "-":
                html_out += "<span class='base gap'>-</span>"
            elif highlight_mutations and i < len(consensus_slice) and char == consensus_slice[i]:
                html_out += "<span class='base match-dot'>.</span>"
            else:
                html_out += (
                    f"<span class='base' style='background-color: {colors.get(char, '#ffffff')};'>"
                    f"{char}</span>"
                )
        html_out += "</div>"

    html_out += "</div>"
    with st.container(border=True):
        st.markdown(html_out, unsafe_allow_html=True)




def barcode_summary_row(sample: str, mode: str) -> dict:
    """One row of cross-barcode comparison metrics."""
    stats = preflight_stats(sample, mode)
    row: dict = {
        "Barcode": sample,
        "Clusters": stats["n_tsv"],
        "FASTAs found": stats["n_fasta"],
        "Missing FASTA": stats["n_missing"],
        "MAFFT raw": "yes" if stats["has_raw"] else "—",
        "MAFFT trimmed": "yes" if stats["has_trimmed"] else "—",
    }
    tsv = blast_tsv_path(sample, mode)
    if not os.path.isfile(tsv):
        return row

    df = load_blast_table(tsv, os.path.getmtime(tsv))
    if df.empty:
        return row

    if "Genus" in df.columns:
        row["Genera"] = int(df["Genus"].nunique())
    if "Percent_Identity" in df.columns:
        row["Mean identity %"] = round(float(df["Percent_Identity"].mean()), 1)
    if "Confidence" in df.columns:
        conf = df["Confidence"].astype(str)
        row["High + medium"] = int(conf.isin(["high", "medium"]).sum())
        row["Genus-only"] = int(conf.eq("genus_only").sum())
        row["Fail"] = int(conf.eq("fail").sum())
    return row


def genus_sets_by_barcode(samples: list[str], mode: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for sample in samples:
        tsv = blast_tsv_path(sample, mode)
        if not os.path.isfile(tsv):
            continue
        df = load_blast_table(tsv, os.path.getmtime(tsv))
        if df.empty or "Genus" not in df.columns:
            continue
        out[sample] = set(df["Genus"].dropna().astype(str))
    return out


def render_compare_tab(all_samples: list[str]) -> None:
    st.caption(
        "Compare cluster counts, MAFFT status, and genus composition across barcodes "
        "(same BLAST strict/full set for each)."
    )

    default_pick = all_samples[: min(4, len(all_samples))]
    picked = st.multiselect(
        "Barcodes to compare",
        all_samples,
        default=default_pick,
        key="mafft_compare_samples",
    )
    if len(picked) < 2:
        st.info("Select at least **2 barcodes** to compare.")
        return

    cmp_mode = st.radio(
        "BLAST set for comparison",
        options=["strict", "full"],
        format_func=lambda m: "Strict (high + medium)"
        if m == "strict"
        else "Full (all BLAST rows)",
        horizontal=True,
        key="mafft_compare_mode",
    )

    summary = pd.DataFrame([barcode_summary_row(s, cmp_mode) for s in picked])
    st.markdown("**Summary**")
    st.dataframe(summary, hide_index=True, width="stretch")

    genus_sets = genus_sets_by_barcode(picked, cmp_mode)
    if not genus_sets:
        st.warning("No genus data available for the selected barcodes.")
        return

    shared = set.intersection(*genus_sets.values()) if len(genus_sets) > 1 else set()
    union = set.union(*genus_sets.values())

    m1, m2, m3 = st.columns(3)
    m1.metric("Barcodes compared", len(picked))
    m2.metric("Shared genera", len(shared))
    m3.metric("Total unique genera", len(union))

    if shared:
        shown = sorted(shared)[:40]
        extra = len(shared) - len(shown)
        suffix = f" … +{extra} more" if extra > 0 else ""
        st.markdown(f"**In all barcodes:** {', '.join(shown)}{suffix}")

    with st.expander("Genera unique to one barcode"):
        for sample in picked:
            if sample not in genus_sets:
                st.markdown(f"**{sample}** — no data")
                continue
            others = set.union(
                *(genus_sets[b] for b in picked if b != sample and b in genus_sets)
            )
            unique = genus_sets[sample] - others
            preview = ", ".join(sorted(unique)[:25]) if unique else "—"
            st.markdown(f"**{sample}** ({len(unique)}): {preview}")

    chart_rows = []
    for sample in picked:
        tsv = blast_tsv_path(sample, cmp_mode)
        if not os.path.isfile(tsv):
            continue
        df = load_blast_table(tsv, os.path.getmtime(tsv))
        if df.empty or "Genus" not in df.columns:
            continue
        for genus, n in df["Genus"].value_counts().head(10).items():
            chart_rows.append(
                {"Barcode": sample, "Genus": str(genus), "Clusters": int(n)}
            )

    if chart_rows:
        st.markdown("**Top genera per barcode**")
        chart_df = pd.DataFrame(chart_rows)
        fig = px.bar(
            chart_df,
            x="Genus",
            y="Clusters",
            color="Barcode",
            barmode="group",
            title=f"Top genera · {cmp_mode} BLAST set",
        )
        fig.update_layout(**PLOTLY_LAYOUT, height=420, xaxis_tickangle=-35)
        st.plotly_chart(fig, width="stretch")

    # Strict vs full cluster counts for same barcodes (optional insight)
    if st.checkbox("Show strict vs full cluster counts", value=False, key="mafft_cmp_dual"):
        dual_rows = []
        for sample in picked:
            dual_rows.append(
                {
                    "Barcode": sample,
                    "strict": preflight_stats(sample, "strict")["n_tsv"],
                    "full": preflight_stats(sample, "full")["n_tsv"],
                }
            )
        dual = pd.DataFrame(dual_rows)
        st.dataframe(dual, hide_index=True, width="stretch")
        melt = dual.melt(id_vars=["Barcode"], var_name="Set", value_name="Clusters")
        fig2 = px.bar(
            melt,
            x="Barcode",
            y="Clusters",
            color="Set",
            barmode="group",
            title="Clusters per barcode · strict vs full",
        )
        fig2.update_layout(**PLOTLY_LAYOUT, height=360)
        st.plotly_chart(fig2, width="stretch")


@st.cache_data(show_spinner=False)
def load_organism_manifest(mtime: float) -> pd.DataFrame:
    del mtime
    if not os.path.isfile(ORGANISM_MANIFEST):
        return pd.DataFrame()
    try:
        return pd.read_csv(ORGANISM_MANIFEST, sep="\t")
    except (pd.errors.ParserError, OSError, ValueError):
        return pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_organism_clusters(tsv_path: str, mtime: float) -> pd.DataFrame:
    del mtime
    if not os.path.isfile(tsv_path):
        return pd.DataFrame()
    return pd.read_csv(tsv_path, sep="\t")


def organism_metadata_lookup(clusters_df: pd.DataFrame) -> dict[str, dict]:
    lookup: dict[str, dict] = {}
    if clusters_df.empty:
        return lookup
    for _, row in clusters_df.iterrows():
        barcode = str(row.get("Barcode", "")).strip()
        cluster = str(row.get("Cluster_Name", "")).strip()
        meta = row.to_dict()
        if barcode and cluster:
            lookup[f"{barcode}_{cluster}"] = meta
        if cluster:
            lookup[cluster] = meta
            if cluster.isdigit():
                lookup[str(int(cluster))] = meta
    return lookup


def organism_alignment_paths(slug_row: pd.Series) -> dict[str, str]:
    return {
        "raw": str(slug_row.get("Raw_alignment", "") or ""),
        "trimmed": str(slug_row.get("Trimmed_alignment", "") or ""),
        "clusters": str(slug_row.get("Cluster_list", "") or ""),
    }


def organism_mafft_is_running() -> bool:
    return os.path.isfile(ORGANISM_MAFFT_RUNNING_FLAG)


def read_organism_mafft_log_tail(n: int = 25) -> str:
    if not os.path.isfile(ORGANISM_MAFFT_LOG):
        return ""
    try:
        with open(ORGANISM_MAFFT_LOG, encoding="utf-8", errors="replace") as f:
            return "".join(f.readlines()[-n:])
    except OSError:
        return ""


def launch_organism_mafft() -> None:
    os.makedirs(ORGANISM_MAFFT_DIR, exist_ok=True)
    with open(ORGANISM_MAFFT_RUNNING_FLAG, "w", encoding="utf-8") as f:
        f.write("running\n")
    shell_cmd = (
        f"python3 {ORGANISM_MAFFT_SCRIPT} >> {ORGANISM_MAFFT_LOG} 2>&1; "
        f"rm -f {ORGANISM_MAFFT_RUNNING_FLAG}"
    )
    with open(ORGANISM_MAFFT_LOG, "w", encoding="utf-8") as log:
        log.write(f"$ {shell_cmd}\n\n")
    subprocess.Popen(
        ["bash", "-c", shell_cmd],
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    load_organism_manifest.clear()
    load_organism_clusters.clear()
    load_alignment.clear()


def render_organism_view_tab(
    fasta_path: str,
    clusters_path: str,
    organism_label: str,
    selected: list[str],
) -> None:
    if not fasta_path or not os.path.isfile(fasta_path):
        st.info(
            "No alignment file for this organism. Build alignments on the **By organism** tab."
        )
        return

    if not selected:
        st.warning("No sequences selected. Adjust barcode filters or sequence list.")
        return

    mtime = os.path.getmtime(fasta_path)
    sequences = load_alignment(fasta_path, mtime)
    in_aln = [s for s in selected if s in sequences]
    missing = [s for s in selected if s not in sequences]

    if missing:
        st.caption(f"{len(missing)} selected IDs not found in the alignment file.")

    st.caption(
        f"`{os.path.basename(fasta_path)}` · **{organism_label}** · "
        f"showing {len(in_aln)} / {len(sequences)} sequences"
    )

    dl1, dl2 = st.columns(2)
    with dl1:
        st.download_button(
            "Download full alignment",
            alignment_fasta_bytes(fasta_path),
            file_name=os.path.basename(fasta_path),
            mime="text/plain",
            width="stretch",
            key=f"org_dl_full_{organism_label}",
        )
    with dl2:
        if in_aln:
            safe_name = re.sub(r"[^\w]+", "_", organism_label)[:60]
            st.download_button(
                "Download selected subset (FASTA)",
                subset_fasta_bytes(sequences, in_aln),
                file_name=f"{safe_name}_subset.fasta",
                mime="text/plain",
                width="stretch",
                key=f"org_dl_sub_{organism_label}",
            )

    clusters_df = (
        load_organism_clusters(clusters_path, os.path.getmtime(clusters_path))
        if clusters_path and os.path.isfile(clusters_path)
        else pd.DataFrame()
    )
    meta_lookup = organism_metadata_lookup(clusters_df)

    extra_labels = st.multiselect(
        "Also show in row labels",
        options=LABEL_EXTRA_FIELDS + ["Barcode"],
        default=["Barcode", "Genus", "Species"],
        format_func=lambda f: "Barcode" if f == "Barcode" else LABEL_FIELD_LABELS.get(f, f),
        key=f"org_view_labels_{organism_label}",
    )
    label_fields = ["Cluster_Name"] + [f for f in extra_labels if f != "Cluster_Name"]
    org_key = re.sub(r"[^\w]+", "_", organism_label)[:80]
    visualize_mafft_alignment(
        sequences,
        in_aln,
        meta_lookup,
        label_fields,
        key_prefix=f"mafft_org_{org_key}",
    )


def render_by_organism_tab() -> None:
    st.caption(
        "Cross-barcode MAFFT alignments per top species (`intermediate_data/mafft/by_organism/`). "
        "One alignment per taxon — sequence IDs are `barcodeXX_cluster`."
    )

    if organism_mafft_is_running():
        st.warning("Organism alignments are running in the background.")
        if st.button("Refresh", key="organism_mafft_refresh_running"):
            st.rerun()
        st.code(read_organism_mafft_log_tail(30) or "(empty)", language="log")
        return

    manifest_mtime = (
        os.path.getmtime(ORGANISM_MANIFEST) if os.path.isfile(ORGANISM_MANIFEST) else 0.0
    )
    manifest = load_organism_manifest(manifest_mtime)
    ok_rows = (
        manifest[manifest["Status"].astype(str) == "ok"] if not manifest.empty else pd.DataFrame()
    )

    c1, c2 = st.columns([1, 3])
    with c1:
        if st.button("Build top-10 alignments", type="primary", key="organism_mafft_run_btn"):
            launch_organism_mafft()
            st.rerun()
    with c2:
        st.caption(
            "MAFFT on all clusters assigned to each of the top 10 species (full BLAST set)."
        )

    if manifest.empty or ok_rows.empty:
        st.info(
            "No organism alignments yet. Click **Build top-10 alignments** to generate them."
        )
        with st.expander("Run log", expanded=False):
            st.code(read_organism_mafft_log_tail(40) or "No log yet.", language="log")
        return

    organism_options = ok_rows["Organism"].astype(str).tolist()
    default_org = organism_options[0]
    rank_map = dict(zip(ok_rows["Organism"].astype(str), ok_rows["Rank"]))

    picked_org = st.selectbox(
        "Organism",
        organism_options,
        format_func=lambda o: f"#{int(rank_map.get(o, 0))} {o}",
        key="mafft_organism_pick",
    )
    row = ok_rows[ok_rows["Organism"].astype(str) == picked_org].iloc[0]
    paths = organism_alignment_paths(row)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Sequences", int(row.get("Sequences_aligned", 0)))
    m2.metric("Expected clusters", int(row.get("Clusters_expected", 0)))
    m3.metric("Missing FASTA", int(row.get("Sequences_missing", 0)))
    m4.metric("Status", str(row.get("Status", "—")))

    clusters_path = paths["clusters"]
    clusters_df = (
        load_organism_clusters(clusters_path, os.path.getmtime(clusters_path))
        if clusters_path and os.path.isfile(clusters_path)
        else pd.DataFrame()
    )

    sub_run, sub_clusters, sub_view = st.tabs(["Summary", "Cluster list", "View alignment"])

    with sub_run:
        st.markdown("**Selected organism**")
        st.dataframe(
            pd.DataFrame([row]).drop(columns=["Raw_alignment", "Trimmed_alignment"], errors="ignore"),
            hide_index=True,
            width="stretch",
        )
        if not manifest.empty:
            st.markdown("**All organism alignments**")
            show_manifest = manifest.copy()
            for col in ("Raw_alignment", "Trimmed_alignment", "Cluster_list"):
                if col in show_manifest.columns:
                    show_manifest[col] = show_manifest[col].apply(
                        lambda p: os.path.basename(str(p)) if pd.notna(p) and str(p) else ""
                    )
            st.dataframe(show_manifest, hide_index=True, width="stretch")
        with st.expander("Run log", expanded=False):
            st.code(read_organism_mafft_log_tail(60) or "No log yet.", language="log")

    with sub_clusters:
        if clusters_df.empty:
            st.warning("Cluster list TSV not found for this organism.")
        else:
            st.caption(f"`{os.path.basename(clusters_path)}` · {len(clusters_df)} clusters")
            show_cols = [
                c
                for c in [
                    "Barcode",
                    "Cluster_Name",
                    "Genus",
                    "Species",
                    "Percent_Identity",
                    "Query_Coverage(%)",
                    "Confidence",
                    "Single_Organism",
                ]
                if c in clusters_df.columns
            ]
            st.dataframe(
                clusters_df[show_cols] if show_cols else clusters_df,
                hide_index=True,
                width="stretch",
            )
            if "Barcode" in clusters_df.columns:
                bc_counts = clusters_df["Barcode"].value_counts().reset_index()
                bc_counts.columns = ["Barcode", "Clusters"]
                fig = px.bar(
                    bc_counts,
                    x="Barcode",
                    y="Clusters",
                    title=f"Clusters per barcode · {picked_org}",
                )
                fig.update_layout(**PLOTLY_LAYOUT, height=320)
                st.plotly_chart(fig, width="stretch", key=f"org_bc_chart_{picked_org}")

    with sub_view:
        align_type = st.radio(
            "Alignment file",
            ["Trimmed (trimAl)", "Raw (MAFFT)"],
            horizontal=True,
            key="mafft_organism_view_type",
        )
        fasta_path = paths["trimmed"] if align_type.startswith("Trimmed") else paths["raw"]
        if not fasta_path or not os.path.isfile(fasta_path):
            fasta_path = paths["raw"]

        all_seq_ids: list[str] = []
        if fasta_path and os.path.isfile(fasta_path):
            all_seq_ids = list(load_alignment(fasta_path, os.path.getmtime(fasta_path)).keys())

        barcodes = sorted(
            {
                str(b)
                for b in clusters_df.get("Barcode", pd.Series(dtype=str)).dropna().astype(str)
            }
        )
        if not barcodes and all_seq_ids:
            barcodes = sorted(
                {m.group(1) for s in all_seq_ids if (m := re.match(r"^(barcode\d+)_", s))}
            )

        f1, f2 = st.columns(2)
        with f1:
            barcode_filter = st.multiselect(
                "Barcodes to show",
                barcodes,
                default=barcodes,
                key=f"org_view_barcodes_{picked_org}",
            )
        with f2:
            max_show = st.number_input(
                "Max sequences in viewer",
                min_value=5,
                max_value=max(5, len(all_seq_ids)),
                value=min(40, len(all_seq_ids)) if all_seq_ids else 40,
                key=f"org_view_max_{picked_org}",
            )

        filtered_ids = []
        for seq_id in all_seq_ids:
            m = re.match(r"^(barcode\d+)_(.+)$", seq_id)
            if m and barcode_filter and m.group(1) not in barcode_filter:
                continue
            filtered_ids.append(seq_id)

        if len(filtered_ids) > max_show:
            st.caption(
                f"Showing first **{max_show}** of **{len(filtered_ids)}** sequences "
                "(large cross-barcode sets — narrow barcodes or download full FASTA)."
            )
            filtered_ids = filtered_ids[: int(max_show)]

        render_organism_view_tab(
            fasta_path,
            clusters_path,
            picked_org,
            filtered_ids,
        )


def show_mafft_page() -> None:
    page_hero(
        "Multiple Sequence Alignment (MAFFT)",
        "BLAST-driven sequence sets per barcode · strict or full · run, curate, and view",
    )

    samples = list_blast_samples()
    if not samples:
        st.warning("No BLAST summaries in `blast_results/`. Run BLAST before MAFFT.")
        return

    top1, top2 = st.columns([2, 1])
    with top1:
        sample = st.selectbox("Sample (barcode)", samples, key="mafft_sample")
        blast_mode = st.radio(
            "BLAST sequence set",
            options=["strict", "full"],
            format_func=lambda m: "Strict (high + medium confidence)"
            if m == "strict"
            else "Full (all BLAST rows)",
            horizontal=True,
            key="mafft_blast_mode",
        )
    with top2:
        stats = preflight_stats(sample, blast_mode)
        st.metric("In alignment", "✓" if stats["has_raw"] else "—")
        st.caption(
            f"{'trimmed available' if stats['has_trimmed'] else 'raw only'}"
            if stats["has_raw"]
            else "not run yet"
        )

    tsv_path = blast_tsv_path(sample, blast_mode)
    df = (
        load_blast_table(tsv_path, os.path.getmtime(tsv_path))
        if os.path.isfile(tsv_path)
        else pd.DataFrame()
    )

    paths = alignment_paths(sample, blast_mode)
    tab_run, tab_sets, tab_view, tab_compare, tab_organism = st.tabs(
        ["Run", "Sequence sets", "View", "Compare barcodes", "By organism"]
    )

    with tab_run:
        render_run_tab(sample, blast_mode)

    with tab_sets:
        selected = render_sets_tab(sample, blast_mode, df)

    with tab_view:
        align_type = st.radio(
            "Alignment file",
            ["Raw (MAFFT)", "Trimmed (trimAl)"],
            horizontal=True,
            key="mafft_view_type",
        )
        fasta_path = paths["raw"] if align_type.startswith("Raw") else paths["trimmed"]
        render_view_tab(
            fasta_path,
            st.session_state.get(selection_key(sample, blast_mode), selected),
            sample,
            blast_mode,
        )

    with tab_compare:
        render_compare_tab(samples)

    with tab_organism:
        render_by_organism_tab()
