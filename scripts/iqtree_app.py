import glob
import os
import shutil
import subprocess

import streamlit as st

INTERMEDIATE_DIR = "/data/intermediate_data"
PHASE5_DIR = f"{INTERMEDIATE_DIR}/phase5"
MAFFT_DIR = f"{INTERMEDIATE_DIR}/mafft"
MAFFT_SCRIPT = "/data/scripts/run_mafft.sh"
PHASE5_SCRIPT = "/data/scripts/run_phase5_representatives.py"
PHASE5_PREP_FLAG = f"{PHASE5_DIR}/representatives_prep_running.flag"
R_SCRIPT = "/data/scripts/plot_tree.R"
IQTREE_THREADS = os.environ.get("IQTREE_THREADS", "AUTO")

REPRESENTATIVES_TSV = f"{PHASE5_DIR}/representatives.tsv"
REPRESENTATIVES_FASTA = f"{PHASE5_DIR}/representatives.fasta"
REPRESENTATIVES_MAFFT = f"{PHASE5_DIR}/representatives_mafft.fasta"
REPRESENTATIVES_MAFFT_TRIM = f"{PHASE5_DIR}/representatives_mafft_trimmed.fasta"
TOP10_IQ_PREFIX = REPRESENTATIVES_MAFFT_TRIM

SCOPE_LABELS = {
    "top10": "Top 10 species (all BLAST calls)",
    "top10_single_species": "Top 10 (single-species QC)",
}


import re

import pandas as pd

from mafft_app import blast_tsv_path, load_blast_table, parse_taxonomy_from_species_name





def tree_blast_context(treefile: str) -> tuple[str, str] | None:
    """Parse sample + strict/full from IQ-TREE output name."""
    base = os.path.basename(treefile)
    match = re.match(
        r"sample_(.+?)_(strict|full)_.+_trimmed\.fasta\.treefile$",
        base,
    )
    if match:
        return match.group(1), match.group(2)
    return None


def _tip_label_text(cluster: str, genus: str, species: str) -> str:
    g = (genus or "").strip() or "Unknown"
    sp = (species or "").strip()
    if sp and sp.lower() not in g.lower():
        return f"{cluster}\n{g} — {sp}"
    return f"{cluster}\n{g}"


def tip_meta_sidecar_path(treefile: str) -> str:
    """Sidecar CSV path aligned with IQ-TREE prefix (strip .treefile suffix)."""
    if treefile.endswith(".treefile"):
        return treefile[: -len(".treefile")] + ".tip_meta.csv"
    return f"{treefile}.tip_meta.csv"


def is_top10_representative_tree(treefile: str) -> bool:
    return os.path.basename(treefile) == "representatives_mafft_trimmed.fasta.treefile"


def _display_tip_label(genus: str, species: str) -> str:
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


def _clean_genus_label(genus: str) -> str:
    g = str(genus or "").strip()
    if not g or g.lower() in ("nan", "na", "unknown", ""):
        return "Unknown"
    return g.replace("_", " ")


def write_top10_tip_meta() -> str | None:
    """Build tip_meta.csv for the top-10 representative tree from representatives.tsv."""
    if not os.path.isfile(REPRESENTATIVES_TSV):
        return None
    try:
        df = pd.read_csv(REPRESENTATIVES_TSV, sep="\t")
    except (pd.errors.ParserError, OSError, ValueError):
        return None
    if df.empty or "Seq_ID" not in df.columns:
        return None

    rows = []
    for _, row in df.iterrows():
        genus = _clean_genus_label(row.get("Genus", ""))
        species = str(row.get("Species", "") or "").strip()
        if species.lower() in ("nan", ""):
            species = ""
        tip_label = str(row.get("Tip_Label", "") or "").strip()
        if not tip_label or tip_label.lower() == "nan":
            tip_label = _display_tip_label(genus, species)
        rows.append(
            {
                "cluster": str(row["Seq_ID"]),
                "genus": genus,
                "species": species,
                "tip_label": tip_label,
                "source_cluster": str(row.get("Source_Cluster", "")),
            }
        )

    if not rows:
        return None

    path = tip_meta_sidecar_path(f"{TOP10_IQ_PREFIX}.treefile")
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def build_tip_metadata_csv(treefile: str) -> str | None:
    """Resolve tip metadata CSV for R/ggtree (species labels, genus colors)."""
    sidecar = tip_meta_sidecar_path(treefile)
    if os.path.isfile(sidecar):
        return sidecar

    if is_top10_representative_tree(treefile):
        return write_top10_tip_meta()

    ctx = tree_blast_context(treefile)
    if not ctx:
        return None
    sample, mode = ctx
    tsv_path = blast_tsv_path(sample, mode)
    if not os.path.isfile(tsv_path):
        return None

    df = load_blast_table(tsv_path, os.path.getmtime(tsv_path))
    if df.empty or "Cluster_Name" not in df.columns:
        return None

    rows = []
    for _, row in df.iterrows():
        cluster = str(row["Cluster_Name"]).strip()
        genus = str(row.get("Genus", "") or "").strip()
        species = str(row.get("Species", "") or "").strip()
        if (not genus or genus == "nan") and "Species_Name" in row:
            parsed = parse_taxonomy_from_species_name(row.get("Species_Name"))
            genus = parsed.get("Genus", genus) or genus
            species = parsed.get("Species", species) or species
        rows.append(
            {
                "cluster": cluster,
                "genus": genus or "Unknown",
                "species": species,
                "tip_label": _tip_label_text(cluster, genus, species),
            }
        )

    if not rows:
        return None

    csv_path = tip_meta_sidecar_path(treefile)
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


def page_hero(title: str, subtitle: str) -> None:
    st.markdown(
        f'<div class="ff-page-hero"><h1>{title}</h1><p>{subtitle}</p></div>',
        unsafe_allow_html=True,
    )


def find_tree_files() -> list[str]:
    patterns = [
        os.path.join(INTERMEDIATE_DIR, "*.treefile"),
        os.path.join(INTERMEDIATE_DIR, "**", "*.treefile"),
    ]
    found: list[str] = []
    for pattern in patterns:
        found.extend(glob.glob(pattern, recursive=True))
    return sorted(set(found))


def tree_asset_paths(treefile: str) -> dict[str, str]:
    base = treefile.replace(".treefile", "")
    return {
        "rect_png": f"{base}_tree_rect.png",
        "circ_png": f"{base}_tree_circ.png",
        "rect_pdf": f"{base}_tree_rect.pdf",
        "circ_pdf": f"{base}_tree_circ.pdf",
        "legacy_png": f"{base}_tree.png",
        "legacy_pdf": f"{base}_tree.pdf",
        "report": f"{base}.iqtree",
        "log": f"{base}.log",
    }


def count_fasta_seqs(path: str) -> int:
    with open(path, encoding="utf-8") as f:
        return sum(1 for line in f if line.startswith(">"))


def render_tree_plot(treefile: str, layout: str) -> str:
    """Run R script; return path to PNG."""
    assets = tree_asset_paths(treefile)
    out_png = assets["circ_png"] if layout == "circular" else assets["rect_png"]
    meta_csv = build_tip_metadata_csv(treefile)
    cmd = ["Rscript", R_SCRIPT, treefile, out_png, layout]
    if meta_csv:
        cmd.append(meta_csv)
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return out_png


def ensure_tree_plots(treefile: str, layouts: tuple[str, ...] = ("rectangular", "circular")) -> None:
    assets = tree_asset_paths(treefile)
    for layout in layouts:
        target = assets["circ_png"] if layout == "circular" else assets["rect_png"]
        if not os.path.isfile(target):
            render_tree_plot(treefile, layout)


def display_tree_image(treefile: str, layout_choice: str) -> None:
    assets = tree_asset_paths(treefile)
    if layout_choice == "Circular":
        img = assets["circ_png"]
        pdf = assets["circ_pdf"]
    else:
        img = assets["rect_png"] if os.path.isfile(assets["rect_png"]) else assets["legacy_png"]
        pdf = assets["rect_pdf"] if os.path.isfile(assets["rect_pdf"]) else assets["legacy_pdf"]

    if not os.path.isfile(img):
        with st.spinner(f"Rendering {layout_choice.lower()} tree…"):
            try:
                render_tree_plot(
                    treefile,
                    "circular" if layout_choice == "Circular" else "rectangular",
                )
            except subprocess.CalledProcessError as exc:
                st.error(f"R/ggtree failed: {exc.stderr or exc}")
                return

    if os.path.isfile(img):
        st.image(img, use_container_width=True)
    else:
        st.warning("Tree image could not be generated.")

    if os.path.isfile(pdf):
        with open(pdf, "rb") as f:
            st.download_button(
                f"Download {layout_choice} PDF",
                f.read(),
                file_name=os.path.basename(pdf),
                mime="application/pdf",
                width="stretch",
            )


def build_iqtree_command(
    alignment: str,
    prefix: str,
    *,
    model: str,
    bootstrap: int,
    alrt: int,
) -> list[str]:
    cmd = [
        "iqtree2" if shutil.which("iqtree2") else "iqtree",
        "-s",
        alignment,
        "-T",
        str(IQTREE_THREADS),
        "-pre",
        prefix,
        "-redo",
    ]
    if model == "MFP":
        cmd.extend(["-m", "MFP"])
    else:
        cmd.extend(["-m", "HKY+F+R5"])

    if bootstrap > 0:
        cmd.extend(["-bb", str(bootstrap), "-bnni"])
    if alrt > 0:
        cmd.extend(["-alrt", str(alrt)])
    return cmd


def iqtree_running(prefix: str) -> bool:
    treefile = f"{prefix}.treefile"
    logfile = f"{prefix}.log"
    return os.path.isfile(logfile) and not os.path.isfile(treefile)


def launch_iqtree(cmd: list[str]) -> None:
    log_path = cmd[cmd.index("-pre") + 1] + ".log"
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.Popen(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )


def representatives_prep_running() -> bool:
    return os.path.isfile(PHASE5_PREP_FLAG)


def read_prep_log_tail(n: int = 20) -> str:
    log_path = f"{PHASE5_DIR}/representatives_prep.log"
    if not os.path.isfile(log_path):
        return ""
    try:
        with open(log_path, encoding="utf-8", errors="replace") as f:
            return "".join(f.readlines()[-n:])
    except OSError:
        return ""


def launch_representatives_prep(scope: str) -> None:
    os.makedirs(PHASE5_DIR, exist_ok=True)
    with open(PHASE5_PREP_FLAG, "w", encoding="utf-8") as f:
        f.write("running\n")
    log_path = f"{PHASE5_DIR}/representatives_prep.log"
    shell_cmd = (
        f"python3 {PHASE5_SCRIPT} --scope {scope} >> {log_path} 2>&1; "
        f"rm -f {PHASE5_PREP_FLAG}"
    )
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"$ {shell_cmd}\n\n")
    subprocess.Popen(
        ["bash", "-c", shell_cmd],
        start_new_session=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def load_representatives_table() -> pd.DataFrame:
    if not os.path.isfile(REPRESENTATIVES_TSV):
        return pd.DataFrame()
    try:
        return pd.read_csv(REPRESENTATIVES_TSV, sep="\t")
    except (pd.errors.ParserError, OSError, ValueError):
        return pd.DataFrame()


def top10_alignment_path() -> str | None:
    if os.path.isfile(REPRESENTATIVES_MAFFT_TRIM):
        return REPRESENTATIVES_MAFFT_TRIM
    if os.path.isfile(REPRESENTATIVES_MAFFT):
        return REPRESENTATIVES_MAFFT
    return None


def render_iqtree_settings() -> tuple[str, int, int]:
    with st.container(border=True):
        st.markdown("**IQ-TREE settings**")
        o1, o2, o3 = st.columns(3)
        model = o1.selectbox(
            "Substitution model",
            ["MFP", "HKY+F+R5"],
            help="MFP = ModelFinder (recommended for ITS).",
            key="iqtree_model",
        )
        bootstrap = o2.selectbox(
            "Ultrafast bootstrap (UFBoot)",
            [1000, 500, 200, 0],
            index=0,
            help="Lower values run faster; 0 skips bootstrap.",
            key="iqtree_bootstrap",
        )
        alrt = o3.selectbox(
            "SH-aLRT replicates",
            [1000, 500, 0],
            index=0,
            key="iqtree_alrt",
        )
    return model, int(bootstrap), int(alrt)


def render_iqtree_run_status(iq_prefix: str, delete_key: str, refresh_key: str) -> bool:
    """Show running / finished state. Returns True if caller should stop."""
    treefile = f"{iq_prefix}.treefile"
    logfile = f"{iq_prefix}.log"

    if iqtree_running(iq_prefix):
        st.warning("IQ-TREE is running…")
        if st.button("Refresh", key=refresh_key):
            st.rerun()
        if os.path.isfile(logfile):
            with open(logfile, encoding="utf-8", errors="replace") as f:
                st.code("".join(f.readlines()[-20:]), language="log")
        return True

    if os.path.isfile(treefile):
        st.success("This run already finished. Pick it in **Tree explorer** above or delete to re-run.")
        if st.button("Delete this run", type="primary", key=delete_key):
            for path in glob.glob(f"{iq_prefix}*"):
                try:
                    os.remove(path)
                except OSError:
                    pass
            st.rerun()
        return True

    return False


def run_iqtree_on_alignment(
    align_input: str,
    iq_prefix: str,
    *,
    model: str,
    bootstrap: int,
    alrt: int,
    trimmed_copy_path: str | None = None,
) -> None:
    if trimmed_copy_path and align_input != trimmed_copy_path:
        if align_input.endswith("_trimmed.fasta") or align_input.endswith("_mafft_trimmed.fasta"):
            shutil.copy(align_input, trimmed_copy_path)
            trimmed = trimmed_copy_path
        else:
            subprocess.run(
                ["trimal", "-in", align_input, "-out", trimmed_copy_path, "-gappyout"],
                check=True,
            )
            trimmed = trimmed_copy_path
    else:
        trimmed = align_input

    cmd = build_iqtree_command(
        trimmed,
        iq_prefix,
        model=model,
        bootstrap=bootstrap,
        alrt=alrt,
    )
    launch_iqtree(cmd)


def render_top10_panel(model: str, bootstrap: int, alrt: int) -> None:
    st.caption(
        "One best ITS consensus per top species (highest identity + coverage), "
        "MAFFT-aligned across all barcodes."
    )

    if representatives_prep_running():
        st.warning("Preparing representatives and MAFFT alignment…")
        if st.button("Refresh", key="iqtree_top10_prep_refresh"):
            st.rerun()
        st.code(read_prep_log_tail(25) or "(empty)", language="log")
        return

    scope_key = st.radio(
        "Representative set",
        list(SCOPE_LABELS.keys()),
        format_func=lambda k: SCOPE_LABELS[k],
        horizontal=True,
        key="iqtree_top10_scope",
    )

    c1, c2 = st.columns([1, 3])
    with c1:
        if st.button("Prepare representatives", key="iqtree_top10_prep"):
            launch_representatives_prep(scope_key)
            st.rerun()
    with c2:
        st.caption(
            "Selects one cluster per species, builds `representatives.fasta`, runs MAFFT + trimAl."
        )

    reps_df = load_representatives_table()
    align_path = top10_alignment_path()

    if not reps_df.empty:
        write_top10_tip_meta()
        show_cols = [
            c
            for c in [
                "Representative_Rank",
                "Organism",
                "Tip_Label",
                "Seq_ID",
                "Source_Cluster",
                "Percent_Identity",
                "Query_Coverage(%)",
                "Confidence",
            ]
            if c in reps_df.columns
        ]
        st.markdown("**Selected representatives**")
        st.dataframe(
            reps_df[show_cols] if show_cols else reps_df,
            hide_index=True,
            width="stretch",
        )
    else:
        st.info("No representatives yet. Click **Prepare representatives** first.")

    m1, m2, m3 = st.columns(3)
    n_reps = len(reps_df) if not reps_df.empty else 0
    n_align = count_fasta_seqs(align_path) if align_path else 0
    m1.metric("Representatives", n_reps)
    m2.metric("Sequences in alignment", n_align)
    m3.metric("Alignment ready", "Yes" if align_path and n_align >= 3 else "No")

    iq_prefix = TOP10_IQ_PREFIX
    st.caption(f"Tree outputs: `{os.path.basename(iq_prefix)}.*`")

    if render_iqtree_run_status(
        iq_prefix,
        delete_key="iqtree_top10_delete",
        refresh_key="iqtree_top10_refresh",
    ):
        return

    if not align_path or n_align < 3:
        st.warning("Need an alignment with ≥3 sequences. Run **Prepare representatives** first.")
        return

    if st.button("Run IQ-TREE on top 10", type="primary", key="iqtree_top10_run"):
        write_top10_tip_meta()
        with st.status("Starting IQ-TREE…", expanded=True) as status:
            st.write(f"Alignment: `{os.path.basename(align_path)}` ({n_align} sequences)")
            try:
                run_iqtree_on_alignment(
                    align_path,
                    iq_prefix,
                    model=model,
                    bootstrap=bootstrap,
                    alrt=alrt,
                )
                status.update(label="IQ-TREE started in background", state="complete")
            except subprocess.CalledProcessError as exc:
                st.error(f"Pipeline failed: {exc}")
                return
        st.rerun()


def render_explorer() -> str | None:
    trees = find_tree_files()
    if not trees:
        return None

    labels = {os.path.basename(p): p for p in trees}
    selected_name = st.selectbox(
        "Select tree run",
        options=sorted(labels.keys()),
        key="iqtree_select_tree",
    )
    treefile = labels[selected_name]
    if os.path.isfile(treefile):
        try:
            ensure_tree_plots(treefile)
        except subprocess.CalledProcessError:
            pass
    assets = tree_asset_paths(treefile)

    c1, c2, c3 = st.columns(3)
    with open(treefile, encoding="utf-8") as f:
        c1.download_button(
            "Download Newick (.treefile)",
            f.read(),
            file_name=os.path.basename(treefile),
            mime="text/plain",
            width="stretch",
        )

    layout_choice = st.radio(
        "Tree layout",
        ["Rectangular", "Circular"],
        horizontal=True,
        key="iqtree_layout_view",
    )

    if c2.button("Regenerate both layouts", width="stretch"):
        for layout in ("rectangular", "circular"):
            try:
                render_tree_plot(treefile, layout)
                st.success(f"Updated {layout} plot.")
            except subprocess.CalledProcessError as exc:
                st.error(exc.stderr or str(exc))

    with st.container(border=True):
        display_tree_image(treefile, layout_choice)

    st.markdown("---")
    st.subheader("Run reports")
    if os.path.isfile(assets["report"]):
        with st.expander("IQ-TREE report (.iqtree)"):
            with open(assets["report"], encoding="utf-8", errors="replace") as f:
                st.text(f.read())
    if os.path.isfile(assets["log"]):
        with st.expander("Execution log (.log)"):
            with open(assets["log"], encoding="utf-8", errors="replace") as f:
                st.text("".join(f.readlines()[-40:]))

    return treefile


def render_generate_panel(*, has_trees: bool) -> None:
    st.markdown("### Generate new tree")

    analysis_mode = st.radio(
        "Scope",
        [
            "Top 10 representatives (cross-barcode)",
            "Single sample (BLAST + MAFFT)",
            "Legacy global alignment",
        ],
        horizontal=False,
        help="Top 10: one best sequence per leading species across all barcodes.",
    )

    model, bootstrap, alrt = render_iqtree_settings()

    if analysis_mode.startswith("Top 10"):
        render_top10_panel(model, bootstrap, alrt)
        return

    selected_sample = ""
    blast_set = "strict"
    if analysis_mode.startswith("Single"):
        consensus_root = "/data/consensus_results"
        if not os.path.isdir(consensus_root):
            st.warning("No `consensus_results/` folder found.")
            return
        samples = sorted(
            d
            for d in os.listdir(consensus_root)
            if os.path.isdir(os.path.join(consensus_root, d))
        )
        if not samples:
            st.warning("No samples found.")
            return
        selected_sample = st.selectbox("Sample (barcode)", samples, key="iqtree_sample")
        blast_set = st.radio(
            "BLAST set",
            ["strict", "full"],
            format_func=lambda m: "Strict (high + medium)" if m == "strict" else "Full",
            horizontal=True,
            key="iqtree_blast_set",
        )

    run_suffix = (
        st.text_input("Run name", value="default", key="iqtree_run_suffix")
        .strip()
        .replace(" ", "_")
    )

    if analysis_mode.startswith("Single"):
        new_prefix = f"{INTERMEDIATE_DIR}/sample_{selected_sample}_{blast_set}_{run_suffix}"
        mafft_raw = f"{MAFFT_DIR}/{selected_sample}_{blast_set}_mafft.fasta"
        mafft_trim = f"{MAFFT_DIR}/{selected_sample}_{blast_set}_mafft_trimmed.fasta"
    else:
        new_prefix = f"{INTERMEDIATE_DIR}/global_{run_suffix}"
        mafft_raw = f"{INTERMEDIATE_DIR}/mafft_alignment.fasta"
        mafft_trim = mafft_raw

    trimmed = f"{new_prefix}_trimmed.fasta"
    iq_prefix = f"{new_prefix}_trimmed.fasta"

    st.caption(f"Outputs: `{os.path.basename(iq_prefix)}.*`")

    if render_iqtree_run_status(
        iq_prefix,
        delete_key="iqtree_delete_run",
        refresh_key="iqtree_refresh_run",
    ):
        return

    if st.button("Run trimAl + IQ-TREE", type="primary", key="iqtree_run_pipeline"):
        align_input = mafft_trim if os.path.isfile(mafft_trim) else mafft_raw

        if analysis_mode.startswith("Single"):
            if not os.path.isfile(align_input):
                with st.spinner("Running MAFFT…"):
                    subprocess.run(
                        ["bash", MAFFT_SCRIPT, selected_sample, blast_set],
                        check=False,
                    )
                align_input = mafft_trim if os.path.isfile(mafft_trim) else mafft_raw
            if not os.path.isfile(align_input):
                st.error("No MAFFT alignment. Run MAFFT on the Consensuses page first.")
                return
            n = count_fasta_seqs(align_input)
            if n < 3:
                st.error(f"Need ≥3 sequences (found {n}).")
                return
        elif not os.path.isfile(align_input):
            st.error(f"Legacy alignment missing: {align_input}")
            return

        with st.status("Building tree…", expanded=True) as status:
            st.write("trimAl (-gappyout)" if align_input == mafft_raw else "Using trimmed alignment")
            try:
                run_iqtree_on_alignment(
                    align_input,
                    iq_prefix,
                    model=model,
                    bootstrap=bootstrap,
                    alrt=alrt,
                    trimmed_copy_path=trimmed,
                )
                status.update(label="IQ-TREE started in background", state="complete")
            except subprocess.CalledProcessError as exc:
                st.error(f"Pipeline failed: {exc}")
                return

        st.rerun()


def show_iqtree_page() -> None:
    page_hero(
        "Phylogenetic trees (IQ-TREE)",
        "Top-10 cross-barcode trees · per-sample alignments · ggtree views",
    )

    trees = find_tree_files()
    tab_explore, tab_build = st.tabs(["Tree explorer", "Generate"])

    with tab_explore:
        if not trees:
            st.info("No trees yet. Use **Generate** to run IQ-TREE.")
        else:
            render_explorer()

    with tab_build:
        render_generate_panel(has_trees=bool(trees))
