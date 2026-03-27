"""
Common marker genes across patients + UMAP visualization.

Workflow
--------
1. Load comprehensive pipeline outputs (*_comprehensive*.h5ad) and group by ``Patient``.
2. For each target patient, merge all GEX samples (inner join on genes) and recompute
   PCA / neighbors / UMAP (same defaults as ``06_b_visualize_umap_by_patientID.py``).
3. Run ``sc.tl.rank_genes_groups`` on ``seurat_clusters`` (Leiden) within that patient.
4. Build each patient's marker set as the **union** of the top ``--top-n`` genes per
   cluster (union across clusters).
5. **Common markers** = intersection of those per-patient sets, restricted to genes
   present in every merged object (so they can be plotted on both sides).
6. Save a CSV listing those genes and a multi-panel PNG: one row per patient, one
   column per gene (UMAP colored by expression).

Requires ``Patient`` and ``seurat_clusters`` in obs (from the comprehensive pipeline).

Example
-------
  conda run -n trans python 06_c_visualise_umap_by_mgenes.py \\
    --input-dir sc_output --output-dir sc_output --patients P01 P03
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

_CODE_DIR = Path(__file__).resolve().parent
_DEFAULT_INPUT = _CODE_DIR / "sc_output"
_DEFAULT_OUTPUT = _CODE_DIR / "sc_output"


def _load_umap_by_patient_module():
    """Load numeric-prefixed sibling module without packaging."""
    path = _CODE_DIR / "06_b_visualize_umap_by_patientID.py"
    spec = importlib.util.spec_from_file_location("umap_by_patient_b", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _markers_union_top_per_cluster(
    adata: sc.AnnData,
    groupby: str,
    top_n: int,
    method: str = "wilcoxon",
) -> Set[str]:
    """
    Marker set = union over clusters of the top ``top_n`` genes by score
    (Wilcoxon rank-sum by default).
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"Missing '{groupby}' in obs")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=False,
    )
    df = sc.get.rank_genes_groups_df(adata, group=None)
    if df is None or len(df) == 0:
        return set()
    markers: Set[str] = set()
    for grp in df["group"].unique():
        sub = df[df["group"] == grp].head(int(top_n))
        markers.update(sub["names"].astype(str).tolist())
    return markers


def _prepare_merged_per_patient(
    by_patient: Dict[str, List[Tuple[str, sc.AnnData]]],
    patients: List[str],
    merge_fn,
    embed_fn,
    random_state: int,
    full_reprocess: bool,
) -> Dict[str, sc.AnnData]:
    out: Dict[str, sc.AnnData] = {}
    for pid in patients:
        if pid not in by_patient or not by_patient[pid]:
            print(f"[WARNING] No samples for patient {pid}, skipping.")
            continue
        print(f"\nMerging + UMAP for {pid} ({len(by_patient[pid])} sample(s))...")
        merged = merge_fn(by_patient[pid])
        print(f"  Shape: {merged.n_obs:,} cells x {merged.n_vars:,} genes")
        embed_fn(merged, random_state=random_state, full_reprocess=full_reprocess)
        out[pid] = merged
    return out


def _intersect_marker_sets(
    per_patient_markers: Dict[str, Set[str]],
    var_names_per_patient: Dict[str, List[str]],
) -> List[str]:
    """Intersection of marker sets, keeping only genes present in all patients' matrices."""
    if len(per_patient_markers) < 2:
        return []
    sets = list(per_patient_markers.values())
    common = set.intersection(*sets) if sets else set()
    vars_intersection = set(var_names_per_patient[list(var_names_per_patient.keys())[0]])
    for v in var_names_per_patient.values():
        vars_intersection &= set(v)
    common &= vars_intersection
    return sorted(common)


def plot_umap_marker_grid(
    adata_by_patient: Dict[str, sc.AnnData],
    patients_order: List[str],
    genes: List[str],
    output_path: Path,
) -> None:
    """One row per patient, one column per gene; UMAP colored by gene expression."""
    patients = [p for p in patients_order if p in adata_by_patient]
    if not patients or not genes:
        print("[ERROR] Nothing to plot (missing patients or genes).")
        return

    nrows = len(patients)
    ncols = len(genes)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.2 * nrows),
        squeeze=False,
    )

    for i, pid in enumerate(patients):
        adata = adata_by_patient[pid]
        for j, gene in enumerate(genes):
            ax = axes[i][j]
            if gene not in adata.var_names:
                ax.text(0.5, 0.5, f"{gene}\nnot in matrix", ha="center", va="center")
                ax.set_axis_off()
                continue
            sc.pl.umap(
                adata,
                color=gene,
                ax=ax,
                show=False,
                frameon=False,
                cmap="viridis",
                title=f"{pid} — {gene}",
                legend_loc="right margin",
            )

    plt.suptitle(
        "Common cluster markers (intersection across patients)",
        fontsize=14,
        y=1.01,
    )
    plt.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n[OK] Saved figure: {output_path}")


def main() -> None:
    mod = _load_umap_by_patient_module()
    load_h5ad_grouped_by_patient = mod.load_h5ad_grouped_by_patient
    merge_patient_samples = mod.merge_patient_samples
    recompute_umap_merged = mod.recompute_umap_merged

    parser = argparse.ArgumentParser(
        description="Intersect per-patient cluster markers and plot UMAPs for common genes."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=_DEFAULT_INPUT,
        help=f"Directory with *_comprehensive*.h5ad (default: {_DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTPUT,
        help=f"Directory for CSV/PNG (default: {_DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        default=["P01", "P03"],
        help="Patients to compare (default: P01 P03)",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        default="seurat_clusters",
        help="obs column for DE (default: seurat_clusters)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Top N genes per cluster included in each patient's marker set (default: 10)",
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=16,
        help="Max common genes to plot (default: 16)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="wilcoxon",
        help="rank_genes_groups method (default: wilcoxon)",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="common_markers_P01_P03",
        help="Prefix for output CSV/PNG names",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--full-reprocess",
        action="store_true",
        help="Passed to 06_b embedding (raw-count-like matrices only)",
    )

    args = parser.parse_args()

    sc.settings.verbosity = 1
    sc.set_figure_params(dpi=150, facecolor="white")

    print("=" * 80)
    print("COMMON MARKER GENES + UMAP (by patient)")
    print("=" * 80)
    print(f"Input directory:  {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Patients:         {args.patients}")
    print(f"groupby:          {args.groupby}")
    print(f"top N / cluster:  {args.top_n}")
    print()

    by_patient = load_h5ad_grouped_by_patient(args.input_dir)
    if not by_patient:
        print("[ERROR] No data loaded.")
        return

    adata_by_patient = _prepare_merged_per_patient(
        by_patient,
        list(args.patients),
        merge_patient_samples,
        recompute_umap_merged,
        random_state=args.seed,
        full_reprocess=args.full_reprocess,
    )

    if len(adata_by_patient) < 2:
        print("[ERROR] Need at least two successfully merged patients.")
        return

    per_patient_markers: Dict[str, Set[str]] = {}
    var_names_per_patient: Dict[str, List[str]] = {}

    for pid, adata in adata_by_patient.items():
        if args.groupby not in adata.obs.columns:
            print(f"[ERROR] '{args.groupby}' missing for patient {pid}.")
            return
        print(f"\nRanked genes for {pid} ({args.groupby})...")
        try:
            mset = _markers_union_top_per_cluster(
                adata,
                groupby=args.groupby,
                top_n=args.top_n,
                method=args.method,
            )
        except Exception as e:
            print(f"[ERROR] rank_genes_groups failed for {pid}: {e}")
            raise
        print(f"  Union of top-{args.top_n} per cluster: {len(mset)} genes")
        per_patient_markers[pid] = mset
        var_names_per_patient[pid] = list(adata.var_names)

    common = _intersect_marker_sets(per_patient_markers, var_names_per_patient)
    print(f"\nCommon markers (intersection): {len(common)} genes")

    if not common:
        print(
            "[WARNING] Empty intersection. Try increasing --top-n, or check that "
            "cluster markers overlap and that gene symbols match across merged objects."
        )
        # Still save per-patient unions for debugging
        dbg = pd.DataFrame(
            {f"markers_{k}": pd.Series(sorted(v)) for k, v in per_patient_markers.items()}
        )
        dbg_path = args.output_dir / f"{args.out_prefix}_per_patient_union_only.csv"
        dbg.to_csv(dbg_path, index=False)
        print(f"  Wrote per-patient union columns for inspection: {dbg_path}")
        return

    common_limited = common[: args.max_genes]
    if len(common) > args.max_genes:
        print(
            f"  Plotting first {args.max_genes} genes (alphabetically); "
            f"{len(common) - args.max_genes} not shown."
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{args.out_prefix}_genes.csv"
    pd.DataFrame({"gene": common}).to_csv(csv_path, index=False)
    print(f"[OK] Wrote gene list: {csv_path}")

    png_path = out_dir / f"{args.out_prefix}_umap_grid.png"
    plot_umap_marker_grid(
        adata_by_patient,
        list(args.patients),
        common_limited,
        png_path,
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
