"""
UMAP by patient ID — side-by-side panels (one embedding per patient).

Loads comprehensive pipeline outputs (*_comprehensive*.h5ad), merges all GEX samples
that share the same Patient (from adata.obs['Patient']), recomputes PCA/neighbors/UMAP
on each merged object (coordinates are not comparable across patients), and saves one figure.

Expects Patient metadata from 06_comprehensive_pipeline_standalone_PN_version.py (Step 16).
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

# Defaults relative to this file
_CODE_DIR = Path(__file__).resolve().parent
_DEFAULT_INPUT = _CODE_DIR / "sc_output"
_DEFAULT_OUTPUT = _CODE_DIR / "sc_output"


def _parse_sample_id(h5ad_path: Path) -> str:
    name = h5ad_path.name
    if "_comprehensive" in name:
        return name.split("_comprehensive")[0]
    return h5ad_path.stem


def load_h5ad_grouped_by_patient(
    input_dir: Path,
) -> Dict[str, List[Tuple[str, sc.AnnData]]]:
    """
    Load all *_comprehensive*.h5ad files and group by obs['Patient'].
    Returns: { patient_id: [(gex_sample_id, adata), ...] }
    """
    input_dir = Path(input_dir)
    by_patient: Dict[str, List[Tuple[str, sc.AnnData]]] = defaultdict(list)

    h5ad_files = sorted(input_dir.glob("*_comprehensive*.h5ad"))
    if not h5ad_files:
        print(f"[WARNING] No *_comprehensive*.h5ad files in {input_dir}")

    for h5ad_path in h5ad_files:
        sample_id = _parse_sample_id(h5ad_path)
        adata = sc.read_h5ad(h5ad_path)
        if "Patient" not in adata.obs.columns:
            print(f"  [SKIP] {h5ad_path.name}: no 'Patient' in obs")
            continue
        patient = str(adata.obs["Patient"].iloc[0])
        if patient in ("", "nan", "None"):
            print(f"  [SKIP] {h5ad_path.name}: empty Patient")
            continue
        by_patient[patient].append((sample_id, adata))
        print(f"  Loaded {sample_id} -> Patient={patient} ({adata.n_obs:,} cells)")

    return dict(by_patient)


def recompute_umap_merged(
    adata: sc.AnnData, random_state: int = 0, full_reprocess: bool = False
) -> None:
    """
    Build a joint embedding on merged cells.

    Default (`full_reprocess=False`): comprehensive pipeline `.h5ad` files already hold
    processed values (log-normalized, HVG-subset, scaled). Re-applying normalize_total +
    seurat_v3 HVG on that matrix is invalid. We run PCA -> neighbors -> UMAP on the
    merged `X` (inner-joined genes across samples).

    Optional `full_reprocess=True`: same recipe as `07_visualize_gene_expression`
    combined plot (for raw-count-like matrices only).
    """
    if full_reprocess:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, n_comps=50, random_state=random_state)
        adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
        n_pcs = 50
    else:
        n_pcs = int(min(50, adata.n_vars, max(adata.n_obs - 1, 1)))
        n_pcs = max(1, n_pcs)
        sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)

    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_pcs, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)


def merge_patient_samples(
    samples: List[Tuple[str, sc.AnnData]],
) -> sc.AnnData:
    """Concatenate all GEX samples for one patient (inner join on genes)."""
    adata_list: List[sc.AnnData] = []
    keys: List[str] = []
    for gex_id, ad in samples:
        ac = ad.copy()
        ac.obs["gex_sample"] = gex_id
        adata_list.append(ac)
        keys.append(gex_id)

    adata_merged = sc.concat(
        adata_list,
        join="inner",
        label="batch",
        keys=keys,
        index_unique="_",
    )
    return adata_merged


def plot_umap_panels_by_patient(
    by_patient: Dict[str, List[Tuple[str, sc.AnnData]]],
    patients_order: List[str],
    output_path: Path,
    color_by: str = "gex_sample",
    random_state: int = 0,
    full_reprocess: bool = False,
    figsize_per_panel: Tuple[float, float] = (6.5, 5.5),
) -> None:
    """
    One UMAP per patient in a single row of axes; saves PNG.
    color_by must exist after concat (e.g. 'gex_sample', 'orig.ident', 'seurat_clusters').
    """
    panels: List[Tuple[str, sc.AnnData]] = []

    for pid in patients_order:
        if pid not in by_patient or not by_patient[pid]:
            print(f"[WARNING] No samples for patient {pid}, skipping panel.")
            continue
        print(f"\nMerging and embedding Patient {pid} ({len(by_patient[pid])} sample(s))...")
        merged = merge_patient_samples(by_patient[pid])
        print(f"  Merged shape: {merged.n_obs:,} cells x {merged.n_vars:,} genes")

        recompute_umap_merged(
            merged, random_state=random_state, full_reprocess=full_reprocess
        )

        if color_by not in merged.obs.columns:
            print(
                f"  [WARNING] '{color_by}' not in obs; falling back to 'orig.ident' if present."
            )
            fallback = "orig.ident" if "orig.ident" in merged.obs.columns else None
            if fallback is None:
                print("  [ERROR] Cannot color UMAP — no suitable column.")
                continue
            use_color = fallback
        else:
            use_color = color_by

        panels.append((pid, merged, use_color))

    if not panels:
        print("[ERROR] Nothing to plot.")
        return

    n = len(panels)
    fig, axes = plt.subplots(
        1,
        n,
        figsize=(figsize_per_panel[0] * n, figsize_per_panel[1]),
        squeeze=False,
    )
    ax_row = axes[0]

    for i, (pid, adata_m, use_color) in enumerate(panels):
        ax = ax_row[i]
        sc.pl.umap(
            adata_m,
            color=use_color,
            ax=ax,
            show=False,
            title=f"Patient {pid}",
            frameon=False,
            legend_loc="right margin",
        )

    plt.suptitle("UMAP by patient (separate embedding per patient)", fontsize=14, y=1.02)
    plt.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n[OK] Saved: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot side-by-side UMAPs per Patient from comprehensive .h5ad files."
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
        help=f"Directory for PNG output (default: {_DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        default=["P01", "P03"],
        help="Patient IDs to plot, in left-to-right order (default: P01 P03)",
    )
    parser.add_argument(
        "--color-by",
        type=str,
        default="gex_sample",
        help="obs column for point color (default: gex_sample)",
    )
    parser.add_argument(
        "--out-name",
        type=str,
        default="umap_by_patient_P01_P03.png",
        help="Output PNG filename (under --output-dir)",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed for PCA/UMAP")
    parser.add_argument(
        "--full-reprocess",
        action="store_true",
        help="Run normalize+log1p+seurat_v3 HVG+scale before PCA (only for raw-count-like X)",
    )

    args = parser.parse_args()

    sc.settings.verbosity = 1
    sc.set_figure_params(dpi=150, facecolor="white")

    print("=" * 80)
    print("UMAP BY PATIENT ID")
    print("=" * 80)
    print(f"Input directory:  {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Patients (order): {args.patients}")
    print(f"Color by:         {args.color_by}")
    print()

    by_patient = load_h5ad_grouped_by_patient(args.input_dir)
    if not by_patient:
        print("[ERROR] No data loaded. Check --input-dir and file names.")
        return

    print(f"\nFound patients in data: {sorted(by_patient.keys())}")

    out_file = args.output_dir / args.out_name
    plot_umap_panels_by_patient(
        by_patient=by_patient,
        patients_order=list(args.patients),
        output_path=out_file,
        color_by=args.color_by,
        random_state=args.seed,
        full_reprocess=args.full_reprocess,
    )


if __name__ == "__main__":
    main()
