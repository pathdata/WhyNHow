"""
UMAP visualization for a fixed list of known marker genes (per patient).

Loads comprehensive ``*_comprehensive*.h5ad`` files, merges by ``Patient`` (same as
``06_b_visualize_umap_by_patientID.py``), recomputes PCA/neighbors/UMAP per merged
object, then plots each gene on the UMAP in a grid (rows = patients, columns
chunked with ``--max-cols``).

Default gene list matches ``07_visualize_gene_expression.py`` T cell / cytotoxic /
exhaustion / proliferation / memory / cytokine / TF panel.

Example
-------
  conda run -n trans python 06_c1_visualise_umap_by_genelist.py \\
    --input-dir sc_output --output-dir sc_output --patients P01 P03
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

_CODE_DIR = Path(__file__).resolve().parent
_DEFAULT_INPUT = _CODE_DIR / "sc_output"
_DEFAULT_OUTPUT = _CODE_DIR / "sc_output"

# Known marker panel (order preserved for plotting)
DEFAULT_GENELIST: List[str] = [
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "CD8B",
    "GZMA",
    "GZMB",
    "PRF1",
    "NKG7",
    "GNLY",
    "PDCD1",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "MKI67",
    "PCNA",
    "TOP2A",
    "CCR7",
    "SELL",
    "IL7R",
    "IFNG",
    "TNF",
    "FOXP3",
    "TBX21",
]


def _load_umap_by_patient_module():
    path = _CODE_DIR / "06_b_visualize_umap_by_patientID.py"
    spec = importlib.util.spec_from_file_location("umap_by_patient_b", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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


def _parse_genelist_file(path: Path) -> List[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    genes: List[str] = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for part in line.replace(",", " ").split():
            g = part.strip().strip("'\"")
            if g:
                genes.append(g)
    return genes


def plot_umap_genelist_grid(
    adata_by_patient: Dict[str, sc.AnnData],
    patients_order: List[str],
    genes: List[str],
    output_path: Path,
    max_cols: int = 8,
) -> None:
    """
    For each patient, plot UMAPs in a sub-grid of ``max_cols`` columns; multiple
    rows per patient if needed.
    """
    patients = [p for p in patients_order if p in adata_by_patient]
    if not patients:
        print("[ERROR] No patients to plot.")
        return
    if not genes:
        print("[ERROR] No genes to plot.")
        return

    max_cols = max(1, int(max_cols))
    n_genes = len(genes)
    n_rows_per_patient = int(np.ceil(n_genes / max_cols))
    nrows = len(patients) * n_rows_per_patient

    fig, axes = plt.subplots(
        nrows,
        max_cols,
        figsize=(3.6 * max_cols, 3.6 * nrows),
        squeeze=False,
    )

    for p_idx, pid in enumerate(patients):
        adata = adata_by_patient[pid]
        base_row = p_idx * n_rows_per_patient
        for g_idx, gene in enumerate(genes):
            r = g_idx // max_cols
            c = g_idx % max_cols
            ax = axes[base_row + r][c]
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
        # hide unused axes in this patient's block
        used = n_genes
        for extra in range(used, n_rows_per_patient * max_cols):
            r = extra // max_cols
            cc = extra % max_cols
            axes[base_row + r][cc].axis("off")

    plt.suptitle(
        "Known marker genes (fixed gene list)",
        fontsize=14,
        y=1.005,
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
        description="Plot fixed gene list on per-patient UMAPs (from comprehensive h5ad)."
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
        help=f"Output directory (default: {_DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        default=["P01", "P03"],
        help="Patients to plot (default: P01 P03)",
    )
    parser.add_argument(
        "--genes-file",
        type=Path,
        default=None,
        help="Optional text file: one gene per line (or comma-separated). "
        "If omitted, uses built-in DEFAULT_GENELIST.",
    )
    parser.add_argument(
        "--max-cols",
        type=int,
        default=8,
        help="Maximum UMAP panels per row within each patient block (default: 8)",
    )
    parser.add_argument(
        "--out-prefix",
        type=str,
        default="umap_known_markers_P01_P03",
        help="Prefix for output PNG/CSV (default: umap_known_markers_P01_P03)",
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

    if args.genes_file is not None:
        genes = _parse_genelist_file(Path(args.genes_file))
        if not genes:
            print(f"[ERROR] No genes read from {args.genes_file}")
            return
    else:
        genes = list(DEFAULT_GENELIST)

    print("=" * 80)
    print("UMAP BY FIXED GENE LIST")
    print("=" * 80)
    print(f"Input directory:  {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Patients:         {args.patients}")
    print(f"Genes ({len(genes)}): {', '.join(genes)}")
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

    if not adata_by_patient:
        print("[ERROR] No merged patients to plot.")
        return

    # Availability table
    rows = []
    for g in genes:
        row = {"gene": g}
        for pid, adata in adata_by_patient.items():
            row[f"in_{pid}"] = g in adata.var_names
        rows.append(row)
    avail = pd.DataFrame(rows)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{args.out_prefix}_gene_availability.csv"
    avail.to_csv(csv_path, index=False)
    print(f"[OK] Wrote availability: {csv_path}")

    png_path = out_dir / f"{args.out_prefix}_umap_grid.png"
    plot_umap_genelist_grid(
        adata_by_patient,
        list(args.patients),
        genes,
        png_path,
        max_cols=args.max_cols,
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
