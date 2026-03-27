# `sc_output` — UMAP visualization scripts

This folder typically holds **comprehensive pipeline outputs** (`*_comprehensive*.h5ad`) from `06_comprehensive_pipeline_standalone_PN_version.py`, plus figures and tables produced by the helper scripts below.

All three scripts live in the parent **`code/`** directory. Run them from that directory (or pass absolute paths). They read `.h5ad` files from `--input-dir` (defaults to `code/sc_output` next to each script) and write PNG/CSV into `--output-dir`.

**Requirements**

- Files named like `*_comprehensive*.h5ad`.
- `Patient` in `adata.obs` (from the comprehensive pipeline metadata step).
- For clustering-based markers (`06_c`), `seurat_clusters` must be present.

**Important:** Each script **merges all GEX samples per patient** (inner join on genes) and **recomputes PCA / neighbors / UMAP** on that merged object. Embeddings are **not comparable** between patients (different panels). Merging several per-sample HVG subsets can **reduce** the number of shared genes; some requested genes may be absent from the matrix.

---

## `06_b_visualize_umap_by_patientID.py`

**Purpose:** One **side-by-side** figure: UMAP for each patient, colored by **metadata** (default: sample ID).

**Typical flow**

1. Load all `*_comprehensive*.h5ad` files and group by `Patient`.
2. For each requested patient, concatenate samples and recompute embedding.
3. Save a single PNG with one panel per patient (left-to-right order = `--patients`).

**Default outputs** (under `--output-dir`)

- `umap_by_patient_P01_P03.png` (or `--out-name`)

**Useful options**

| Option | Default | Notes |
|--------|---------|--------|
| `--input-dir` | `WhyNHow` | Directory containing `.h5ad` files |
| `--output-dir` | `sc_output` | Where to save the PNG |
| `--patients` | `P01 P03` | Panel order |
| `--color-by` | `gex_sample` | `obs` column for point colors (e.g. `orig.ident`, `seurat_clusters` if present) |
| `--out-name` | `umap_by_patient_P01_P03.png` | Output filename |
| `--seed` | `0` | Random seed for PCA/UMAP |
| `--full-reprocess` | off | Full normalize → HVG → scale → PCA (only if `X` is raw-count-like) |

**Example**

```bash
conda run -n trans python 06_b_visualize_umap_by_patientID.py --input-dir WhyNHow --output-dir sc_output
```

---

## `06_c_visualise_umap_by_mgenes.py`

**Purpose:** Find **marker genes** that appear in **both** patients (data-driven), then plot those genes on each patient’s UMAP in a **grid** (rows = patients, columns = genes).

**Definition of “common markers”**

1. For each patient, run `sc.tl.rank_genes_groups` on **`--groupby`** (default: `seurat_clusters`).
2. For each patient, build a set of genes: **union** of the top **`--top-n`** genes per cluster.
3. **Common markers** = **intersection** of those sets, restricted to genes present in **each** patient’s merged matrix.

**Default outputs**

- `{out_prefix}_genes.csv` — full sorted list of common genes (default prefix: `common_markers_P01_P03`).
- `{out_prefix}_umap_grid.png` — UMAP grid; only the first **`--max-genes`** genes are plotted (alphabetically).

If the intersection is empty, the script may write `{out_prefix}_per_patient_union_only.csv` for debugging.

**Useful options**

| Option | Default | Notes |
|--------|---------|--------|
| `--top-n` | `10` | Top genes per cluster per patient, before intersection |
| `--max-genes` | `16` | Cap on genes plotted |
| `--groupby` | `seurat_clusters` | Cluster column for DE |
| `--method` | `wilcoxon` | `rank_genes_groups` method |
| `--out-prefix` | `common_markers_P01_P03` | Prefix for CSV/PNG filenames |
| `--seed`, `--full-reprocess` | | Same meaning as in `06_b` |

**Example**

```bash
conda run -n trans python 06_c_visualise_umap_by_mgenes.py --input-dir WhyNHow --output-dir sc_output --top-n 15
```

---

## `06_c1_visualise_umap_by_genelist.py`

**Purpose:** Plot a **fixed list** of known marker genes (no differential expression). Same per-patient merge + UMAP as `06_b`, then **UMAP colored by each gene** in a grid.

**Default gene list**

CD3D, CD3E, CD4, CD8A, CD8B, GZMA, GZMB, PRF1, NKG7, GNLY, PDCD1, HAVCR2, LAG3, TIGIT, MKI67, PCNA, TOP2A, CCR7, SELL, IL7R, IFNG, TNF, FOXP3, TBX21

**Default outputs**

- `{out_prefix}_gene_availability.csv` — whether each gene is present in each patient’s merged matrix.
- `{out_prefix}_umap_grid.png` — multi-row layout: each patient gets a block of rows; **`--max-cols`** panels per row (default `8`).

**Useful options**

| Option | Default | Notes |
|--------|---------|--------|
| `--genes-file` | (none) | Optional file: one gene per line or comma-separated; `#` starts comments |
| `--max-cols` | `8` | Panels per row within each patient block |
| `--out-prefix` | `umap_known_markers_P01_P03` | Prefix for CSV/PNG |
| `--patients`, `--seed`, `--full-reprocess` | | Same idea as above |

**Example**

```bash
conda run -n trans python 06_c1_visualise_umap_by_genelist.py --input-dir WhyNHow --output-dir sc_output
```

---

## Relationship between the scripts

| Script | Question it answers |
|--------|---------------------|
| `06_b` | How are samples (or clusters) arranged in UMAP **within** each patient? |
| `06_c` | Which **DE-derived** marker genes are **shared** between patients, and how do they look on UMAP? |
| `06_c1` | How does a **predefined** gene panel look on UMAP **per patient**? |

Scripts `06_c` and `06_c1` import embedding helpers from `06_b` (`06_b_visualize_umap_by_patientID.py`) so behaviour stays consistent.
