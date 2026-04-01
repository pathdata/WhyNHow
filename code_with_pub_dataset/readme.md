# TCGA BRCA Multimodal Integration

**RNA-seq + Clinical Data Analysis Pipeline**

Prepared by **Dr Priya Lakshmi Narayanan**  
Institute of Cancer Research  
GitHub: [pathdata/LearningCurve](https://github.com/pathdata/LearningCurve) ·[pathdata/WhyNHow](https://github.com/pathdata/WhyNHow) 

## Overview

| Modality | File | Content |
|---|---|---|
| Transcriptomic | `tcga_brca_counts_table.txt` | Raw RNA-seq counts — 60,488 genes × 32 patients |
| Molecular / Survival | `tcga_brca_clinical_RNASeq.txt` | PAM50 subtype, ER/PR/HER2, mutation count, OS |
| Clinical | `tcga_brca_clinical.txt` | Receptor status, age, OS from 1,098-patient cohort |

---

## Requirements

### Python version
Python 3.8 or higher

### Dependencies

Install all required packages with:

```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn lifelines umap-learn
```

| Package | Version tested | Purpose |
|---|---|---|
| pandas | ≥ 1.5 | Data loading and manipulation |
| numpy | ≥ 1.23 | Numerical operations |
| matplotlib | ≥ 3.6 | Plotting |
| seaborn | ≥ 0.12 | Statistical visualisation |
| scipy | ≥ 1.9 | Statistical tests |
| scikit-learn | ≥ 1.1 | PCA, Random Forest, Logistic Regression |
| lifelines | ≥ 0.27 | Kaplan–Meier survival analysis |
| umap-learn | ≥ 0.5 | UMAP dimensionality reduction (optional) |

> `umap-learn` and `lifelines` are optional — the script will skip those sections gracefully if they are not installed and print an informative message.
>
> ---

## Input Files

Place all three data files in the **same directory** as `TCGA_BRCA_Multimodal.py` before running:

```
your_directory/
├── TCGA_BRCA_Multimodal.py
├── tcga_brca_counts_table.txt
├── tcga_brca_clinical_RNASeq.txt
└── tcga_brca_clinical.txt
```

### File formats

**`tcga_brca_counts_table.txt`**  
Tab-separated. Rows = Ensembl gene IDs, columns = patient IDs (`TCGA.XX.XXXX.01` format). Values are raw integer read counts.

**`tcga_brca_clinical_RNASeq.txt`**  
Tab-separated. Rows = patient IDs (`TCGA.XX.XXXX.01`). Columns include PAM50 subtype, ER/PR/HER2 status, OS months, OS status, tumour stage, mutation count, fraction genome altered, and clustering variables.

**`tcga_brca_clinical.txt`**  
Tab-separated. Rows = patients from the full TCGA-BRCA cohort (n=1,098). Patient IDs use dash format (`TCGA-XX-XXXX`). Columns include ER, PR, HER2 receptor status, age, OS in days, and OS status.

---

## Running the Script

```bash
python TCGA_BRCA_Multimodal.py
```

The script prints progress for each of its 8 sections and saves figures automatically. No arguments or configuration are required.

---

## Pipeline: 8 Analysis Sections

### Section 1 — Patient ID Harmonisation
Loads all three files and reconciles mismatched ID formats:

- `TCGA.AR.A255.01` (counts + RNASeq files) → converted to `TCGA-AR-A255`
- Matched against `TCGA-AR-A255` style IDs in the clinical file
- Result: 32 patients matched across all three modalities

### Section 2 — RNA-seq Preprocessing
Transforms raw integer counts into analysis-ready expression values:

1. **CPM normalisation** — divides each sample's counts by its library size × 10⁶, correcting for sequencing depth differences between patients
2. **Log₂(CPM + 1) transformation** — stabilises variance and compresses the dynamic range; the `+1` pseudocount avoids log(0)
3. **Variance filter** — retains the top 5,000 most variable genes across patients, removing uninformative lowly-expressed genes

Output figure: `fig1_normalisation.png`

### Section 3 — PCA (Principal Component Analysis)
Reduces 5,000 gene dimensions to a compact representation:

- Genes standardised (zero mean, unit variance) before PCA
- 20 PCs computed; scree plot shows variance explained per component
- Scatter plots of PC1 vs PC2 and PC1 vs PC3, coloured by PAM50 subtype and ER status

Output figures: `fig2_pca_scree.png`, `fig3_pca_scatter.png`

### Section 4 — UMAP
Applies non-linear dimensionality reduction to capture structure that PCA misses. Patients are plotted in 2D UMAP space coloured by PAM50 subtype and by overall survival status. Requires `umap-learn`; skipped gracefully if not installed.

Output figure: `fig4_umap.png`

### Section 5 — Multimodal Feature Matrix
Assembles a unified 19-feature design matrix per patient from all three modalities:

| Feature group | Features | Count |
|---|---|---|
| Transcriptomic | PC1 – PC10 from RNA-seq PCA | 10 |
| Molecular | ER, PR, HER2 positivity; node positivity; metastasis; log mutation count; fraction genome altered | 7 |
| Clinical | Age; tumour stage (ordinal 1–4) | 2 |

A Spearman correlation heatmap visualises cross-modality relationships, with black lines marking modality boundaries.

Output figure: `fig5_correlation_heatmap.png`

### Section 6 — PAM50 Subtype Classification

**Why only 17 of 32 samples are used here**

All 32 samples are used in Sections 1–5 (preprocessing, PCA, multimodal feature building) and Section 8 (OS prediction), since those analyses require only RNA-seq counts and OS status, which are present for all 32 patients.

Classification requires a ground-truth label. In `tcga_brca_clinical_RNASeq.txt`, 15 of the 32 samples have `PAM50_SUBTYPE = NaN`. No proxy label is available for those 15 in any other column (`INTEGRATED_CLUSTERS_WITH_PAM50` and `SIGCLUST_INTRINSIC_MRNA` are also NaN). They are therefore excluded from classification only — not from the rest of the pipeline.

**Leave-One-Out cross-validation (LOO-CV) on 17 labelled samples**

A fixed train/test split (e.g. 80/20) would leave only ~3 test samples, which is too few for meaningful evaluation. LOO-CV avoids this problem:

```
Fold 1:  train on samples 2–17 (n=16),  test on sample 1
Fold 2:  train on samples 1,3–17 (n=16), test on sample 2
...
Fold 17: train on samples 1–16 (n=16),  test on sample 17
```

Every sample is the test sample exactly once. There is no permanently held-out set. The reported accuracy is the mean over all 17 folds. LOO-CV is the recommended strategy when n < 30.

**Statistical testing of classification accuracy**

For each modality, a **chi-square goodness-of-fit test** assesses whether the LOO-CV accuracy is significantly better than chance:

- H₀: the classifier performs at chance level (1/3 ≈ 0.333 for 3 classes)
- Observed: [number correct, number wrong]
- Expected: [n × 0.333, n × 0.667]

The figure reports χ², p-value, and significance level for each modality. Two reference lines are shown:

- Dashed line: uniform chance level (0.333)
- Dotted line: majority-class baseline (0.647) — the accuracy achieved by always predicting Luminal A, the most common class (11/17 samples). Any model at or below this threshold offers no improvement over a trivial classifier.

| Modality | LOO Accuracy | Correct/n | χ² | p | vs baseline |
|---|---|---|---|---|---|
| Transcriptomic (PC1–10) | 0.706 | 12/17 | 10.62 | 0.0011 ** | above |
| Molecular features | 0.706 | 12/17 | 10.62 | 0.0011 ** | above |
| Clinical features | 0.294 | 5/17 | 0.12 | 0.73 ns | below chance |
| Multimodal (all) | 0.647 | 11/17 | 7.53 | 0.0061 ** | at baseline |

Output figures: `fig6_classification_accuracy.png`, `fig7_feature_importance.png`

### Section 7 — Survival Analysis
Three complementary survival analyses:

**Kaplan–Meier by PAM50 subtype** — survival curves for Luminal A, Luminal B, and Basal-like patients with 95% confidence intervals and a multivariate log-rank test p-value.

**Kaplan–Meier by PC1 split** — patients are split at the PC1 median into PC1-High and PC1-Low groups; log-rank test assesses whether transcriptomic variation stratifies survival.

**PC1 gene loadings** — bar chart of the 15 genes with the highest and lowest PC1 loadings, identifying the transcriptomic drivers of the primary axis of variation. Requires `lifelines`.

Output figures: `fig8_km_pam50.png`, `fig9_km_pc1.png`, `fig10_pc1_loadings.png`

### Section 8 — Multimodal OS Prediction
Logistic regression (L2 regularisation, C=0.1) predicts overall survival status (LIVING vs DECEASED) from all 19 multimodal features:

- LOO-CV used throughout to avoid overfitting given the small sample (n=32)
- AUC-ROC reported as the primary performance metric
- ROC curve plotted alongside a coefficient plot showing which features and modalities drive prediction

Output figure: `fig11_os_prediction.png`

---

## Output Files

All figures are saved as PNG files in the working directory:

| File | Section | Description |
|---|---|---|
| `fig1_normalisation.png` | 2 | Raw counts vs log₂(CPM+1) distribution |
| `fig2_pca_scree.png` | 3 | Variance explained per PC |
| `fig3_pca_scatter.png` | 3 | PC1/2/3 coloured by PAM50 and ER status |
| `fig4_umap.png` | 4 | UMAP coloured by subtype and OS status |
| `fig5_correlation_heatmap.png` | 5 | Spearman correlation across all features |
| `fig6_classification_accuracy.png` | 6 | LOO-CV accuracy by modality |
| `fig7_feature_importance.png` | 6 | Random Forest feature importances |
| `fig8_km_pam50.png` | 7 | Kaplan–Meier by PAM50 subtype |
| `fig9_km_pc1.png` | 7 | Kaplan–Meier by PC1 split |
| `fig10_pc1_loadings.png` | 7 | Top genes driving PC1 |
| `fig11_os_prediction.png` | 8 | ROC curve and logistic regression coefficients |

---

## Sample Usage Across the Pipeline

| Section | Samples used | Reason |
|---|---|---|
| 1–5 (preprocessing, PCA, feature matrix) | All **32** | No labels required |
| 6 (PAM50 classification) | **17** with PAM50 label | 15 have `PAM50_SUBTYPE = NaN` in source; no proxy available |
| 7 (survival KM) | All **32** with `OS_MONTHS` | PAM50 label not required for KM curves |
| 8 (OS prediction) | All **32** | `OS_STATUS` available for all patients |

The 15 unlabelled samples are not wasted — they contribute to the PCA embedding, the multimodal feature matrix, and the OS prediction model. They are only excluded from the PAM50 classification task where a ground-truth label is required.

---

## Interpreting the Results

**PCA and UMAP** — Basal-like tumours typically separate from Luminal subtypes along PC1. Strong separation indicates that gene expression alone carries subtype signal.

**LOO-CV accuracy** — With n=32, individual LOO folds have high variance. Treat accuracy values as indicative rather than definitive. A larger cohort (full TCGA-BRCA n=1,098) is needed for robust benchmarking.

**PC1 as a survival biomarker** — A significant log-rank p-value for the PC1 split suggests the dominant axis of transcriptomic variation is clinically relevant, consistent with known associations between PAM50 subtype and prognosis.

**OS prediction AUC** — The cohort has only 6 deceased patients out of 32, making this a class-imbalanced problem. AUC is more informative than accuracy in this setting, but results should be interpreted cautiously.

---

## Limitations and Next Steps

| Limitation | Recommendation |
|---|---|
| Small sample (n=32) | Validate on full TCGA-BRCA cohort (n=1,098) available from GDC portal |
| No batch correction | Apply ComBat or similar before cross-cohort comparisons |
| Gene IDs are Ensembl | Map to gene symbols using `pyensembl` or a GTF annotation file for biological interpretation |
| No imaging modality | Extend with H&E or multiplex IHC spatial features for true multimodal fusion |
| Logistic regression for OS | Consider Cox proportional hazards model for time-to-event outcomes |

---

## Data Source

Data are derived from **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** project.

- Data portal: https://portal.gdc.cancer.gov/projects/TCGA-BRCA
- Publication: Cancer Genome Atlas Network. *Nature* 490, 61–70 (2012)
- Access: Open-access tier (no dbGaP approval required for the files used here)

---

## Related Resources

-  Narayanan et al., *Nature* 2021 — Deep learning to unmask immune microecology in breast cancer
- [pathdata/LearningCurve](https://github.com/pathdata/LearningCurve) — COMP0188 Deep Representation Learning teaching notebooks (UCL)
- [pathdata/WhyNHow](https://github.com/pathdata/WhyNHow) — Single-cell RNA-seq and TCR analytic workflow
- [pathdata/https://github.com/pathdata/WhyNHow/tree/master/code_with_pub_dataset] Toy dataset for data harmonisation and integration
- [lifelines documentation](https://lifelines.readthedocs.io/)
- [UMAP documentation](https://umap-learn.readthedocs.io/)

---


