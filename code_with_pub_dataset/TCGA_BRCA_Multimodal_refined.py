"""
TCGA BRCA Multimodal Integration
=================================
RNA-seq + Clinical Data (matched by patient ID)

Run with:   python TCGA_BRCA_Multimodal_refined.py
Requires:   tcga_brca_counts_table.txt
            tcga_brca_clinical_RNASeq.txt
            tcga_brca_clinical.txt
            (all three files in the same directory as this script)

Install dependencies:
    pip install pandas numpy matplotlib seaborn scipy scikit-learn lifelines umap-learn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import chisquare
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_score, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

sns.set_theme(style="whitegrid", palette="Set2")
plt.rcParams.update({'figure.dpi': 120})

print("=" * 60)
print("TCGA BRCA Multimodal Analysis")
print("=" * 60)

# ──────────────────────────────────────────────────────────────
# SECTION 1 – Load & Harmonise Patient IDs
# ──────────────────────────────────────────────────────────────
print("\n[1/8] Loading data and harmonising patient IDs...")

counts_raw = pd.read_csv('tcga_brca_counts_table.txt',      sep='\t', index_col=0)
clin_rna   = pd.read_csv('tcga_brca_clinical_RNASeq.txt',  sep='\t', index_col=0)
clin       = pd.read_csv('tcga_brca_clinical.txt',          sep='\t', index_col=0)

print(f"  RNA-seq counts :  {counts_raw.shape[0]:,} genes x {counts_raw.shape[1]} patients")
print(f"  Clinical RNASeq: {clin_rna.shape}")
print(f"  Clinical:        {clin.shape}")

# Convert TCGA.AR.A255.01  →  TCGA-AR-A255
def dot_to_dash(pid):
    return '-'.join(pid.split('.')[:3])

clin_rna['barcode_short'] = [dot_to_dash(x) for x in clin_rna.index]
clin['barcode_upper']     = clin['barcode'].str.upper()

meta = clin_rna.merge(
    clin[['barcode_upper', 'age', 'ER', 'PR', 'HER2', 'OS_years', 'status']],
    left_on='barcode_short', right_on='barcode_upper', how='left'
)
meta.index = clin_rna.index

shared_ids = list(set(counts_raw.columns) & set(meta.index))
meta   = meta.loc[shared_ids].copy()
counts = counts_raw[shared_ids].copy()

print(f"  Matched patients across all 3 files: {len(shared_ids)}")

# ──────────────────────────────────────────────────────────────
# SECTION 2 – RNA-seq Pre-processing
# ──────────────────────────────────────────────────────────────
print("\n[2/8] Pre-processing RNA-seq counts...")

# CPM normalisation
lib_sizes = counts.sum(axis=0)
cpm       = counts.div(lib_sizes, axis=1) * 1e6

# Log2(CPM + 1)
log_cpm = np.log2(cpm + 1)

# Variance filter: top 5,000 genes
gene_var       = log_cpm.var(axis=1).sort_values(ascending=False)
top_genes      = gene_var.index[:5000]
log_cpm_filt   = log_cpm.loc[top_genes]
print(f"  Genes retained after variance filter: {len(top_genes):,} / {log_cpm.shape[0]:,}")

# Plot distribution before vs after
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
sid = counts.columns[0]
axes[0].hist(counts[sid], bins=100, color='#E07B54', log=True)
axes[0].set(title=f'Raw counts – {sid}', xlabel='Counts', ylabel='Frequency (log)')
axes[1].hist(log_cpm[sid], bins=60, color='#54A0E0')
axes[1].set(title=f'Log2(CPM+1) – {sid}', xlabel='Log2(CPM+1)', ylabel='Frequency')
plt.suptitle('RNA-seq: Before vs After Normalisation', fontsize=12)
plt.tight_layout()
plt.savefig('fig1_normalisation.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Saved: fig1_normalisation.png")

# ──────────────────────────────────────────────────────────────
# SECTION 3 – PCA
# ──────────────────────────────────────────────────────────────
print("\n[3/8] Running PCA on filtered log-CPM matrix...")

X_rna    = log_cpm_filt.T          # (32 patients x 5000 genes)
scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X_rna)

pca       = PCA(n_components=20, random_state=42)
pca_coords = pca.fit_transform(X_scaled)
pca_df    = pd.DataFrame(pca_coords, index=X_rna.index,
                          columns=[f'PC{i+1}' for i in range(20)])
var_exp   = pca.explained_variance_ratio_ * 100
print(f"  PCs 1–5 explain {var_exp[:5].sum():.1f}% of variance")

# Scree plot
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(range(1, 21), var_exp, color='#54A0E0', edgecolor='white')
ax.plot(range(1, 21), np.cumsum(var_exp), 'r-o', markersize=4, label='Cumulative')
ax.axhline(80, color='grey', linestyle='--', linewidth=0.8, label='80% threshold')
ax.set(xlabel='Principal Component', ylabel='Variance Explained (%)',
       title='PCA Scree Plot – RNA-seq (top 5,000 variable genes)')
ax.legend()
plt.tight_layout()
plt.savefig('fig2_pca_scree.png', dpi=150, bbox_inches='tight')
plt.show()

# PCA scatter coloured by PAM50 / ER status
pam50 = meta.loc[pca_df.index, 'PAM50_SUBTYPE'].fillna('Unknown')
pam50_palette = {
    'Luminal A': '#2ecc71', 'Luminal B': '#3498db',
    'Basal-like': '#e74c3c', 'HER2-enriched': '#9b59b6',
    'Normal-like': '#f39c12', 'Unknown': '#aaaaaa'
}

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, (xpc, ypc) in zip(axes, [('PC1', 'PC2'), ('PC1', 'PC3')]):
    for subtype, grp in pca_df.groupby(pam50):
        color = pam50_palette.get(subtype, '#aaaaaa')
        ax.scatter(grp[xpc], grp[ypc], label=subtype, color=color,
                   s=80, edgecolors='white', linewidth=0.5)
    xi, yi = int(xpc[2:]) - 1, int(ypc[2:]) - 1
    ax.set(xlabel=f'{xpc} ({var_exp[xi]:.1f}%)',
           ylabel=f'{ypc} ({var_exp[yi]:.1f}%)',
           title='PCA – PAM50 Subtype')
    ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('fig3_pca_scatter.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Saved: fig2_pca_scree.png, fig3_pca_scatter.png")

# ──────────────────────────────────────────────────────────────
# SECTION 4 – UMAP
# ──────────────────────────────────────────────────────────────
print("\n[4/8] Running UMAP...")
try:
    import umap
    reducer    = umap.UMAP(n_components=2, random_state=42, n_neighbors=10, min_dist=0.3)
    umap_coords = reducer.fit_transform(X_scaled)
    umap_df    = pd.DataFrame(umap_coords, index=X_rna.index,
                               columns=['UMAP1', 'UMAP2'])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for subtype, grp in umap_df.groupby(pam50):
        axes[0].scatter(grp['UMAP1'], grp['UMAP2'], label=subtype,
                        color=pam50_palette.get(subtype, '#aaaaaa'),
                        s=80, edgecolors='white', linewidth=0.5)
    axes[0].set(title='UMAP – PAM50 Subtype', xlabel='UMAP1', ylabel='UMAP2')
    axes[0].legend(fontsize=8)

    os_status  = meta.loc[umap_df.index, 'OS_STATUS'].fillna('Unknown')
    os_palette = {'LIVING': '#2ecc71', 'DECEASED': '#e74c3c', 'Unknown': '#aaaaaa'}
    for status, grp in umap_df.groupby(os_status):
        axes[1].scatter(grp['UMAP1'], grp['UMAP2'], label=status,
                        color=os_palette.get(status, '#aaaaaa'),
                        s=80, edgecolors='white', linewidth=0.5)
    axes[1].set(title='UMAP – OS Status', xlabel='UMAP1', ylabel='UMAP2')
    axes[1].legend()
    plt.tight_layout()
    plt.savefig('fig4_umap.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved: fig4_umap.png")
except ImportError:
    print("  umap-learn not installed – skipping UMAP (pip install umap-learn)")

# ──────────────────────────────────────────────────────────────
# SECTION 5 – Build Multimodal Feature Matrix
# ──────────────────────────────────────────────────────────────
print("\n[5/8] Building multimodal feature matrix...")

feat = pd.DataFrame(index=pca_df.index)

# Transcriptomic: PCs 1–10
feat = feat.join(pca_df[[f'PC{i}' for i in range(1, 11)]])

# Molecular features
def binary_encode(series, pos_label='Positive'):
    return (series == pos_label).astype(int)

feat['ER_pos']   = binary_encode(meta.loc[feat.index, 'ER_STATUS'])
feat['PR_pos']   = binary_encode(meta.loc[feat.index, 'PR_STATUS'])
feat['HER2_pos'] = binary_encode(meta.loc[feat.index, 'HER2_STATUS'])
feat['node_pos'] = binary_encode(meta.loc[feat.index, 'NODE_CODED'])
feat['mets_pos'] = binary_encode(meta.loc[feat.index, 'METASTASIS_CODED'])

mc = pd.to_numeric(meta.loc[feat.index, 'MUTATION_COUNT'], errors='coerce')
feat['log_mutation_count'] = np.log1p(mc.fillna(mc.median()))
feat['FGA'] = pd.to_numeric(meta.loc[feat.index, 'FRACTION_GENOME_ALTERED'], errors='coerce')
feat['FGA'] = feat['FGA'].fillna(feat['FGA'].median())

# Clinical features
feat['age'] = pd.to_numeric(meta.loc[feat.index, 'AGE'], errors='coerce')
feat['age'] = feat['age'].fillna(feat['age'].median())
stage_map = {
    'Stage I': 1, 'Stage IIA': 2, 'Stage IIB': 2,
    'Stage IIIA': 3, 'Stage IIIB': 3, 'Stage IIIC': 3,
    'Stage IV': 4, 'No_Conversion': np.nan
}
feat['tumour_stage_num'] = meta.loc[feat.index, 'CONVERTED_STAGE'].map(stage_map)
feat['tumour_stage_num'] = feat['tumour_stage_num'].fillna(feat['tumour_stage_num'].median())

print(f"  Feature matrix: {feat.shape[0]} patients x {feat.shape[1]} features")
print(f"    Transcriptomic (PCs 1-10): 10")
print(f"    Molecular (ER/PR/HER2/node/mets/mut/FGA): 7")
print(f"    Clinical (age, stage): 2")

# Correlation heatmap
fig, ax = plt.subplots(figsize=(13, 11))
corr = feat.corr(method='spearman')
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='coolwarm',
            center=0, vmin=-1, vmax=1, ax=ax,
            annot_kws={'size': 7}, linewidths=0.3)
ax.set_title("Spearman Correlation – Multimodal Feature Matrix", fontsize=12)
for boundary in [10, 17, 19]:
    ax.axhline(boundary, color='black', linewidth=1.5)
    ax.axvline(boundary, color='black', linewidth=1.5)
plt.tight_layout()
plt.savefig('fig5_correlation_heatmap.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Saved: fig5_correlation_heatmap.png")

# ──────────────────────────────────────────────────────────────
# SECTION 6 – PAM50 Classification
# ──────────────────────────────────────────────────────────────
print("\n[6/8] PAM50 classification with Random Forest (LOO-CV)...")

# NOTE ON SAMPLE COUNTS
# ---------------------
# All 32 samples are used in Sections 1–5 (preprocessing, PCA, feature matrix)
# and Section 8 (OS prediction), since those analyses require only RNA-seq and
# OS_STATUS data, both of which are available for all 32 patients.
#
# Here in Section 6, we can only classify samples that have a PAM50 ground-truth
# label. In the source file (tcga_brca_clinical_RNASeq.txt), 15 of the 32 samples
# have PAM50_SUBTYPE = NaN, with no usable proxy label in any other column
# (INTEGRATED_CLUSTERS_WITH_PAM50 and SIGCLUST_INTRINSIC_MRNA are also NaN for
# those 15). They are therefore excluded from classification only — not from the
# rest of the pipeline.
#
# LOO-CV EXPLAINED (n=17 labelled samples)
# -----------------------------------------
# A fixed train/test split (e.g. 80/20) would leave only ~3 test samples —
# too few for meaningful evaluation. Leave-One-Out cross-validation avoids this:
#
#   Fold 1:  train on samples 2–17 (n=16),  test on sample 1
#   Fold 2:  train on samples 1,3–17 (n=16), test on sample 2
#   ...
#   Fold 17: train on samples 1–16 (n=16),  test on sample 17
#
# Every sample is the test sample exactly once → 17 folds total.
# No data is held out permanently; the reported accuracy is the mean over all folds.

pam50_labels = meta.loc[feat.index, 'PAM50_SUBTYPE'].dropna()
pam50_labels = pam50_labels[pam50_labels != '']
feat_pam     = feat.loc[pam50_labels.index]
le           = LabelEncoder()
y_pam        = le.fit_transform(pam50_labels)

print(f"  Total samples in pipeline : 32")
print(f"  Samples with PAM50 label  : {len(y_pam)}  (used for classification only)")
print(f"  Samples missing PAM50     : {32 - len(y_pam)}  (used in all other sections)")
print(f"  LOO-CV folds              : {len(y_pam)}  (each fold: {len(y_pam)-1} train, 1 test)")
print(f"  Classes: {dict(zip(le.classes_, range(len(le.classes_))))}")

loo = LeaveOneOut()
rf  = RandomForestClassifier(n_estimators=200, random_state=42)

modalities = {
    'Transcriptomic (PC1-10)': [c for c in feat_pam.columns if c.startswith('PC')],
    'Molecular features':      ['ER_pos', 'PR_pos', 'HER2_pos', 'node_pos', 'mets_pos',
                                 'log_mutation_count', 'FGA'],
    'Clinical features':       ['age', 'tumour_stage_num'],
    'Multimodal (all)':        feat_pam.columns.tolist(),
}

from scipy.stats import chisquare
from sklearn.model_selection import cross_val_predict

n_classes        = len(le.classes_)
chance_level     = 1.0 / n_classes                          # 0.333 (3-class uniform)
majority_baseline = np.bincount(y_pam).max() / len(y_pam)  # always-predict-majority-class

print(f"  Chance level (uniform):         {chance_level:.3f}")
print(f"  Majority-class baseline:        {majority_baseline:.3f}  "
      f"(always predict '{le.classes_[np.bincount(y_pam).argmax()]}')")
print()
print(f"  {'Modality':<30}  {'Acc':>6}  {'Correct':>9}  {'chi2':>8}  {'p':>8}  {'Sig':>5}")
print("  " + "-" * 72)

results      = {}
chi2_results = {}
for name, cols in modalities.items():
    X_mod  = feat_pam[cols].fillna(0)
    # cross_val_predict gives the actual prediction for each LOO test sample
    y_pred = cross_val_predict(rf, X_mod, y_pam, cv=loo)
    correct = (y_pred == y_pam).sum()
    acc     = correct / len(y_pam)

    # Chi-square goodness-of-fit:
    #   H0: classifier performs at chance level (1/n_classes correct)
    #   Observed: [correct, wrong]
    #   Expected: [n * chance, n * (1 - chance)]
    exp_correct = len(y_pam) * chance_level
    exp_wrong   = len(y_pam) * (1 - chance_level)
    chi2_stat, p_val = chisquare(f_obs=[correct, len(y_pam) - correct],
                                  f_exp=[exp_correct, exp_wrong])
    sig = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else
          ("*"   if p_val < 0.05  else "ns"))

    results[name]      = acc
    chi2_results[name] = (chi2_stat, p_val, correct, sig)
    print(f"  {name:<30}  {acc:>6.3f}  {correct:>5}/{len(y_pam):<3}  "
          f"{chi2_stat:>8.3f}  {p_val:>8.4f}  {sig:>5}")

print()
print("  Chi-square H0: accuracy = chance level (1/3 = 0.333)")
print("  Significance: * p<0.05, ** p<0.01, *** p<0.001, ns = not significant")
print()
print(f"  NOTE: Majority-class baseline = {majority_baseline:.3f}. Any model at or below")
print(f"  this threshold offers no improvement over always predicting Luminal A.")

# Bar chart with chi-square p-values and majority baseline
fig, ax = plt.subplots(figsize=(10, 4))
names      = list(results.keys())
accs       = list(results.values())
colors_bar = ['#3498db', '#e67e22', '#2ecc71', '#e74c3c']
bars       = ax.barh(names, accs, color=colors_bar, edgecolor='white')

for bar, name, acc in zip(bars, names, accs):
    chi2_stat, p_val, correct, sig = chi2_results[name]
    label = (f'{acc:.3f}  ({correct}/{len(y_pam)})  '
             f'χ²={chi2_stat:.2f}, p={p_val:.4f} {sig}')
    ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
            label, va='center', fontsize=8.5)

# Reference lines
ax.axvline(chance_level,      color='grey',  linestyle='--', linewidth=1,
           label=f'Chance level (1/3 = {chance_level:.2f})')
ax.axvline(majority_baseline, color='black', linestyle=':',  linewidth=1.2,
           label=f'Majority-class baseline ({majority_baseline:.2f})')
ax.set(xlabel='LOO-CV Accuracy',
       title='PAM50 Subtype Classification by Modality\n'
             r'(χ² goodness-of-fit vs chance; dashed = chance, dotted = majority-class baseline)',
       xlim=(0, 1.35))
ax.legend(fontsize=8, loc='lower right')
plt.tight_layout()
plt.savefig('fig6_classification_accuracy.png', dpi=150, bbox_inches='tight')
plt.show()

# Feature importance
def get_modality(f):
    if f.startswith('PC'): return 'Transcriptomic'
    if f in ['ER_pos', 'PR_pos', 'HER2_pos', 'node_pos', 'mets_pos',
             'log_mutation_count', 'FGA']: return 'Molecular'
    return 'Clinical'

mod_colors = {'Transcriptomic': '#3498db', 'Molecular': '#e67e22', 'Clinical': '#2ecc71'}
X_all = feat_pam[modalities['Multimodal (all)']].fillna(0)
rf.fit(X_all, y_pam)
imp_df = pd.DataFrame({'Feature': X_all.columns,
                       'Importance': rf.feature_importances_}).sort_values('Importance')
imp_df['Modality'] = imp_df['Feature'].apply(get_modality)

fig, ax = plt.subplots(figsize=(9, 7))
ax.barh(imp_df['Feature'], imp_df['Importance'],
        color=[mod_colors[m] for m in imp_df['Modality']], edgecolor='white')
ax.set(xlabel='Gini Importance', title='Multimodal Random Forest – Feature Importance')
handles = [mpatches.Patch(color=c, label=m) for m, c in mod_colors.items()]
ax.legend(handles=handles, loc='lower right')
plt.tight_layout()
plt.savefig('fig7_feature_importance.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Saved: fig6_classification_accuracy.png, fig7_feature_importance.png")

# ──────────────────────────────────────────────────────────────
# SECTION 7 – Survival Analysis
# ──────────────────────────────────────────────────────────────
print("\n[7/8] Survival analysis (Kaplan-Meier)...")
try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test, logrank_test

    surv_idx = meta.index[meta['OS_MONTHS'].notna()]
    T = pd.to_numeric(meta.loc[surv_idx, 'OS_MONTHS'], errors='coerce')
    E = (meta.loc[surv_idx, 'OS_STATUS'] == 'DECEASED').astype(int)

    # KM by PAM50
    fig, ax = plt.subplots(figsize=(10, 6))
    colors_km = {'Luminal A': '#2ecc71', 'Luminal B': '#3498db', 'Basal-like': '#e74c3c'}
    subtypes_present = []
    for subtype, color in colors_km.items():
        mask = meta.loc[surv_idx, 'PAM50_SUBTYPE'] == subtype
        if mask.sum() < 3:
            continue
        subtypes_present.append(subtype)
        kmf = KaplanMeierFitter()
        kmf.fit(T[mask], E[mask], label=f'{subtype} (n={mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)

    pam_mask = meta.loc[surv_idx, 'PAM50_SUBTYPE'].isin(subtypes_present)
    lr = multivariate_logrank_test(
        T[pam_mask], meta.loc[surv_idx, 'PAM50_SUBTYPE'][pam_mask], E[pam_mask])
    ax.set(title=f'Overall Survival by PAM50 Subtype (log-rank p={lr.p_value:.3f})',
           xlabel='Time (months)', ylabel='Survival Probability')
    plt.tight_layout()
    plt.savefig('fig8_km_pam50.png', dpi=150, bbox_inches='tight')
    plt.show()

    # KM by PC1 split
    pc1       = pca_df.loc[surv_idx, 'PC1']
    pc1_group = (pc1 >= pc1.median()).map({True: 'PC1-High', False: 'PC1-Low'})

    fig, ax = plt.subplots(figsize=(9, 5))
    for grp, color in [('PC1-High', '#e74c3c'), ('PC1-Low', '#3498db')]:
        mask = pc1_group == grp
        kmf  = KaplanMeierFitter()
        kmf.fit(T[mask], E[mask], label=f'{grp} (n={mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)

    lr2 = logrank_test(T[pc1_group == 'PC1-High'], T[pc1_group == 'PC1-Low'],
                       E[pc1_group == 'PC1-High'], E[pc1_group == 'PC1-Low'])
    ax.set(title=f'Survival by PC1 split (log-rank p={lr2.p_value:.3f})',
           xlabel='Time (months)', ylabel='Survival Probability')
    plt.tight_layout()
    plt.savefig('fig9_km_pc1.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved: fig8_km_pam50.png, fig9_km_pc1.png")

    # PC1 gene loadings
    loadings    = pd.Series(pca.components_[0], index=top_genes, name='PC1_loading')
    top_loadings = pd.concat([loadings.nlargest(15), loadings.nsmallest(15)]).sort_values()
    fig, ax = plt.subplots(figsize=(9, 8))
    colors_load = ['#e74c3c' if v > 0 else '#3498db' for v in top_loadings]
    ax.barh(top_loadings.index, top_loadings.values, color=colors_load, edgecolor='white')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set(xlabel='PC1 Loading', title='Top 30 Genes Driving PC1')
    handles = [mpatches.Patch(color='#e74c3c', label='Positive loading'),
               mpatches.Patch(color='#3498db', label='Negative loading')]
    ax.legend(handles=handles)
    plt.tight_layout()
    plt.savefig('fig10_pc1_loadings.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved: fig10_pc1_loadings.png")

except ImportError:
    print("  lifelines not installed – skipping survival (pip install lifelines)")

# ──────────────────────────────────────────────────────────────
# SECTION 8 – Multimodal OS Prediction
# ──────────────────────────────────────────────────────────────
print("\n[8/8] Multimodal logistic regression for OS prediction...")

os_binary = (meta.loc[feat.index, 'OS_STATUS'] == 'DECEASED').astype(int)
X_os      = feat.fillna(0)
loo       = LeaveOneOut()
lr_model  = LogisticRegression(max_iter=1000, C=0.1, random_state=42)

y_pred_proba = np.zeros(len(os_binary))
for train_idx, test_idx in loo.split(X_os):
    lr_model.fit(X_os.iloc[train_idx], os_binary.iloc[train_idx])
    y_pred_proba[test_idx] = lr_model.predict_proba(X_os.iloc[test_idx])[:, 1]

auc = roc_auc_score(os_binary, y_pred_proba)
print(f"  LOO-CV AUC-ROC (OS prediction): {auc:.3f}")
print(f"  Deceased: {os_binary.sum()} / {len(os_binary)} patients")

fpr, tpr, _ = roc_curve(os_binary, y_pred_proba)
fig, axes   = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(fpr, tpr, color='#e74c3c', lw=2, label=f'Multimodal (AUC={auc:.3f})')
axes[0].plot([0, 1], [0, 1], '--', color='grey', linewidth=0.8, label='Random')
axes[0].set(xlabel='False Positive Rate', ylabel='True Positive Rate',
            title='ROC – OS Prediction (LOO-CV)')
axes[0].legend()

lr_model.fit(X_os, os_binary)
coef_df = pd.DataFrame({'Feature': X_os.columns,
                         'Log_Odds': lr_model.coef_[0]}).sort_values('Log_Odds')
coef_df['Modality'] = coef_df['Feature'].apply(get_modality)
axes[1].barh(coef_df['Feature'], coef_df['Log_Odds'],
             color=[mod_colors[m] for m in coef_df['Modality']], edgecolor='white')
axes[1].axvline(0, color='black', linewidth=0.8)
axes[1].set(xlabel='Log-Odds (positive = DECEASED)',
            title='Logistic Regression Coefficients')
handles = [mpatches.Patch(color=c, label=m) for m, c in mod_colors.items()]
axes[1].legend(handles=handles, loc='lower right')
plt.tight_layout()
plt.savefig('fig11_os_prediction.png', dpi=150, bbox_inches='tight')
plt.show()
print("  Saved: fig11_os_prediction.png")

print("\n" + "=" * 60)
print("Analysis complete.")
print("Figures saved as fig1_normalisation.png ... fig11_os_prediction.png")
print("=" * 60)
