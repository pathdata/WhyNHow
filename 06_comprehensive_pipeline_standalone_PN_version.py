"""
Comprehensive Single cell RNA-seq Pipeline - Standalone Version
Handles GEX/TCR mapping correctly and doesn't require external modules
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
from pathlib import Path
import warnings
from typing import Optional, Tuple, Dict, List
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')


def load_metadata_files_standalone(supplement_dir: str) -> Tuple[set, pd.DataFrame, pd.DataFrame]:
    """
    Load metadata files from the supplement directory
    
    Parameters:
    -----------
    supplement_dir : str
        Path to supplement directory containing metadata files
    
    Returns:
    --------
    exclude_genes : set
        Set of genes to exclude
    cart_index : DataFrame
        Sample metadata
    barcode_metadata : DataFrame
        Cell-level metadata
    """
    print(f"\nLoading metadata files from: {supplement_dir}")
    
    # Load exclude genes
    exclude_genes_path = os.path.join(supplement_dir, 'excludeGenes.tsv')
    if os.path.exists(exclude_genes_path):
        exclude_df = pd.read_csv(exclude_genes_path, sep='\t', header=None)
        exclude_genes = set(exclude_df[0].tolist())
        print(f"  Loaded {len(exclude_genes)} genes to exclude")
    else:
        print(f"  [WARNING] excludeGenes.tsv not found at {exclude_genes_path}")
        exclude_genes = set()
    
    # Load sample metadata
    cart_index_path = os.path.join(supplement_dir, 'CARPALL_scRNAseq_CART_metadata.txt')
    if os.path.exists(cart_index_path):
        cart_index = pd.read_csv(cart_index_path, sep='\t')
        print(f"  Loaded sample metadata for {len(cart_index)} samples")
    else:
        print(f"  [WARNING] CARPALL_scRNAseq_CART_metadata.txt not found at {cart_index_path}")
        cart_index = pd.DataFrame()
    
    # Load cell-level metadata
    barcode_path = os.path.join(supplement_dir, 'CARTcell_barcodes.txt')
    if os.path.exists(barcode_path):
        barcode_metadata = pd.read_csv(barcode_path, sep='\t')
        print(f"  Loaded cell metadata for {len(barcode_metadata)} cells")
    else:
        print(f"  [WARNING] CARTcell_barcodes.txt not found at {barcode_path}")
        barcode_metadata = pd.DataFrame()
    
    return exclude_genes, cart_index, barcode_metadata


def load_vdj_data(vdj_path: str, sample_id: str) -> Optional[pd.DataFrame]:
    """Load and process VDJ (TCR/BCR) contig annotations"""
    if not os.path.exists(vdj_path):
        return None
    
    try:
        print(f"   Loading VDJ data from: {vdj_path}")
        vdj_df = pd.read_csv(vdj_path)
        
        # Filter for productive, high-confidence contigs
        vdj_df = vdj_df[
            (vdj_df['productive'] == True) & 
            (vdj_df['high_confidence'] == True) &
            (vdj_df['is_cell'] == True)
        ].copy()
        
        # Add sample prefix to barcodes
        vdj_df['barcode_full'] = vdj_df['barcode'].apply(
            lambda x: f"{sample_id}_{x}" if not x.startswith(sample_id) else x
        )
        
        print(f"   [OK] Loaded {len(vdj_df)} productive contigs")
        return vdj_df
        
    except Exception as e:
        print(f"   [WARNING] Error loading VDJ data: {str(e)}")
        return None


def process_vdj_per_cell(vdj_df: pd.DataFrame) -> pd.DataFrame:
    """Process VDJ data to get one row per cell with TCR/BCR information"""
    if vdj_df is None or len(vdj_df) == 0:
        return pd.DataFrame()
    
    cell_data = []
    
    for barcode in vdj_df['barcode_full'].unique():
        cell_contigs = vdj_df[vdj_df['barcode_full'] == barcode]
        
        # Get TCR alpha chain
        tra_contigs = cell_contigs[cell_contigs['chain'] == 'TRA']
        tra_info = None
        if len(tra_contigs) > 0:
            tra_best = tra_contigs.loc[tra_contigs['umis'].idxmax()]
            tra_info = {
                'TRA_v_gene': str(tra_best['v_gene']) if pd.notna(tra_best['v_gene']) else '',
                'TRA_j_gene': str(tra_best['j_gene']) if pd.notna(tra_best['j_gene']) else '',
                'TRA_c_gene': str(tra_best['c_gene']) if pd.notna(tra_best['c_gene']) else '',
                'TRA_cdr3': str(tra_best['cdr3']) if pd.notna(tra_best['cdr3']) else '',
                'TRA_productive': True
            }
        
        # Get TCR beta chain
        trb_contigs = cell_contigs[cell_contigs['chain'] == 'TRB']
        trb_info = None
        if len(trb_contigs) > 0:
            trb_best = trb_contigs.loc[trb_contigs['umis'].idxmax()]
            trb_info = {
                'TRB_v_gene': str(trb_best['v_gene']) if pd.notna(trb_best['v_gene']) else '',
                'TRB_d_gene': str(trb_best['d_gene']) if pd.notna(trb_best['d_gene']) else '',
                'TRB_j_gene': str(trb_best['j_gene']) if pd.notna(trb_best['j_gene']) else '',
                'TRB_c_gene': str(trb_best['c_gene']) if pd.notna(trb_best['c_gene']) else '',
                'TRB_cdr3': str(trb_best['cdr3']) if pd.notna(trb_best['cdr3']) else '',
                'TRB_productive': True
            }
        
        # Get clonotype
        clonotype = str(cell_contigs['raw_clonotype_id'].iloc[0]) if len(cell_contigs) > 0 and pd.notna(cell_contigs['raw_clonotype_id'].iloc[0]) else ''
        
        # Combine information
        cell_info = {'CellID': barcode}
        if tra_info:
            cell_info.update(tra_info)
        if trb_info:
            cell_info.update(trb_info)
        if clonotype:
            cell_info['clonotype_id'] = clonotype
        
        cell_info['has_TCR'] = (tra_info is not None) and (trb_info is not None)
        
        cell_data.append(cell_info)
    
    vdj_cells = pd.DataFrame(cell_data)
    
    # Ensure consistent types for all columns
    for col in vdj_cells.columns:
        if col in ['TRA_productive', 'TRB_productive', 'has_TCR']:
            vdj_cells[col] = vdj_cells[col].fillna(False).astype(bool)
        elif col != 'CellID':
            vdj_cells[col] = vdj_cells[col].fillna('').astype(str)
    
    return vdj_cells


def load_protein_data(h5_path: str) -> Tuple[Optional[np.ndarray], Optional[List[str]]]:
    """Load protein/ADT data from CITE-seq h5 file"""
    try:
        adata_full = sc.read_10x_h5(h5_path, gex_only=False)
        
        if 'feature_types' in adata_full.var.columns:
            protein_mask = adata_full.var['feature_types'] == 'Antibody Capture'
            n_proteins = protein_mask.sum()
            
            if n_proteins > 0:
                protein_data = adata_full[:, protein_mask].X
                protein_names = adata_full.var_names[protein_mask].tolist()
                
                if hasattr(protein_data, 'toarray'):
                    protein_data = protein_data.toarray()
                
                return protein_data, protein_names
        
        return None, None
        
    except Exception as e:
        return None, None


def clr_normalize_protein(protein_data: np.ndarray) -> np.ndarray:
    """CLR (Centered Log Ratio) normalization for protein data"""
    protein_data = protein_data + 1
    geom_mean = np.exp(np.mean(np.log(protein_data), axis=1, keepdims=True))
    protein_clr = np.log(protein_data / geom_mean)
    return protein_clr


def process_sample_comprehensive(
    gex_path: str,
    sample_id: str,
    supplement_dir: str,
    vdj_path: Optional[str] = None,
    nFeature_RNA_min: int = 300,
    percentmt: float = 10,
    nCountRNA_min: int = 1000,
    npcs1: int = 75,
    clustres: float = 1.0,
    finalMinDist: float = 0.5,
    finalNN: int = 50,
    subCAR: bool = True,
    save_h5ad: bool = True,
    output_dir: Optional[str] = None
) -> sc.AnnData:
    """
    Process a single sample through comprehensive CITE-seq pipeline
    
    Parameters:
    -----------
    gex_path : str
        Path to GEX h5 file
    sample_id : str
        Sample identifier (GEX sample ID)
    supplement_dir : str
        Path to directory containing metadata files
    vdj_path : str, optional
        Path to VDJ CSV file
    output_dir : str, optional
        Directory to save output (defaults to supplement_dir)
    """
    
    if output_dir is None:
        output_dir = supplement_dir
    
    # Load metadata files
    exclude_genes, cart_index, barcode_metadata = load_metadata_files_standalone(supplement_dir)
    
    print(f"\n{'='*80}")
    print(f"COMPREHENSIVE CITE-SEQ PIPELINE: {sample_id}")
    print(f"{'='*80}\n")
    
    # Load RNA data
    print("Step 1: Loading RNA data...")
    adata = sc.read_10x_h5(gex_path, gex_only=True)
    adata.var_names_make_unique()
    adata.obs['orig.ident'] = sample_id
    print(f"   Initial: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Load protein data
    print("\nStep 2: Checking for protein/ADT data...")
    protein_data, protein_names = load_protein_data(gex_path)
    has_protein = protein_data is not None
    if has_protein:
        print(f"   [OK] Found {len(protein_names)} proteins")
    else:
        print(f"   No protein data found")
    
    # Load VDJ data
    print("\nStep 3: Loading VDJ data...")
    vdj_df = None
    vdj_cells = None
    if vdj_path and os.path.exists(vdj_path):
        vdj_df = load_vdj_data(vdj_path, sample_id)
        if vdj_df is not None and len(vdj_df) > 0:
            vdj_cells = process_vdj_per_cell(vdj_df)
            print(f"   [OK] Processed VDJ for {len(vdj_cells)} cells")
    else:
        print(f"   No VDJ data available")
    
    # QC metrics
    print("\nStep 4: Calculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    if 'total_counts' not in adata.obs.columns:
        adata.obs['total_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)
    if 'n_genes' not in adata.obs.columns:
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    
    adata.obs['percent.mt'] = adata.obs['pct_counts_mt']
    adata.obs['CellID'] = [f"{sample_id}_{barcode}" for barcode in adata.obs.index]
    
    # QC filtering
    print("\nStep 5: Applying QC filters...")
    initial_cells = adata.n_obs
    adata = adata[adata.obs['total_counts'] > nCountRNA_min, :].copy()
    adata = adata[adata.obs['n_genes'] > nFeature_RNA_min, :].copy()
    adata = adata[adata.obs['percent.mt'] <= percentmt, :].copy()
    print(f"   Retained: {adata.n_obs:,} / {initial_cells:,} cells")
    
    # CAR T-cell subsetting
    if subCAR and len(barcode_metadata) > 0:
        print("\nStep 6: Subsetting for CAR T-cells...")
        car_cellids = set(barcode_metadata['CellID'].values)
        adata_car_mask = adata.obs['CellID'].isin(car_cellids)
        n_car = adata_car_mask.sum()
        if n_car > 0:
            adata = adata[adata_car_mask, :].copy()
            print(f"   Retained: {adata.n_obs:,} CAR T-cells")
    
    # Gene filtering
    print("\nStep 7: Removing unwanted genes...")
    if len(exclude_genes) > 0:
        keep_features = [gene for gene in adata.var_names if gene not in exclude_genes]
        adata = adata[:, keep_features].copy()
    print(f"   Retained: {adata.n_vars:,} genes")
    
    # RNA normalization
    print("\nStep 8: Normalizing RNA data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', subset=True)
    sc.pp.scale(adata, max_value=10)
    
    # Protein processing
    if has_protein and protein_data.shape[0] == len(adata.obs):
        print("\nStep 9: Processing protein data...")
        adata.obsm['protein'] = protein_data[adata.obs.index.isin(adata.obs.index), :]
        adata.uns['protein_names'] = protein_names
        protein_clr = clr_normalize_protein(adata.obsm['protein'])
        adata.obsm['protein_clr'] = protein_clr
    
    # PCA
    print(f"\nStep 10: Running PCA ({npcs1} components)...")
    sc.tl.pca(adata, n_comps=npcs1)
    
    # Neighbors and clustering
    print(f"\nStep 11: Building neighbor graph...")
    sc.pp.neighbors(adata, n_neighbors=finalNN, n_pcs=npcs1)
    
    print(f"\nStep 12: Clustering (resolution {clustres})...")
    sc.tl.leiden(adata, resolution=clustres, key_added='seurat_clusters')
    
    # UMAP
    print(f"\nStep 13: Running UMAP...")
    sc.tl.umap(adata, min_dist=finalMinDist)
    
    # Cell cycle scoring
    print("\nStep 14: Cell cycle scoring...")
    s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6']
    g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2']
    
    s_genes_present = [g for g in s_genes if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]
    
    if len(s_genes_present) > 0 and len(g2m_genes_present) > 0:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
    
    # Add VDJ annotations
    if vdj_cells is not None and len(vdj_cells) > 0:
        print("\nStep 15: Adding VDJ annotations...")
        vdj_dict = vdj_cells.set_index('CellID').to_dict('index')
        for col in vdj_cells.columns:
            if col != 'CellID':
                adata.obs[col] = adata.obs['CellID'].map(
                    lambda x: vdj_dict.get(x, {}).get(col, np.nan)
                )
                
                # Convert boolean columns to categorical to avoid h5ad save issues
                if col in ['TRA_productive', 'TRB_productive', 'has_TCR']:
                    adata.obs[col] = adata.obs[col].fillna(False).astype(bool)
                # Convert string columns with NaN to empty strings
                elif adata.obs[col].dtype == object:
                    adata.obs[col] = adata.obs[col].fillna('').astype(str)
    
    # Add sample metadata
    print("\nStep 16: Adding sample metadata...")
    if len(cart_index) > 0:
        if 'GEX_Sample' in cart_index.columns:
            cart_index = cart_index.rename(columns={'GEX_Sample': 'orig.ident'})
        sample_meta = cart_index[cart_index['orig.ident'] == sample_id]
        if len(sample_meta) > 0:
            for col in ['Patient', 'TCR_Sample', 'Sort', 'Source', 'Timepoint', 'Timepoint_bin']:
                if col in sample_meta.columns:
                    value = sample_meta[col].values[0]
                    # Convert to string if not already
                    adata.obs[col] = str(value) if pd.notna(value) else ''
    
    # Add cell type annotations
    if len(barcode_metadata) > 0:
        barcode_dict = barcode_metadata.set_index('CellID').to_dict('index')
        for col in ['CellType1', 'CellType2', 'CARTcell']:
            if col in barcode_metadata.columns:
                adata.obs[col] = adata.obs['CellID'].map(
                    lambda x: barcode_dict.get(x, {}).get(col, '')
                )
                # Ensure string type
                adata.obs[col] = adata.obs[col].fillna('').astype(str)
    
    # Save
    if save_h5ad:
        # Clean up data types before saving to avoid h5ad compatibility issues
        print("\nStep 17: Preparing data for saving...")
        
        # Convert all object dtype columns to strings (except for special ones)
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object:
                # Check if it's actually categorical or should be string
                try:
                    # Try to convert to string, handling any remaining NaN
                    adata.obs[col] = adata.obs[col].fillna('').astype(str)
                except:
                    pass
        
        # Ensure boolean columns are proper booleans
        for col in ['TRA_productive', 'TRB_productive', 'has_TCR']:
            if col in adata.obs.columns:
                adata.obs[col] = adata.obs[col].fillna(False).astype(bool)
        
        output_filename = f"{sample_id}_comprehensive.{nFeature_RNA_min}.{percentmt}.{nCountRNA_min}.{npcs1}.h5ad"
        output_path = os.path.join(output_dir, output_filename)
        print(f"Saving to {output_path}")
        adata.write(output_path)
        print(f"[OK] Saved successfully")
    
    print(f"\n{'='*80}")
    print(f"[OK] COMPLETED: {sample_id}")
    print(f"   Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")
    print(f"{'='*80}\n")
    
    return adata


def process_multiple_samples(
    sample_ids: List[str],
    data_dir: str,
    supplement_dir: str,
    metadata_path: str,
    output_dir: Optional[str] = None,
    **kwargs
) -> Dict[str, sc.AnnData]:
    """
    Process multiple samples with correct GEX/TCR mapping
    
    Parameters:
    -----------
    sample_ids : list
        List of GEX sample IDs
    data_dir : str
        Path to data directory
    supplement_dir : str
        Path to supplement directory with metadata files
    metadata_path : str
        Path to CARPALL metadata file
    output_dir : str, optional
        Directory to save outputs
    """
    
    if output_dir is None:
        output_dir = supplement_dir
    
    # Load metadata and create GEX to TCR mapping
    print("="*80)
    print("LOADING METADATA")
    print("="*80)
    metadata = pd.read_csv(metadata_path, sep='\t')
    gex_to_tcr = dict(zip(metadata['GEX_Sample'], metadata['TCR_Sample']))
    print(f"Loaded metadata for {len(metadata)} samples")
    print(f"Created {len(gex_to_tcr)} GEX->TCR mappings\n")
    
    adata_dict = {}
    
    for gex_sample_id in sample_ids:
        print(f"\n{'='*80}")
        print(f"PROCESSING: {gex_sample_id}")
        print(f"{'='*80}")
        
        # Get GEX file path
        gex_dir = os.path.join(data_dir, 'GEX', gex_sample_id)
        h5_file = os.path.join(gex_dir, 'filtered_feature_bc_matrix.h5')
        
        if not os.path.exists(h5_file):
            print(f"[WARNING] GEX file not found: {h5_file}\n")
            continue
        
        print(f"GEX file: {h5_file}")
        
        # Get corresponding TCR sample ID and VDJ path
        tcr_sample_id = gex_to_tcr.get(gex_sample_id)
        vdj_path = None
        
        if tcr_sample_id:
            print(f"TCR sample ID: {tcr_sample_id}")
            vdj_dir = os.path.join(data_dir, 'VDJ', tcr_sample_id)
            vdj_file = os.path.join(vdj_dir, 'filtered_contig_annotations.csv')
            
            if os.path.exists(vdj_file):
                vdj_path = vdj_file
                print(f"VDJ file: {vdj_file}")
            else:
                print(f"[WARNING] VDJ file not found: {vdj_file}")
        else:
            print(f"[WARNING] No TCR sample ID found for {gex_sample_id}")
        
        # Process sample
        try:
            adata = process_sample_comprehensive(
                gex_path=h5_file,
                sample_id=gex_sample_id,
                supplement_dir=supplement_dir,
                vdj_path=vdj_path,
                output_dir=output_dir,
                **kwargs
            )
            adata_dict[gex_sample_id] = adata
            
        except Exception as e:
            print(f"\n[ERROR] Failed to process {gex_sample_id}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    return adata_dict


if __name__ == "__main__":
    
    # MODIFY THESE PATHS FOR YOUR ENVIRONMENT
    base_dir = Path(r"D:\CAT_CART_paper_2023-v.1.0.0")
    data_dir = base_dir / 'data'
    supplement_dir = base_dir / 'code' / 'supplement'
    metadata_path = supplement_dir / 'CARPALL_scRNAseq_CART_metadata.txt'
    output_dir = supplement_dir / 'outputs'
    
    # Create output directory
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print("="*80)
    print("COMPREHENSIVE CITE-SEQ PIPELINE - STANDALONE")
    print("="*80)
    print(f"\nData directory: {data_dir}")
    print(f"Supplement directory: {supplement_dir}")
    print(f"Metadata file: {metadata_path}")
    print(f"Output directory: {output_dir}\n")

    #To process specific sample use sample_ids populated in the list
    
    # Sample IDs (GEX sample IDs)
    sample_ids = [
        'CAR-T10191179',
        'CAR-T10191180',
        'CAR-T10191181',
        'CAR-T9422250',
        'CAR-T9422251'
    ]
    # To process all sample use below two lines
    # gex_dir = data_dir / 'GEX'
    # sample_ids = [f.name for f in gex_dir.glob("CAR-T*")]
    
    # Process samples
    adata_dict = process_multiple_samples(
        sample_ids=sample_ids,
        data_dir=str(data_dir),
        supplement_dir=str(supplement_dir),
        metadata_path=str(metadata_path),
        output_dir=str(output_dir),
        nFeature_RNA_min=300,
        percentmt=10,
        nCountRNA_min=1000,
        npcs1=75,
        clustres=1.0,
        finalMinDist=0.5,
        finalNN=50,
        subCAR=True,
        save_h5ad=True
    )
    
    print(f"\n{'='*80}")
    print(f"[OK] Processed {len(adata_dict)} / {len(sample_ids)} samples")
    print(f"{'='*80}\n")
    
    for sample_id, adata in adata_dict.items():
        print(f"  {sample_id}: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
