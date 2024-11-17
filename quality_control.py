import scanpy as sc

def perform_quality_control(adata):
    """
    Filters cells with high mitochondrial gene percentage and low gene counts.
    """
    # Calculate QC metrics
    adata.obs['pct_counts_mt'] = (adata[:, adata.var_names.str.startswith('MT-')].X.sum(1) / adata.X.sum(1)) * 100
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['pct_counts_mt'] < 5, :]
    
    return adata