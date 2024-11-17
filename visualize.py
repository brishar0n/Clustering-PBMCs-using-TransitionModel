import scanpy as sc

def plot_data(adata):
    """
    Generates UMAP plots for visualization.
    """
    # Compute UMAP embedding
    sc.tl.umap(adata)
    
    # Plot UMAP
    sc.pl.umap(adata, color=['leiden', 'pct_counts_mt', 'n_genes'])
