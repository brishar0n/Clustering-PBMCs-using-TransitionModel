import scanpy as sc

def perform_clustering(adata):
    """
    Performs PCA, computes neighborhood graph, and clusters cells.
    """
    # Normalize and log-transform the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes and scale the data
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    
    # Perform PCA
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Compute the neighborhood graph and cluster cells
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata)
    
    return adata
