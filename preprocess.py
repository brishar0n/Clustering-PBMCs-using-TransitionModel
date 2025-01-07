import scanpy as sc

def load_data(data_path):
    """
    Loads the PBMC dataset from 10X Genomics.
    """
    adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    return adata
