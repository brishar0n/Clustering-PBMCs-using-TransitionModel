from preprocess import load_data
from quality_control import perform_quality_control
from clustering import perform_clustering
from visualize import plot_data

def main():
    print("Loading data...")
    adata = load_data('data/filtered_gene_bc_matrices/hg19/')
    
    print("Performing quality control...")
    adata = perform_quality_control(adata)
    
    print("Performing clustering...")
    adata = perform_clustering(adata)
    
    print("Visualizing results...")
    plot_data(adata)

if __name__ == "__main__":
    main()