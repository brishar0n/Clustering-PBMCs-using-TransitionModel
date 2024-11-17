# Clustering 3K PBMCs with Scanpy

This project demonstrates the clustering of 3K Peripheral Blood Mononuclear Cells (PBMCs) using Scanpy, a Python package for analyzing single-cell gene expression data, and without the usage of machine learning algorithms.

## Requirements

Make sure you have the following Python packages installed:
- `scanpy`
- `anndata`

You can install these with:
```bash
pip install scanpy anndata
```

## Data Preparation

The dataset used in this project is located in the `data/filtered_gene_bc_matrices/hg19/` directory. Ensure the following files are present:
- `matrix.mtx.gz` (compressed gene count matrix, needs decompressing)
- `genes.tsv` (gene names and IDs)
- `barcodes.tsv` (cell barcodes)

**Note:** The `matrix.mtx.gz` file in this repository is compressed. You need to decompress it for the code to run successfully.

### Decompress the Matrix File

Run the following command to decompress `matrix.mtx.gz`:

```bash
gunzip data/filtered_gene_bc_matrices/hg19/matrix.mtx.gz
```

## Running the Project

To run the project, follow these steps:

1. **Navigate to the Project Directory**  
   Open a terminal and change to the project directory:

   ```bash
   cd path/to/Clustering-3kPBMCs-using-Scanpy

2. **Run the Main Script**
   Execute the main Python script to start the analysis:

   ```bash
   python3 main.py
   ```

## Project Outline
- data/filtered_gene_bc_matrices/hg19/matrix.mtx.gz: Gene expression matrix (compressed version)
- data/filtered_gene_bc_matrices/hg19/genes.tsv: Gene metadata
- data/filtered_gene_bc_matrices/hg19/barcodes.tsv: Barcode metadata
- main.py: Main script to load data and perform clustering
- preprocess.py: Contains functions for data loading and preprocessing
- quality_control.py: Conducts quality control operations and normalizes the data
- clustering.py: Implements clustering methods to classify cells based on expression profiles.
- visualize.py: Generates visualizations to help interpret clustering results.

### Reference
[Clustering 3K PBMCs with Scanpy Tutorial by Galaxy Training](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html)
