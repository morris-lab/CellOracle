import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [6, 4.5]


def main():
    test_tutorial_scRNAseq_data_processing()

def test_tutorial_scRNAseq_data_processing():
    # Start analysis
    import time
    start = time.time()



    # 1. Load data
    print("\n1. Load data")
    # Download dataset. You can change the code blow to use your data.
    adata = sc.datasets.paul15()

    # 2. Filtering
    print("\n2. Filtering")
    # Only consider genes with more than 1 count
    sc.pp.filter_genes(adata, min_counts=1)


    # 3. Normalizatioin
    print("\n3. Normalization")
    # Normalize gene expression matrix with total UMI count per cell
    sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

    # 4. Identificatin of highly variable genes
    print("\n4. Identificatin of highly variable genes")
    # Select top 2000 highly-variable genes
    filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                                  flavor='cell_ranger',
                                                  n_top_genes=2000,
                                                  log=False)
    # Subset the genes
    adata = adata[:, filter_result.gene_subset]
    # Renormalize after filtering
    sc.pp.normalize_per_cell(adata)

    # 5. Log transformation
    print("\n5. Log transformation")
    # keep raw cont data before log transformation
    adata.raw = adata
    adata.layers["raw_count"] = adata.raw.X.copy()
    # Log transformation and scaling
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    # 6. PCA and neighbor calculations
    print("\n6. PCA and neighbor calculations")
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    # Diffusion map
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
    sc.tl.diffmap(adata)
    # Calculate neihbors again based on diffusionmap
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

    # 7. Cell clustering
    print("\n7. Cell clustering")
    sc.tl.louvain(adata, resolution=0.8)

    # 8. Dimensionality reduction using PAGA and force-directed graphs
    print("\n8. Dimensionality reduction using PAGA and force-directed graphs")
    # PAGA graph construction
    sc.tl.paga(adata, groups='louvain')
    plt.rcParams["figure.figsize"] = [6, 4.5]
    sc.pl.paga(adata)
    sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
    sc.pl.draw_graph(adata, color='louvain', legend_loc='on data')

    # 9. Checkdata
    print("\n9. Check data")
    plt.rcParams["figure.figsize"] = [4.5, 4.5]
    markers = {"Erythroids":["Gata1", "Klf1", "Gypa", "Hba-a2"],
               "Megakaryocytes":["Itga2b", "Pbx1", "Sdpr", "Vwf"],
                "Granulocytes":["Elane", "Cebpe", "Ctsg", "Mpo", "Gfi1"],
                "Monocytes":["Irf8", "Csf1r", "Ctsg", "Mpo"],
                "Mast_cells":["Cma1", "Gzmb", "Kit"],
                "Basophils":["Mcpt8", "Prss34"]
                }
    for cell_type, genes in markers.items():
        print(f"marker gene of {cell_type}")
        sc.pl.draw_graph(adata, color=genes, use_raw=False, ncols=2)
        plt.show()

    adata.write_h5ad("Paul_etal_15.h5ad")

    # 6. Remove files
    print("Remove files")
    for i in ["Paul_etal_15.h5ad"]:
        os.remove(i)
    os.remove("data/paul15/paul15.h5")
    os.removedirs("data/paul15")

    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
