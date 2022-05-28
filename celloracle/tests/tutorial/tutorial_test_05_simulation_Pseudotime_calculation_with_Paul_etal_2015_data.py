import copy
import glob
import time
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.auto import tqdm

import celloracle as co
from celloracle.applications import Pseudotime_calculator
print("co version: ", co.__version__)

plt.rcParams['figure.figsize'] = [5, 5]
plt.rcParams["savefig.dpi"] = 300

def main():
    test_tutorial_pseudotime()
    
def test_tutorial_pseudotime():
    # Start analysis
    import time
    start = time.time()

    # 1.0 Load data
    print("\n1.0 Load data")
    # Option1. Load oracle data and make pt from oracle object
    oracle = co.data.load_tutorial_oracle_object()
    pt = Pseudotime_calculator(oracle_object=oracle)
    # Option2. Load anndata and make pt from anndata
    adata = co.data.load_Paul2015_data()
    pt = Pseudotime_calculator(adata=adata,
                               obsm_key="X_draw_graph_fa",
                               cluster_column_name="louvain_annot")

    # 2. Pseudotime calculation
    print("\n2. Pseudotime calculation")
    print("Clustering name: ", pt.cluster_column_name)
    print("Cluster list", pt.cluster_list)
    # Check data
    pt.plot_cluster(fontsize=8)
    # Here, clusters can be classified into either MEP lineage or GMP lineage
    clusters_in_ME_lineage = ['Ery_0', 'Ery_1', 'Ery_2', 'Ery_3', 'Ery_4', 'Ery_5',
                              'Ery_6', 'Ery_7', 'Ery_8', 'Ery_9', 'MEP_0', 'Mk_0']
    clusters_in_GM_lineage = ['GMP_0', 'GMP_1', 'GMP_2', 'GMPl_0', 'GMPl_1', 'Gran_0',
                              'Gran_1', 'Gran_2', 'Gran_3', 'Mo_0', 'Mo_1', 'Mo_2']
    # Make a dictionary
    lineage_dictionary = {"Lineage_ME": clusters_in_ME_lineage,
               "Lineage_GM": clusters_in_GM_lineage}
    # Input lineage information into pseudotime object
    pt.set_lineage(lineage_dictionary=lineage_dictionary)
    # Visualize lineage information
    pt.plot_lineages()
    # Set root cell
    # Estimated root cell name for each lineage
    root_cells = {"Lineage_ME": "1539", "Lineage_GM": "2244"}
    pt.set_root_cells(root_cells=root_cells)
    # Check root cell and lineage
    pt.plot_root_cells()
    # Check diffusion map data.
    "X_diffmap" in pt.adata.obsm
    # Calculate pseudotime
    pt.get_pseudotime_per_each_lineage()
    # Check results
    pt.plot_pseudotime(cmap="rainbow")
    # Check result
    pt.adata.obs[["Pseudotime"]].head()

    # 3. Save data
    print("\n3. Save data")

    # If you started calculation with oracle object,
    #  Add calculated pseudotime data to the oracle object
    oracle.adata.obs = pt.adata.obs
    oracle.to_hdf5("oracle_with_pt.celloracle.oracle")
    # If you started calculation with anndata object,
    #  Add calculated pseudotime data to the oracle object
    adata.obs = pt.adata.obs
    adata.write_h5ad("anndata_with_pt.h5ad")

    # Remove files
    print("Remove files")
    for i in ["oracle_with_pt.celloracle.oracle", "anndata_with_pt.h5ad"]:
        os.remove(i)

    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
