import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


import celloracle as co
print("co version: ", co.__version__)

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

def main():
    test_tutorial_network_analysis()

def test_tutorial_network_analysis():
    use_reduced_data = True

    # Start analysis
    import time
    start = time.time()

    save_folder = "figures"
    os.makedirs(save_folder, exist_ok=True)


    # 1. Load data
    print("\n1. Load data")
    adata = co.data.load_Paul2015_data()

    # Testing with whole data will take long time. We can test function with subsampled data.
    if use_reduced_data:
        cells = adata.obs.index[adata.obs.louvain_annot.isin(["MEP_0", "GMPl_0"])]
        adata = adata[cells, :]

    # (Option) Down sampling
    print(f"Cell number is :{adata.shape[0]}")
    print(f"Gene number is :{adata.shape[1]}")
    # Random downsampling into 30K cells if the anndata object include more than 30 K cells.
    n_cells_downsample = 30000
    if adata.shape[0] > n_cells_downsample:
        # Let's dowmsample into 30K cells
        sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)
    print(f"Cell number is :{adata.shape[0]}")
    # Load TF info which was made from mouse cell atlas dataset.
    base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()
    # Check data
    base_GRN.head()

    # 2. Make Oracle object
    print("\n2. Make Oracle object")
    # Instantiate Oracle object
    oracle = co.Oracle()
    # Check data in anndata
    print("Metadata columns :", list(adata.obs.columns))
    print("Dimensional reduction: ", list(adata.obsm.keys()))
    # In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
    adata.X = adata.layers["raw_count"].copy()
    # Instantiate Oracle object.
    oracle.import_anndata_as_raw_count(adata=adata,
                                       cluster_column_name="louvain_annot",
                                       embedding_name="X_draw_graph_fa")
    # You can load TF info dataframe with the following code.
    oracle.import_TF_data(TF_info_matrix=base_GRN)

    ## 2.3. (Optional) Add TF-target gene pair manually
    print("\n (Optional) Add TF-target gene pair manually")
    # Download TF-target gene pair data
    from celloracle.utility.data_download_from_web import download_demo_data
    for i in ["TF_data_in_Paul15.csv"]:
        download_demo_data(file=i)
    # Load the TF and target gene information from Paul et al. (2015).
    Paul_15_data = pd.read_csv("TF_data_in_Paul15.csv")
    Paul_15_data
    # Make dictionary: dictionary key is TF and dictionary value is list of target genes.
    TF_to_TG_dictionary = {}
    for TF, TGs in zip(Paul_15_data.TF, Paul_15_data.Target_genes):
        # convert target gene to list
        TG_list = TGs.replace(" ", "").split(",")
        # store target gene list in a dictionary
        TF_to_TG_dictionary[TF] = TG_list
    # We invert the dictionary above using a utility function in celloracle.
    TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
    # Add TF information
    oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

    # 3. KNN imputation
    print("\n3. KNN imputation")
    # Perform PCA
    oracle.perform_PCA()
    # Select important PCs
    plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    plt.axvline(n_comps, c="k")
    plt.show()
    print(n_comps)
    n_comps = min(n_comps, 50)
    # KNN imputation
    n_cell = oracle.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = int(0.025*n_cell)
    print(f"Auto-selected k is :{k}")
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                          b_maxl=k*4, n_jobs=4)

    # 4. Save and Load.
    print("\n4. Save and Load.")
    # Save oracle object.
    oracle.to_hdf5("Paul_15_data.celloracle.oracle")
    # Load file.
    oracle = co.load_hdf5("Paul_15_data.celloracle.oracle")

    # 5. GRN calculation
    print("\n5. GRN calculation")
    # Check clustering data
    sc.pl.draw_graph(oracle.adata, color="louvain_annot")
    # Calculate GRN for each population in "louvain_annot" clustering unit.
    # This step may take some time.(~30 minutes)
    links = oracle.get_links(cluster_name_for_GRN_unit="louvain_annot",
                             alpha=10,
                             verbose_level=1)
    # (Optional) Export GRNs
    links.links_dict.keys()

    if use_reduced_data:
        pass
    else:
        # 5.3. (Optional) Change order
        # Show the contents of pallete
        links.palette
        # Change the order of pallete
        order = ['MEP_0', 'Mk_0', 'Ery_0',
                 'Ery_1', 'Ery_2', 'Ery_3', 'Ery_4', 'Ery_5',
                 'Ery_6', 'Ery_7', 'Ery_8', 'Ery_9',
                 'GMP_0', 'GMP_1', 'GMP_2', 'GMPl_0', 'GMPl_1',
                 'Mo_0', 'Mo_1', 'Mo_2',
                 'Gran_0', 'Gran_1', 'Gran_2', 'Gran_3']
        links.palette = links.palette.loc[order]
        links.palette
    # Save Links object.
    links.to_hdf5(file_path="links.celloracle.links")

    # 6. Network preprocessing
    print("\n6. Network preprocessing")
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

    print("\n")
    plt.rcParams["figure.figsize"] = [9, 4.5]
    links.plot_degree_distributions(plot_model=True,
                                    save=f"{save_folder}/degree_distribution/",
                                    )
    plt.rcParams["figure.figsize"] = [6, 4.5]
    # Calculate network scores.
    links.get_network_score()
    links.merged_score.head()
    # Save Links object.
    links.to_hdf5(file_path="links.celloracle.links")
    # You can load files with the following command.
    links = co.load_hdf5(file_path="links.celloracle.links")

    # 7. Network analysis; Network score for each gene
    print("\n7. Network analysis; Network score for each gene")
    # Check cluster name
    links.cluster
    # Visualize top n-th genes with high scores.
    links.plot_scores_as_rank(cluster="MEP_0", n_gene=30,
                              save=f"{save_folder}/ranked_score")
    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="eigenvector_centrality",
                                   cluster1="MEP_0", cluster2="GMPl_0",
                                   percentile=98,
                                   save=f"{save_folder}/score_comparison")
    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="betweenness_centrality",
                                   cluster1="MEP_0", cluster2="GMPl_0",
                                   percentile=98,
                                   save=f"{save_folder}/score_comparison")
    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="degree_centrality_all",
                                   cluster1="MEP_0", cluster2="GMPl_0",
                                   percentile=98,
                                   save=f"{save_folder}/score_comparison")
    # Visualize Gata2 network score dynamics
    links.plot_score_per_cluster(goi="Gata2",
                                 save=f"{save_folder}/network_score_per_gene/")
    # Visualize Cebpa network score dynamics
    links.plot_score_per_cluster(goi="Cebpa")
    #
    

    # 8. Network analysis; network score distribution
    print("\n8. Network analysis; network score distribution")
    plt.rcParams["figure.figsize"] = [6, 4.5]
    # Plot degree_centrality
    plt.subplots_adjust(left=0.15, bottom=0.3)
    plt.ylim([0,0.040])
    links.plot_score_discributions(values=["degree_centrality_all"],
                                   method="boxplot",
                                   save=f"{save_folder}",
                                  )
    # Plot eigenvector_centrality
    plt.subplots_adjust(left=0.15, bottom=0.3)
    plt.ylim([0, 0.28])
    links.plot_score_discributions(values=["eigenvector_centrality"],
                                   method="boxplot",
                                   save=f"{save_folder}")
    ## Distribution of netowrk entropy
    plt.subplots_adjust(left=0.15, bottom=0.3)
    links.plot_network_entropy_distributions(save=f"{save_folder}")

    # 6. Remove files
    print("Remove files")
    for i in ["TF_data_in_Paul15.csv", "Paul_15_data.celloracle.oracle",
              "links.celloracle.links"]:
        os.remove(i)
    os.system("rm -r figures")

    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
