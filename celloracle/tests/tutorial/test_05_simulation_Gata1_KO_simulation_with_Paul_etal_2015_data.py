import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


import celloracle as co
print("co version: ", co.__version__)

plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams["savefig.dpi"] = 600


def main():
    test_tutorial_simulation()
    
def test_tutorial_simulation():
    # Start analysis
    import time
    start = time.time()

    # Make folder to save plots
    save_folder = "figures"
    os.makedirs(save_folder, exist_ok=True)

    # 1. Load data
    print("\n1. Load data")
    oracle = co.data.load_tutorial_oracle_object()
    oracle
    # Here, we load demo links object for the training purpose.
    links = co.data.load_tutorial_links_object()

    # 2. Make predictive models for simulation
    print("\n2. Make predictive models for simulation")
    links.filter_links()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True,
                                  verbose_level=1)

    # 3. In silico TF perturbation analysis
    print("\n3. In silico TF perturbation analysis")
    # Check gene expression
    goi = "Gata1"
    sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                     layer="imputed_count", use_raw=False, cmap="viridis")
    # Plot gene expression in histogram
    sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
    plt.show()
    # Enter perturbation conditions to simulate signal propagation after the perturbation.
    oracle.simulate_shift(perturb_condition={goi: 0.0},
                          n_propagation=3)
    # Get transition probability
    oracle.estimate_transition_prob(n_neighbors=200,
                                    knn_random=True,
                                    sampled_fraction=1)
    # Calculate embedding
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # 4. Visualization
    print("\n4. Visualization")
    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    scale = 25
    # Show quiver plot
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.show()
    # Vector field graph visualization
    # n_grid = 40 is a good starting value.
    n_grid = 40
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    # Search for best min_mass.
    oracle.suggest_mass_thresholds(n_suggestion=12)
    min_mass = 0.01
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
    # Plot vector fields
    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    scale_simulation = 0.5
    # Show quiver plot
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.show()
    # Plot vector field with cell cluster
    fig, ax = plt.subplots(figsize=[8, 8])
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation,
                                        ax=ax, show_background=False)

    # 5. Compare simulation vector with development vectors
    print("\n5. Compare simulation vector with development vectors")
    # Prepare pseudotime data
    # Visualize pseudotime
    fig, ax = plt.subplots(figsize=[6,6])
    sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name,
                    ax=ax, cmap="rainbow",
                    color=["Pseudotime"])
    # Calculation with Gradient calculator object
    from celloracle.applications import Gradient_calculator
    # Instantiate Gradient calculator object
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
    gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3},
                                     plot=True)
    # Calculate graddient
    gradient.calculate_gradient()
    # Show results
    scale_dev = 40
    gradient.visualize_results(scale=scale_dev, s=5)
    # Visualize results
    fig, ax = plt.subplots(figsize=[6, 6])
    gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)
    # Save gradient object if you want.
    gradient.to_hdf5("Paul_etal.celloracle.gradient")
    # Calculate PS
    from celloracle.applications import Oracle_development_module
    # Make Oracle_development_module to compare two vector field
    dev = Oracle_development_module()
    # Load development flow
    dev.load_differentiation_reference_data(gradient_object=gradient)
    # Load simulation result
    dev.load_perturb_simulation_data(oracle_object=oracle)
    # Calculate inner produc scores
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)
    # Let's visualize the results
    dev.visualize_development_module_layout_0(s=5,
                                              scale_for_simulation=scale_simulation,
                                              s_grid=50,
                                              scale_for_pseudotime=scale_dev,
                                              vm=0.02)
    # Show perturbation scores
    fig, ax = plt.subplots(figsize=[6, 6])
    dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)
    # Show perturbation scores with perturbation simulation vector field
    fig, ax = plt.subplots(figsize=[6, 6])
    dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)
    dev.plot_simulation_flow_on_grid(scale=scale_simulation,
                                     show_background=False, ax=ax)

    # 6. Focus on a single development lineage to interpret the results in detail
    print("\n6. Focus on a single development lineage to interpret the results in detail")
    # Get cell index list for the cells of interest
    clusters = ['Ery_0', 'Ery_1', 'Ery_2', 'Ery_3', 'Ery_4', 'Ery_5', 'Ery_6',
                'Ery_7', 'Ery_8', 'Ery_9', 'MEP_0', 'Mk_0']
    cell_idx = np.where(oracle.adata.obs["louvain_annot"].isin(clusters))[0]
    # Check index
    print(cell_idx)
    # Calculate PS
    dev = Oracle_development_module()
    # Load development flow
    dev.load_differentiation_reference_data(gradient_object=gradient)

    # Load simulation result
    dev.load_perturb_simulation_data(oracle_object=oracle,
                                 cell_idx_use=cell_idx, # Enter cell id list
                                 name="Lineage_MEP" # Name of this cell group. You can enter any name.
                                 )
    # Calculation
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)
    # Let's visualize the results
    dev.visualize_development_module_layout_0(s=5,
                                              scale_for_simulation=scale_simulation,
                                              s_grid=50,
                                              scale_for_pseudotime=scale_dev,
                                              vm=0.03)
    # OLD VISUALIZATION
    from celloracle.visualizations.config import CONFIG
    CONFIG["cmap_ps"] = "coolwarm"
    dev.visualize_development_module_layout_0(s=5,
                                              scale_for_simulation=scale_simulation,
                                              s_grid=50,
                                              scale_for_pseudotime=scale_dev,
                                              vm=0.03)

    # Remove files
    print("Remove files")
    for i in ["Paul_etal.celloracle.gradient"]:
        os.remove(i)
        #pass
    os.removedirs("figures")
    elapsed = time.time() - start

    print("Success")
    print(f"Total time: {elapsed/60:.3g} minutes")


if __name__ == "__main__":
    main()
