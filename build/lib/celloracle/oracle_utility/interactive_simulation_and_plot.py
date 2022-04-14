# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import h5py


from .utility import Oracle_data_strage
from .development_analysis import Oracle_development_module

from ..trajectory.oracle_core import Oracle

#import scanpy as sc
#import seaborn as sns
#from ipywidgets import interactive
import h5py


DEFAULT_PARAMETERS = {"n_neighbors": 200,
                      "quiver_scale": 30, "quiver_scale_grid": 1.5,
                      "min_mass": 1.0, "n_grid": 40,
                      "min_magnitude": None,
                      "sampled_fraction": 0.5,
                      "n_propagation": 3,
                      "n_steps_for_mc": 20,
                      "gene": "Gata1"}


COMMON_ATTRS = ["embedding", "colorandum", "cluster_column_name"]


class Oracle_extended(Oracle_data_strage, Oracle_development_module):

    def __init__(self, oracle, hdf_path, mode, obsm_key="X_umap"):
        self.oracle = oracle
        self.gene = None
        self.n_neighbors = None
        self.n_grid = None
        self.names = []
        self.obsm_key = obsm_key

        self.set_hdf_path(path=hdf_path, create_if_not_exist=True)

        if mode == "write":
            self.oracle.embedding = self.oracle.adata.obsm[self.oracle.embedding_name]
            self.save_data(oracle=self.oracle, place=f"common", attributes=COMMON_ATTRS)

        elif mode == "read":
            self.load_data(oracle=self.oracle, place=f"common", attributes=COMMON_ATTRS)
            self.oracle.cluster_column_name = str(self.oracle.cluster_column_name)

    # Module for perturb simulation
    def load_simulation(self, gene=DEFAULT_PARAMETERS["gene"]):
        print("Loading data ..")
        self.load_data(oracle=self.oracle, place=f"simulation/{gene}",
                       attributes=["delta_embedding", "delta_embedding_random"])

        self.load_dfs(oracle=self.oracle, place=f"mcmc/{gene}",
                      attributes=["mcmc_transition", "mcmc_transition_random"])

        self.oracle.corrcoef_random = "dummy"
        self.gene = gene
        print(f"Data loaded. Perturbed gene: {gene}")

    def interactive_simulation(self, gene,
                               n_neighbors=DEFAULT_PARAMETERS["n_neighbors"],
                               sampled_fraction=DEFAULT_PARAMETERS["sampled_fraction"],
                               n_propagation=DEFAULT_PARAMETERS["n_propagation"],
                               n_steps_for_mc=DEFAULT_PARAMETERS["n_steps_for_mc"],
                               save=True):
        """


        """
        print("Calculating ...")
        # Enter perturbation conditions to simulate signal propagation after the perturbation.
        self.oracle.simulate_shift(perturb_condition={gene: 0},
                                   ignore_warning=True,
                                   n_propagation=n_propagation)

        # Get transition probability
        self.oracle.estimate_transition_prob(n_neighbors=n_neighbors, knn_random=True, sampled_fraction=sampled_fraction)

        # Calculate embedding
        self.oracle.calculate_embedding_shift(sigma_corr = 0.05)

        # Do Markov simulation
        self.oracle.run_markov_chain_simulation(n_steps=n_steps_for_mc, n_duplication=1)
        self.oracle.get_mcmc_cell_transition_table(end=n_steps_for_mc)


        self.gene = gene
        self.n_neighbors = n_neighbors

        if save:
            self.save_data(oracle=self.oracle, place=f"simulation/{gene}",
                           attributes=["delta_embedding", "delta_embedding_random"])


            self.save_dfs(oracle=self.oracle, place=f"mcmc/{gene}",
                          attributes=["mcmc_transition", "mcmc_transition_random"])

            print("Results were saved to hdf5 file.")

        print("finished")


    def interactive_plot_quiver(self, quiver_scale=DEFAULT_PARAMETERS["quiver_scale"]):

        plt.figure(None,(6,6))
        plt.scatter(self.oracle.embedding[:, 0], self.oracle.embedding[:, 1],
                c="0.8", alpha=0.2, s=38, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

        quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                       linewidths=0.35, width=0.0045,edgecolors="k",
                       color=self.oracle.colorandum, alpha=1)
        plt.quiver(self.oracle.embedding[:, 0], self.oracle.embedding[:, 1],
                   self.oracle.delta_embedding[:, 0], self.oracle.delta_embedding[:, 1],
                   scale=quiver_scale, **quiver_kwargs)

        plt.title(f"Knockout simulation: {self.gene}")
        plt.axis("off")

    def interactive_plot_grid(self,
                              quiver_scale=DEFAULT_PARAMETERS["quiver_scale_grid"],
                              min_mass=DEFAULT_PARAMETERS["min_mass"],
                              #min_magnitude=DEFAULT_PARAMETERS["min_magnitude"],
                              n_grid=DEFAULT_PARAMETERS["n_grid"],
                              plot_random=False):
        # Plot whole graph
        if plot_random:
            plt.figure(None,(13,6))
        else:
            plt.figure(None,(6,6))

        if self.n_grid is None:
            self.oracle.calculate_grid_arrows(smooth=0.8, steps=(n_grid, n_grid), n_neighbors=200)
        else:
            if n_grid != self.n_grid:
                self.oracle.calculate_grid_arrows(smooth=0.8, steps=(n_grid, n_grid), n_neighbors=200)

        plt.title(f"Perturb simulation: {self.gene}")

        self.oracle.plot_grid_arrows(quiver_scale=quiver_scale,
                                scatter_kwargs_dict={"alpha":0.35, "lw":0.35,
                                                      "edgecolor":"0.4", "s":38,
                                                      "rasterized":True},
                                min_mass=min_mass,
                                min_magnitude=None,
                                angles='xy', scale_units='xy',
                                headaxislength=2.75,
                                headlength=5, headwidth=4.8, minlength=1.5,
                                plot_random=plot_random, scale_type="relative")




    def _plot_quiver_for_a_cluster(self, cluster_column_name, cluster, quiver_scale, plot_whole_cells=False):

        ix_choice = _get_ix_for_a_cluster(self.oracle, cluster_column_name, cluster)

        if plot_whole_cells:

            plt.scatter(self.oracle.embedding[:, 0], self.oracle.embedding[:, 1],
                        c="lightgray", alpha=1, s=38, lw=0.3, rasterized=True)

        plt.scatter(self.oracle.embedding[ix_choice, 0], self.oracle.embedding[ix_choice, 1],
                    c="0.8", alpha=0.2, s=38, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)



        quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                           linewidths=0.35, width=0.0045,edgecolors="k",
                           color=self.oracle.colorandum[ix_choice], alpha=1)
        plt.quiver(self.oracle.embedding[ix_choice, 0], self.oracle.embedding[ix_choice, 1],
                   self.oracle.delta_embedding[ix_choice, 0],
                   self.oracle.delta_embedding[ix_choice, 1],
                   scale=quiver_scale, **quiver_kwargs)

        plt.axis("off")

    def interactive_plot_quiver_for_a_cluster(self, cluster_column_name, cluster,
                                              quiver_scale=DEFAULT_PARAMETERS["quiver_scale"]):
        plt.figure(figsize=(13, 6))

        plt.title(f"Knockout simulation: {self.gene}")
        plt.subplot(1, 2, 1)
        self._plot_quiver_for_a_cluster(cluster_column_name=cluster_column_name,
                                        cluster=cluster, quiver_scale=quiver_scale,
                                        plot_whole_cells=False)

        plt.subplot(1, 2, 2)
        self._plot_quiver_for_a_cluster(cluster_column_name=cluster_column_name,
                                        cluster=cluster, quiver_scale=quiver_scale,
                                        plot_whole_cells=True)

    # Development analysis module
    def save_development_analysis_results(self, gene, cluster_column_name, cluster):

        self.save_data(oracle=self.oracle_dev,
                       place=f"dev_analysis/{cluster_column_name}/{cluster}/{gene}",
                       attributes=["embedding", "pseudotime", "mass_filter",
                                   "flow_grid", "flow", "flow_norm_rndm",
                                   "new_pseudotime", "gradient", "inner_product", "stage", "stage_grid"])

        self.save_dfs(oracle=self.oracle_dev,
                      place=f"inner_product/{cluster_column_name}/{cluster}/{gene}",
                      attributes=["inner_product_stats", "inner_product_stats_grouped"])


        print("Results were saved to hdf5 file.")

    def load_development_analysis_results(self, gene, cluster_column_name, cluster):

        if not hasattr(self, "oracle_dev"):
            setattr(self, "oracle_dev", Oracle())

        print("Loading data ..")
        self.load_data(oracle=self.oracle_dev,
                       place=f"dev_analysis/{cluster_column_name}/{cluster}/{gene}",
                       attributes=["embedding", "pseudotime", "mass_filter",
                                   "flow_grid", "flow", "flow_norm_rndm",
                                   "new_pseudotime", "gradient",
                                   "inner_product", "stage", "stage_grid"])

        self.load_dfs(oracle=self.oracle_dev,
                      place=f"inner_product/{cluster_column_name}/{cluster}/{gene}",
                      attributes=["inner_product_stats", "inner_product_stats_grouped"])

        self.gene = gene
        self.oracle_dev.cluster_loaded = cluster
        self.oracle_dev.cluster_column_name_loaded = cluster_column_name

        print(f"Data loaded. Gene: {gene}")


def _get_ix_for_a_cluster(oracle, cluster_column_name, cluster):
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix
