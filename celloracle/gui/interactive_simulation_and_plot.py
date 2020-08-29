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

#import scanpy as sc
#import seaborn as sns
#from ipywidgets import interactive

DEFAULT_PARAMETERS = {"n_neighbors": 200,
                      "quiver_scale": 30, "quiver_scale_grid": 2.0,
                      "min_mass": 0.015, "n_grid": 40,
                      "gene": "Gata1"}


class Interactive_celloracle_simulator():
    def __init__(self, oracle=None, ods=None):
        self.oracle = oracle
        self.ods = ods
        self.gene = None
        self.n_neighbors = None
        self.n_grid = None




    def load_calculated_simulation(self, gene=DEFAULT_PARAMETERS["gene"]):
        print("Loading data ..")
        self.ods.load_one_data(gene)

        self.ods.get_back_data(self.oracle, gene)
        self.oracle.embedding = self.oracle.adata.obsm[self.oracle.embedding_name]
        self.gene = gene
        print("Done")

    def interactive_simulation(self, goi, n_neighbors=DEFAULT_PARAMETERS["n_neighbors"]):
        """


        """
        print("Calculating ...")
        # Enter perturbation conditions to simulate signal propagation after the perturbation.
        self.oracle.simulate_shift(perturb_condition={goi: 0},
                              n_propagation=3)

        # Get transition probability
        self.oracle.estimate_transition_prob(n_neighbors=n_neighbors, knn_random=True, sampled_fraction=0.5)

        # Calculate embedding
        self.oracle.calculate_embedding_shift(sigma_corr = 0.05)

        self.gene = goi
        self.n_neighbors = n_neighbors

        print("finished")


    def interactive_plot_quiver(self, quiver_scale=DEFAULT_PARAMETERS["quiver_scale"]):


        embedding = self.oracle.adata.obsm[self.oracle.embedding_name]

        plt.figure(None,(6,6))
        plt.scatter(embedding[:, 0], embedding[:, 1],
                c="0.8", alpha=0.2, s=38, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

        quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                       linewidths=0.35, width=0.0045,edgecolors="k",
                       color=self.oracle.colorandum, alpha=1)
        plt.quiver(embedding[:, 0], embedding[:, 1],
                   self.oracle.delta_embedding[:, 0], self.oracle.delta_embedding[:, 1],
                   scale=quiver_scale, **quiver_kwargs)

        plt.title(f"Knockout simulation: {self.gene}")
        plt.axis("off")

    def interactive_plot_grid(self,
                              quiver_scale=DEFAULT_PARAMETERS["quiver_scale_grid"],
                              min_mass=DEFAULT_PARAMETERS["min_mass"],
                              n_grid=DEFAULT_PARAMETERS["n_grid"],
                              plot_random=False):
        # Plot whole graph
        if plot_random:
            plt.figure(None,(13,6))
        else:
            plt.figure(None,(6,6))

        if self.n_grid is None:
            self.oracle.calculate_grid_arrows(smooth=0.8, steps=(n_grid, n_grid), n_neighbors=300)
        else:
            if n_grid != self.n_grid:
                self.oracle.calculate_grid_arrows(smooth=0.8, steps=(n_grid, n_grid), n_neighbors=300)

        plt.title(f"Knockout simulation: {self.gene}")
        self.oracle.plot_grid_arrows(quiver_scale=quiver_scale,
                                scatter_kwargs_dict={"alpha":0.35, "lw":0.35,
                                                      "edgecolor":"0.4", "s":38,
                                                      "rasterized":True},
                                min_mass=min_mass, angles='xy', scale_units='xy',
                                headaxislength=2.75,
                                headlength=5, headwidth=4.8, minlength=1.5,
                                plot_random=plot_random, scale_type="relative")



    def _plot_quiver_for_a_cluster(self, cluster, quiver_scale, plot_whole_cells=False):

        ix_choice = _get_ix_for_a_cluster(self.oracle, cluster)

        embedding = self.oracle.adata.obsm[self.oracle.embedding_name]

        if plot_whole_cells:

            plt.scatter(embedding[:, 0], embedding[:, 1],
                        c="lightgray", alpha=1, s=38, lw=0.3, rasterized=True)

        plt.scatter(embedding[ix_choice, 0], embedding[ix_choice, 1],
                    c="0.8", alpha=0.2, s=38, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)



        quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                           linewidths=0.35, width=0.0045,edgecolors="k",
                           color=self.oracle.colorandum[ix_choice], alpha=1)
        plt.quiver(embedding[ix_choice, 0], embedding[ix_choice, 1],
                   self.oracle.delta_embedding[ix_choice, 0],
                   self.oracle.delta_embedding[ix_choice, 1],
                   scale=quiver_scale, **quiver_kwargs)

        plt.axis("off")

    def interactive_plot_quiver_for_a_cluster(self, cluster,
                                              quiver_scale=DEFAULT_PARAMETERS["quiver_scale"]):
        plt.figure(figsize=(13, 6))

        plt.title(f"Knockout simulation: {self.gene}")
        plt.subplot(1, 2, 1)
        self._plot_quiver_for_a_cluster(cluster, quiver_scale, plot_whole_cells=False)

        plt.subplot(1, 2, 2)
        self._plot_quiver_for_a_cluster(cluster, quiver_scale, plot_whole_cells=True)


def _get_ix_for_a_cluster(oracle, cluster, cluster_column_name=None):
    if cluster_column_name is None:
        cluster_column_name = oracle.cluster_column_name
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix
