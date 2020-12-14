# -*- coding: utf-8 -*-

import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import scanpy as sc
from sklearn.neighbors import KNeighborsRegressor

from ..trajectory.oracle_utility import _adata_to_color_dict

class Pseudotime_calculator():
    """
    This is a custom object for pseudotime calculation.

    """
    def __init__(self, oracle_object=None, adata=None, obsm_key=None, cluster_column_name=None, lineage_dictionary=None):
        self.adata = None
        self.adata_raw = None
        self.obsm_key = None
        self.cluster_column_name = None
        self.lineage_dictionary = None

        if adata is not None:
            if (obsm_key is not None) & (cluster_column_name is not None):
                self.load_anndata(adata=adata, obsm_key=obsm_key, cluster_column_name=cluster_column_name)
            else:
                print("Please set obsm_key and cluster_column_name to instantiate Pseudotime_calculator with anndata.")

        if oracle_object is not None:
            self.load_oracle_object(oracle_object=oracle_object)

        if lineage_dictionary is not None:
            self.set_lineage(lineage_dictionary=lineage_dictionary)

    def load_anndata(self, adata, obsm_key, cluster_column_name):
        self.adata = adata.copy()
        self.adata_raw = adata.copy()
        self.obsm_key = obsm_key
        self.cluster_column_name = cluster_column_name
        self.cluster_list = sorted(list(adata.obs[cluster_column_name].unique()))

    def load_oracle_object(self, oracle_object):
        self.load_anndata(adata=oracle_object.adata,
                          obsm_key=oracle_object.embedding_name,
                          cluster_column_name=oracle_object.cluster_column_name)

    def reset(self):
        self.adata = self.adata_raw.copy()

    def set_lineage(self, lineage_dictionary):
        self.lineage_dictionary = lineage_dictionary.copy()

        for lineage, clusters in lineage_dictionary.items():
            self.adata.obs[lineage] = self.adata.obs[self.cluster_column_name].isin(clusters)
            self.adata.obs[lineage] = self.adata.obs[lineage].astype("str").astype("category")
            self.adata.uns[f"{lineage}_colors"] = ["#D0D3D4", "#EC7063"]

    def re_calculate_diffmap(self, n_neighbors=50, use_rep=None):
        del self.adata.obsm["X_diffmap"]
        if use_rep is not None:
            use_rep = self.obsm_key
        sc.pp.neighbors(self.adata, n_neighbors=50, use_rep=use_rep)
        sc.tl.diffmap(self.adata)

    def set_root_cells(self, root_cells):
        """
        Args:
            root_cells (dictionary): key is the name of lineage, the value is name of cell
        """
        self.root_cells = root_cells.copy()

        # check if the format is correct
        for lineage, root_cell in root_cells.items():
            if lineage not in self.adata.obs.columns:
                print(f"The key, {lineaeg} is not in the lineage data. Make sure the root_cells has correct format.")
            if root_cell not in self.adata.obs.index[self.adata.obs[lineage] == "True"]:
                print(f"The cell, {root_cell} in the {lineage} was not found in your adata. Please make sure the format is correct.")

    def get_pseudotime_per_each_lineage(self, root_cells=None, plot=False):

        if root_cells is not None:
            self.set_root_cells(root_cells=root_cells)

        pseudotime_estimation_for_each_lineage(adata=self.adata, root_cells=self.root_cells)

        # Impute inf or nan value with real value
        main_pseudotime = self.adata.obs["Pseudotime"]
        main_pseudotime = main_pseudotime.replace(np.inf, np.nan)

        if (main_pseudotime.isna().sum() > 0):
            if main_pseudotime.isna().mean() >= 0.5:
                raise ValueError("Found too many cells that failed to get pseudotime. Please check lineage dictionary.")
            else:
                # Replace na and inf with NN's value
                y_filled =  _fill_inf_and_na(X=self.adata.obsm[self.obsm_key],
                                    y_with_inf=self.adata.obs["Pseudotime"])
                self.adata.obs["Pseudotime"] = y_filled

        if plot:
            #sc.pl.umap(self.adata, color=[i for i in self.adata.obs.columns if "Pseudotime" in i], ncols=3, size=20, cmap="rainbow")
            self.plot_pseudotime()

    def plot_cluster(self, s=10, fontsize=8):

        embedding = self.adata.obsm[self.obsm_key]

        color_dict = _adata_to_color_dict(adata=self.adata, cluster_use=self.cluster_column_name)

        fig = plt.figure()
        for cluster, color in color_dict.items():
            idx = np.where(self.adata.obs[self.cluster_column_name] == cluster)[0]
            plt.scatter(embedding[idx, 0], embedding[idx, 1], color=color, s=s, label=cluster)

        plt.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0, fontsize=fontsize)
        plt.title(self.cluster_column_name)
        plt.axis("off")
        plt.show()

    def plot_lineages(self, s=10):

        embedding = self.adata.obsm[self.obsm_key]

        for lineage in self.lineage_dictionary.keys():
            color = ["#EC7063" if i=="True" else "#D0D3D4" for i in self.adata.obs[lineage].values]

            fig = plt.figure()
            plt.scatter(embedding[:, 0], embedding[:, 1], color=color, s=s)
            plt.title(lineage)
            plt.axis("off")
            plt.show()


    def plot_root_cells(self, s=10):

        embedding = self.adata.obsm[self.obsm_key]

        for lineage, root_cell in self.root_cells.items():
            color = ["#EC7063" if i=="True" else "#D0D3D4" for i in self.adata.obs[lineage].values]
            root_cell_idx = np.where(self.adata.obs.index == root_cell)[0]
            fig = plt.figure()
            plt.scatter(embedding[:, 0], embedding[:, 1], color=color, s=s)
            plt.scatter(embedding[root_cell_idx, 0], embedding[root_cell_idx, 1], color="black", s=s*10, label="root_cell")
            plt.legend()
            plt.title(lineage)
            plt.axis("off")
            plt.show()


    def plot_pseudotime(self, s=10, cmap="rainbow"):

        embedding = self.adata.obsm[self.obsm_key]

        pseudotime_list = [f"Pseudotime_{i}" for i in self.root_cells.keys()] + ["Pseudotime"]
        for pseudotime in pseudotime_list:

            color = self.adata.obs[pseudotime]
            idx = np.where(~color.isna())[0]

            fig = plt.figure()
            plt.scatter(embedding[idx, 0], embedding[idx, 1], c=color[idx], s=s, cmap=cmap)
            plt.title(pseudotime)
            plt.axis("off")
            plt.show()


def pseudotime_estimation_for_each_lineage(adata, root_cells):
    """
    Calculate pseudotime using dpt method in scanpy. https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dpt.html#scanpy.tl.dpt
    This function consists of four steps below.

    (1) Split scRNA-seq data into subgroup (lineage).
    (2) Calculate pseudotime by dpt
    (3) Re-aggregate cells into original adata
    (4) Fill inf and nan value in the pseudotime by KNN regression.

    Args:
        adata (anndata): scRNA-seq data
        root_cells (dictionary): The key is the name of lineage. This information should exist in the adata.obs
            The value of this dictionary is the name of root cell.

    Returns:
        None: Results will be stored in the adata.obs: (i)The pseudotime for each lineage and (ii)the pseudotime for whole population.

    """

    pseudotime_columns = []
    for lineage, root_cell in root_cells.items():
        # Make subset of adata for single lineage
        bool_ = adata.obs[lineage].values

        if isinstance(bool_[0], str):
            bool_ = [{"True": True, "False": False}[i] for i in bool_]
        cells_in_the_lineage = adata.obs.index[bool_].values
        adata_ = adata[cells_in_the_lineage, :]
        # Set root
        adata_.uns['iroot'] = np.where(adata_.obs.index == root_cell)[0][0]

        # Calculate pseudotime
        sc.tl.dpt(adata_)

        # Save calculated data into original adata
        adata.obs[f'Pseudotime_{lineage}'] = np.nan
        adata.obs.loc[adata_.obs.index.values, f'Pseudotime_{lineage}'] = adata_.obs["dpt_pseudotime"]

        pseudotime_columns.append(f'Pseudotime_{lineage}')

    # Calculate pseudotime average
    adata.obs["Pseudotime"] = adata.obs[pseudotime_columns].mean(axis=1, skipna=True)



def _fill_inf_and_na(X, y_with_inf):
    # Make bool list for the cells to be replaced
    na_cell_bool = y_with_inf.replace(np.inf, np.nan).isna().values

    # Make KNN model and fitting, prediction
    knn = KNeighborsRegressor(n_neighbors=30)
    knn.fit(X[~na_cell_bool, :], y_with_inf[~na_cell_bool])
    y_filled = y_with_inf.copy()
    y_filled[na_cell_bool] = knn.predict(X[na_cell_bool, :])

    # check
    assert(y_filled.replace(np.inf, np.nan).isna().sum() == 0)

    return y_filled
