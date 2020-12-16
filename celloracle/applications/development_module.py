# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys
from copy import deepcopy

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
#import h5py

from scipy.stats import wilcoxon

from ..trajectory.oracle_core import Oracle
from .scatter_to_grid import scatter_value_to_grid_value
from .utility import Data_strage

from ..visualizations.config import CONFIG
from ..visualizations.development_module_visualization import (\
                            plot_cluster_whole,
                            plot_cluster_cells_use,
                            plot_background,
                            plot_pseudotime,
                            plot_pseudotime_on_grid,
                            plot_reference_flow_on_grid,
                            plot_simulation_flow_on_grid,
                            plot_inner_product_on_grid,
                            plot_inner_product_on_pseudotime,
                            plot_inner_product_as_box,
                            plot_quiver,
                            visualize_development_module_layout_0,
                            visualize_development_module_layout_1,
                            visualize_development_module_layout_2)

class Oracle_development_module(Data_strage):

    def __init__(self, oracle_object=None, gradient_object=None, name=None):
        super().__init__()
        self._exemptions_when_del_attrs = []

        self.name = name
        self.cell_idx_use = None

        if gradient_object is not None:
            self.load_differentiation_reference_data(gradient_object=gradient_object)

        if oracle_object is not None:
            self.load_perturb_simulation_data(oracle_object=oracle_object)

    def copy(self):
        return deepcopy(self)

    def set_hdf_path(self, path):
        key = "" # placeholder
        self._set_hdf_path(path=path, key=key, create_if_not_exist=True)

    def dump_hdf5(self, gene, misc):

        if hasattr(self, "_path"):
            pass
        else:
            raise ValueError("Run set_hdf_path to set file path first.")

        self._dump_hdf5(place=f"{gene}/{misc}")


    def load_hdf5(self, gene, misc, specify_attributes=None):

        self._load_hdf5(place=f"{gene}/{misc}", specify_attributes=specify_attributes)

    def del_attrs(self, exemptions=[]):
        # Delete all attributes
        self._del_attrs(exemptions=exemptions)

    def get_hdf5_info(self):

        dic = {}
        # Gene list
        dic["gene_list"] = np.unique([i.split("/")[0] for i in self._names])
        # Misc list
        dic["misc_list"] = np.unique([i.split("/")[1] for i in self._names if len(i.split("/")) >= 2])
        #
        unique_keys = np.unique(["/".join(i.split("/")[:2]) for i in self._names if len(i.split("/")) >= 2])
        dic["gene_misc_lists"] = list(map(lambda x: x.split("/"), unique_keys))

        # Get dictionary; key = misc, value = list of gene
        misc_gene_dictionary = {}
        for i in dic["misc_list"]:
            misc_gene_dictionary[i] = [gene for gene, misc in dic["gene_misc_lists"] if misc == i]
        dic["misc_gene_dictionary"] = misc_gene_dictionary

        return dic


    def load_differentiation_reference_data(self, gradient_object):

        self.n_grid = gradient_object.n_grid
        self.min_mass = gradient_object.min_mass
        self.smooth = gradient_object.smooth
        self.n_neighbors = gradient_object.n_neighbors

        self.pseudotime = gradient_object.pseudotime.copy() # shape = (n_cell, )
        self.pseudotime_on_grid = gradient_object.pseudotime_on_grid.copy() # shape = (n_grid*n_grid, )
        self.mass_filter_reference = gradient_object.mass_filter.copy() # shape = (n_grid*n_grid, )
        self.mass_filter_whole_reference = gradient_object.mass_filter_whole.copy() # shape = (n_grid*n_grid, )
        self.gridpoints_coordinates = gradient_object.gridpoints_coordinates.copy() # shape = (n_grid*n_grid, 2)

        self.ref_flow = gradient_object.ref_flow.copy() # shape = (n_grid*n_grid, 2)

        if gradient_object.cell_idx_use is not None:
            self.cell_idx_use = gradient_object.cell_idx_use.copy()

    def load_perturb_simulation_data(self, oracle_object, cell_idx_use=None, name=None, min_mass=None, n_neighbors=None, smooth=None):

        if not hasattr(self, "pseudotime"):
            raise ValueError("Please load differentiation reference data first.")

        # 0. Store data
        for i in ["delta_embedding", "embedding", "colorandum"]:
            setattr(self, i, getattr(oracle_object, i))

        if cell_idx_use is not None:
            self.cell_idx_use = np.array(cell_idx_use)

        if name is not None:
            self.name = name

        # 1. Subset oracle data if cell_idx_use is specified
        if cell_idx_use is None:
            oracle_object_ = oracle_object
        else:
            oracle_object_ = subset_oracle_for_development_analysiis(oracle_object=oracle_object,
                                                          cell_idx_use=cell_idx_use)

        # 2. Grid calculation
        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        if smooth is None:
            smooth = self.smooth
        x_min, y_min = self.embedding.min(axis=0)
        x_max, y_max = self.embedding.max(axis=0)
        xylim = ((x_min, x_max), (y_min, y_max))

        oracle_object_.calculate_grid_arrows(smooth=smooth, steps=(self.n_grid, self.n_grid), n_neighbors=n_neighbors, xylim=xylim)


        # 3. Store result of grid calculation
        for i in ["flow", "flow_rndm"]:
            setattr(self, i, getattr(oracle_object_, i))

        if not (oracle_object_.flow_grid == self.gridpoints_coordinates).all():
            raise ValueError("Grid point cordinates are not mached.")

        # 4. mass_filter calculation
        if min_mass is None:
            min_mass = self.min_mass
        self.mass_filter_simulation = (oracle_object_.total_p_mass < min_mass)


    def calculate_inner_product(self):

        # Calculate inner product between the pseudotime-gradient and the perturb-gradient
        self.inner_product = np.array([np.dot(i, j) for i, j in zip(self.flow, self.ref_flow)])


    def calculate_digitized_ip(self, n_bins=10):

        inner_product_df = pd.DataFrame({"score": self.inner_product[~self.mass_filter_simulation],
                                                "pseudotime": self.pseudotime_on_grid[~self.mass_filter_simulation]})


        bins = _get_bins(inner_product_df.pseudotime, n_bins)
        inner_product_df["pseudotime_id"] = np.digitize(inner_product_df.pseudotime, bins) - 1

        self.inner_product_df = inner_product_df

    def get_p_inner_product(self, method="wilcoxon"):
        return get_p_inner_product(inner_product_df=self.inner_product_df, method=method)

    def get_sum_of_negative_ips(self):
        return get_sum_of_negative_ips(inner_product_df=self.inner_product_df)

    def get_sum_of_positive_ips(self):
        return get_sum_of_positive_ips(inner_product_df=self.inner_product_df)

    # Visualization
    def plot_cluster_whole(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):
        plot_cluster_whole(self=self, ax=ax, s=s, args=args)

    def plot_cluster_cells_use(self, ax=None, s=CONFIG["s_scatter"], color=None, show_background=True, args=CONFIG["default_args"]):
        plot_cluster_cells_use(self=self, ax=ax, s=s, color=color, show_background=show_background, args=args)

    def plot_background(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):
        plot_background(self=self, ax=ax, s=s, args=args)

    def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_pseudotime(self=self, ax=ax, s=s, show_background=show_background, args=args)

    def plot_background_on_grid(self, ax=None, s=CONFIG["s_grid"], args={}):
        plot_background_on_grid(self=self, ax=ax, s=s, args=args)

    def plot_pseudotime_on_grid(self, ax=None, s=CONFIG["s_grid"], show_background=True, args={}):
        plot_pseudotime_on_grid(self=self, ax=ax, s=s, show_background=show_background, args=args)

    def plot_reference_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_reference_flow_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_simulation_flow_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_simulation_flow_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_inner_product_on_grid(self, ax=None, vm=1,s=CONFIG["s_grid"], show_background=True, args={}):
        plot_inner_product_on_grid(self=self, ax=ax, vm=vm, s=s, show_background=show_background, args=args)

    def plot_inner_product_on_pseudotime(self, ax=None, vm=1, s=CONFIG["s_grid"], args={}):
        plot_inner_product_on_pseudotime(self=self, ax=ax, vm=vm, s=s, args=args)

    def plot_inner_product_as_box(self, ax=None, vm=1, args={}):
        plot_inner_product_as_box(self=self, ax=ax, vm=vm, args=args)

    def plot_quiver(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_quiver(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args)

    def visualize_development_module_layout_2(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):

        visualize_development_module_layout_2(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)

    def visualize_development_module_layout_1(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):

        visualize_development_module_layout_1(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)

    def visualize_development_module_layout_0(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):

        visualize_development_module_layout_0(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)


def subset_oracle_for_development_analysiis(oracle_object, cell_idx_use):
    """
    Make a subset of oracle object by specifying of cluster.
    This function pick up some of attributes that needed for development analysis rather than whole attributes.

    """

    # Create new oracle object and transfer data
    oracle_ = Oracle()
    for i in ["embedding", "delta_embedding", "delta_embedding_random", "corrcoef_random", "adata"]:
        setattr(oracle_, i, getattr(oracle_object, i))

    for i in ["embedding", "delta_embedding", "delta_embedding_random", "corrcoef_random"]:
        setattr(oracle_, i, getattr(oracle_, i)[cell_idx_use])

    #cells_of_interest = oracle_object.adata.obs.index.values[cell_idx_use]
    #oracle_.adata = oracle_.adata[cells_of_interest, :]

    return oracle_

def _get_bins(array, n_bins):
    min_ = array.min()
    max_ = array.max()
    width = (max_ - min_)/(n_bins)
    return np.arange(min_, max_ + width, width)


from scipy.stats import wilcoxon

def get_p_inner_product(inner_product_df, method="wilcoxon"):
    """


    """
    dizitized_pseudotimes = np.sort(inner_product_df["pseudotime_id"].unique())

    li = []
    for i in dizitized_pseudotimes:
        pseudotime_ = inner_product_df[inner_product_df["pseudotime_id"]==i].score.values

        if method == "wilcoxon":
            _, p_ts = wilcoxon(x=pseudotime_, alternative="two-sided")
            _, p_greater = wilcoxon(x=pseudotime_, alternative="greater")
            _, p_less = wilcoxon(x=pseudotime_, alternative="less")
            mean_ = pseudotime_.mean()
            median_ = np.median(pseudotime_)
        #print(i, p)
        li.append([mean_, median_, p_ts, p_greater, p_less])

    inner_product_summary = \
        pd.DataFrame(np.array(li),
                     columns=["ip_mean", "ip_median", "p_twosided", "p_less", "p_greater"],
                     index=dizitized_pseudotimes)

    return inner_product_summary


def get_sum_of_positive_ips(inner_product_df):

    df = inner_product_df[["pseudotime_id", "score"]][inner_product_df.score > 0]
    df = df.groupby(by="pseudotime_id").sum()
    df = df.reset_index(drop=False)

    return df

def get_sum_of_negative_ips(inner_product_df):

    df = inner_product_df[["pseudotime_id", "score"]][inner_product_df.score < 0]
    df = df.groupby(by="pseudotime_id").sum()
    df = df.reset_index(drop=False)

    return df
