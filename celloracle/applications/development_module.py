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
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
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
                            plot_selected_pseudotime_on_grid,
                            plot_reference_flow_on_grid,
                            plot_simulation_flow_on_grid,
                            plot_simulation_flow_random_on_grid,
                            plot_inner_product_on_grid,
                            plot_inner_product_random_on_grid,
                            plot_inner_product_on_pseudotime,
                            plot_inner_product_as_box,
                            plot_quiver,
                            plot_quiver_random,
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
        for i in ["delta_embedding", "delta_embedding_random", "embedding", "colorandum"]:
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
        self.inner_product_random = np.array([np.dot(i, j) for i, j in zip(self.flow_rndm, self.ref_flow)])


    def calculate_digitized_ip(self, n_bins=10):

        inner_product_df = pd.DataFrame({"score": self.inner_product[~self.mass_filter_simulation],
                                         "score_randomized": self.inner_product_random[~self.mass_filter_simulation],
                                         "pseudotime": self.pseudotime_on_grid[~self.mass_filter_simulation]})


        bins = _get_bins(inner_product_df.pseudotime, n_bins)
        inner_product_df["pseudotime_id"] = np.digitize(inner_product_df.pseudotime, bins)

        self.inner_product_df = inner_product_df

    def get_negative_PS_p_value(self, pseudotime=None, return_ps_sum=False, plot=False):

        df = self.inner_product_df.copy()

        if pseudotime is not None:
            pseudotime = [i for i in pseudotime if i in list("0123456789")]
            df = df[df.pseudotime_id.isin(pseudotime)]

        x = df["score"]
        y = df["score_randomized"]

        # Clipping positive value to focus on negative IP
        x = np.clip(x, -np.inf, 0)
        y = np.clip(y, -np.inf, 0)

        if plot:
            sns.distplot(x)
            sns.distplot(y)

        # Paired non-parametric test with Wilcoxon's runk sum test
        s, p = wilcoxon(x, y, alternative="less")

        if return_ps_sum:
            return p, -x.sum(), -y.sum()

        return p

    def get_positive_PS_p_value(self, pseudotime=None, return_ps_sum=False, plot=False):
        df = self.inner_product_df.copy()

        if pseudotime is not None:
            pseudotime = [i for i in pseudotime if i in list("0123456789")]
            df = df[df.pseudotime_id.isin(pseudotime)]

        x = df["score"]
        y = df["score_randomized"]

        # Clipping negative value to focus on positive IP
        x = np.clip(x, 0, np.inf)
        y = np.clip(y, 0, np.inf)

        if plot:
            sns.distplot(x)
            sns.distplot(y)

        # Paired non-parametric test with Wilcoxon's runk sum test
        s, p = wilcoxon(x, y, alternative="greater")

        if return_ps_sum:
            return p, x.sum(), y.sum()

        return p

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

    def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"], show_background=True, cmap="rainbow", args=CONFIG["default_args"]):
        plot_pseudotime(self=self, ax=ax, s=s, show_background=show_background, cmap=cmap, args=args)

    def plot_background_on_grid(self, ax=None, s=CONFIG["s_grid"], args={}):
        plot_background_on_grid(self=self, ax=ax, s=s, args=args)

    def plot_pseudotime_on_grid(self, ax=None, s=CONFIG["s_grid"], show_background=True, cmap="rainbow", args={}):
        plot_pseudotime_on_grid(self=self, ax=ax, s=s, show_background=show_background, cmap=cmap, args=args)

    def plot_selected_pseudotime_on_grid(self, ax=None, pseudotime_selected=[], s=CONFIG["s_grid"], show_background=True, args={}):
        plot_selected_pseudotime_on_grid(self=self, ax=ax, pseudotime_selected=pseudotime_selected, s=s, show_background=show_background, args=args)

    def plot_reference_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_reference_flow_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_simulation_flow_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_simulation_flow_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_simulation_flow_random_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_simulation_flow_random_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_inner_product_on_grid(self, ax=None, vm=1, s=CONFIG["s_grid"], show_background=True, vmin=None, vmax=None, cmap=None, args={}):
        plot_inner_product_on_grid(self=self, ax=ax, vm=vm, s=s, show_background=show_background, vmin=vmin, vmax=vmax, cmap=cmap, args=args)

    def plot_inner_product_random_on_grid(self, ax=None, vm=1, s=CONFIG["s_grid"], show_background=True, vmin=None, vmax=None, cmap=None, args={}):
        plot_inner_product_random_on_grid(self=self, ax=ax, vm=vm, s=s, show_background=show_background, vmin=vmin, vmax=vmax, cmap=cmap, args=args)

    def plot_inner_product_on_pseudotime(self, ax=None, vm=1, s=CONFIG["s_grid"], vmin=None, vmax=None, cmap=None, args={}):
        plot_inner_product_on_pseudotime(self=self, ax=ax, vm=vm, s=s, vmin=vmin, vmax=vmax, cmap=cmap, args=args)

    def plot_inner_product_as_box(self, ax=None, vm=1, vmin=None, vmax=None, args={}):
        plot_inner_product_as_box(self=self, ax=ax, vm=vm, vmin=vmin, vmax=vmax, args=args)

    def plot_quiver(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_quiver(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args)

    def plot_quiver_random(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_quiver_random(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args)

    def visualize_development_module_layout_2(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True, return_fig=False):

        fig = visualize_development_module_layout_2(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)

        if return_fig:
            return fig

    def visualize_development_module_layout_1(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True, return_fig=False):

        fig = visualize_development_module_layout_1(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)

        if return_fig:
            return fig

    def visualize_development_module_layout_0(self, scale_for_pseudotime=CONFIG["scale_dev"],
        scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True, return_fig=False):

        fig = visualize_development_module_layout_0(self=self, scale_for_pseudotime=scale_for_pseudotime,
            scale_for_simulation=scale_for_simulation, s=s, s_grid=s_grid, vm=vm, show_background=show_background)

        if return_fig:
            return fig


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

"""def _get_bins(array, n_bins):
    min_ = array.min()
    max_ = array.max()
    width = (max_ - min_)/(n_bins-1)
    return np.arange(min_, max_ + width, width)"""

def _get_bins(array, n_bins):
    min_ = array.min()
    max_ = array.max()
    width = (max_ - min_)/(n_bins)
    return np.arange(min_, max_ + width, width)[1:-1]


from scipy.stats import wilcoxon



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
