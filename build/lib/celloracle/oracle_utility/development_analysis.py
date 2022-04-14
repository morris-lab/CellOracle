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
#import h5py

from ..trajectory.oracle_core import Oracle
from .scatter_to_grid import scatter_value_to_grid_value


def subset_oracle_for_development_analysiis(oracle_object, cluster_column_name, cluster):
    """
    Make a subset of oracle object by specifying of cluster.
    This function pick up some of attributes that needed for development analysis rather than whole attributes.

    """

    cells_of_interest = oracle_object.adata.obs[oracle_object.adata.obs[cluster_column_name] == cluster].index.values

    # Create new oracle object and transfer data
    oracle_ = Oracle()
    for i in ["embedding", "delta_embedding", "delta_embedding_random", "corrcoef_random", "adata"]:
        setattr(oracle_, i, getattr(oracle_object, i))

    index_use = np.where(oracle_.adata.obs.index.isin(cells_of_interest))[0]
    for i in ["embedding", "delta_embedding", "delta_embedding_random", "corrcoef_random"]:
        setattr(oracle_, i, getattr(oracle_, i)[index_use])
    oracle_.adata = oracle_.adata[cells_of_interest, :]

    return oracle_

from scipy.stats import wilcoxon

def get_stat_for_inner_product(oracle_object, n_bins=10):

    # Prepare data
    inner_product_stats = pd.DataFrame({"score": oracle_object.inner_product[~oracle_object.mass_filter],
                                        "pseudotime": oracle_object.new_pseudotime[~oracle_object.mass_filter]})


    bins = _get_bins(inner_product_stats.pseudotime, n_bins)
    inner_product_stats["pseudotime_id"] = np.digitize(inner_product_stats.pseudotime, bins) - 1

    try:
        inner_product_stats["stage"] = oracle_object.stage_grid[~oracle_object.mass_filter]

    except:
        print("stage_grid not calculated")


    # stat test
    ps = []
    for i in np.sort(inner_product_stats["pseudotime_id"].unique()):
        pseudotime_ = inner_product_stats[inner_product_stats["pseudotime_id"]==i].score.values
        stat, p = wilcoxon(x=pseudotime_, alternative="less")
        #print(i, p)
        ps.append(p)

    inner_product_stats_grouped = \
        pd.DataFrame({"ip_score_median":  inner_product_stats.groupby("pseudotime_id").median()["score"].values,
                      "ip_score_mean":  inner_product_stats.groupby("pseudotime_id").mean()["score"],
                      "p-value_negative":  ps},
                      index=np.sort(inner_product_stats["pseudotime_id"].unique()))

    return inner_product_stats, inner_product_stats_grouped


def _get_bins(array, n_bins):
    min_ = array.min()
    max_ = array.max()
    width = (max_ - min_)/(n_bins-1)
    return np.arange(min_, max_ + width, width)


class Oracle_development_module():
    def __init__(self):
        pass

    def extract_data_from_oracle(self, oracle_object, min_mass):

        self.oracle_dev = Oracle()
        ## 1. Extract perturb simulation results as grid matrix
        # 1.1.grid matrix and embedding
        for i in ["flow_grid", "flow", "flow_norm_rndm", "embedding"]:
            setattr(self.oracle_dev, i, getattr(oracle_object, i))

        # 1.2. mass_filter for grid matrix
        self.oracle_dev.mass_filter = (oracle_object.total_p_mass < min_mass)

        ## 2. Extract pseudotime data
        self.oracle_dev.pseudotime = oracle_object.adata.obs["pseudotime"].values

        try:
            self.oracle_dev.stage = np.array(oracle_object.adata.obs["Stage"].values)
        except:
            print("Stage not in data")




    def transfer_data_into_grid(self, args={}):

        if not args:
            args = {"method": "knn",
                    "n_knn": 30}

        self.oracle_dev.new_pseudotime = scatter_value_to_grid_value(embedding=self.oracle_dev.embedding,
                                                          grid=self.oracle_dev.flow_grid,
                                                          value=self.oracle_dev.pseudotime,
                                                          **args)
        try:
            self.oracle_dev.stage_grid = scatter_value_to_grid_value(embedding=self.oracle_dev.embedding,
                                                     grid=self.oracle_dev.flow_grid,
                                                     value=self.oracle_dev.stage,
                                                     **{"method": "knn_class",
                                                     "n_knn": 30})
        except:
            print("Stage not in data")

    def calculate_gradient_and_inner_product(self, scale_factor="l2_norm_mean", normalization=None):

        # Gradient calculation
        gradient = get_gradient(value_on_grid=self.oracle_dev.new_pseudotime.copy())

        if normalization == "sqrt":
            gradient = normalize_gradient(gradient, method="sqrt")

        if scale_factor == "l2_norm_mean":
            # divide gradient by the mean of l2 norm.
            l2_norm = np.linalg.norm(gradient, ord=2, axis=1)
            scale_factor = 1 / l2_norm.mean()

        self.oracle_dev.gradient = gradient * scale_factor

        # Calculate inner product between the pseudotime-gradient and the perturb-gradient
        self.oracle_dev.inner_product = np.array([np.dot(i, j) for i, j in zip(self.oracle_dev.flow, self.oracle_dev.gradient)])


    def calculate_stats(self, n_bins=10):
        inner_product_stats, inner_product_stats_grouped, = \
            get_stat_for_inner_product(self.oracle_dev, n_bins)

        self.oracle_dev.inner_product_stats = inner_product_stats
        self.oracle_dev.inner_product_stats_grouped = inner_product_stats_grouped



class Gradient_based_trajecory():
    def __init__(self, adata=None, obsm_key=None, pseudotime_key="pseudotime", cluster_column_name=None, cluster=None, gt=None):

        if adata is not None:
            self.load_adata(adata=adata, obsm_key=obsm_key,
            pseudotime_key=pseudotime_key,cluster_column_name=cluster_column_name,
            cluster=cluster)
        elif gt is not None:
            self.embedding = gt.embedding_whole.copy()
            self.embedding_whole = gt.embedding_whole.copy()
            self.mass_filter = gt.mass_filter_whole.copy()
            self.mass_filter_whole = gt.mass_filter_whole.copy()
            self.gridpoints_coordinates = gt.gridpoints_coordinates.copy()
            self.pseudotime = gt.pseudotime_whole.copy()

    def load_adata(self, adata, obsm_key, pseudotime_key, cluster_column_name=None, cluster=None):

        self.embedding = adata.obsm[obsm_key]
        self.pseudotime = adata.obs[pseudotime_key].values
        self.embedding_whole = self.embedding.copy()
        self.pseudotime_whole = self.pseudotime.copy()

        if (cluster_column_name is not None) & (cluster is not None):
            cells_ix = np.where(adata.obs[cluster_column_name] == cluster)[0]
            self.embedding = self.embedding[cells_ix, :]
            self.pseudotime = self.pseudotime[cells_ix]





    def calculate_mass_filter(self, min_mass=0.01, smooth=0.8, steps=(40, 40), n_neighbors=200, n_jobs=4):

        x_min, y_min = self.embedding_whole.min(axis=0)
        x_max, y_max = self.embedding_whole.max(axis=0)
        xylim = ((x_min, x_max), (y_min, y_max))

        total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        total_p_mass_whole, _ = calculate_p_mass(self.embedding_whole, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        self.total_p_mass = total_p_mass
        self.mass_filter = (total_p_mass < min_mass)
        self.mass_filter_whole = (total_p_mass_whole < min_mass)
        self.gridpoints_coordinates = gridpoints_coordinates

    def transfer_data_into_grid(self, args={}):

        if not args:
            args = {"method": "knn",
                    "n_knn": 30}

        self.pseudotime_on_grid = scatter_value_to_grid_value(embedding=self.embedding,
                                                          grid=self.gridpoints_coordinates,
                                                          value=self.pseudotime,
                                                          **args)
    def calculate_gradient(self, scale_factor=60, normalization=None):

        # Gradient calculation
        gradient = get_gradient(value_on_grid=self.pseudotime_on_grid.copy())

        if normalization == "sqrt":
            gradient = normalize_gradient(gradient, method="sqrt")

        if scale_factor == "l2_norm_mean":
            # divide gradient by the mean of l2 norm.
            l2_norm = np.linalg.norm(gradient, ord=2, axis=1)
            scale_factor = 1 / l2_norm.mean()

        self.gradient = gradient * scale_factor

    def visualize_dev_flow(self, scale_for_pseudotime=30, s=10, s_grid=30):
        visualize_dev_flow(self, scale_for_pseudotime=scale_for_pseudotime, s=s, s_grid=s_grid)

def aggregate_GT_object(list_GT_object, base_gt=None):

    pseudotime_stack = [i.pseudotime_on_grid for i in list_GT_object]
    gradient_stack = [i.gradient for i in list_GT_object]
    mass_filter_stack = [i.mass_filter for i in list_GT_object]

    new_pseudotime, new_gradient, new_mass_filter = _aggregate_gradients(pseudotime_stack=pseudotime_stack,
                                                    gradient_stack=gradient_stack,
                                                    mass_filter_stack=mass_filter_stack)

    if base_gt is None:
        gt = Gradient_based_trajecory(gt=list_GT_object[0])
        gt.pseudotime_on_grid = new_pseudotime
        gt.gradient = new_gradient

    else:
        gt = base_gt
        gt.pseudotime_on_grid[~new_mass_filter] = new_pseudotime[~new_mass_filter]
        gt.gradient[~new_mass_filter, :] = new_gradient[~new_mass_filter, :]

    return gt

def _aggregate_gradients(pseudotime_stack, gradient_stack, mass_filter_stack):

    new_pseudotime = np.zeros_like(pseudotime_stack[0])
    new_pseudotime_count = np.zeros_like(pseudotime_stack[0])
    new_gradient = np.zeros_like(gradient_stack[0])
    gradient_count = np.zeros_like(gradient_stack[0])
    for fil, pt, gra in zip(mass_filter_stack, pseudotime_stack, gradient_stack):
        new_pseudotime[~fil] += pt[~fil]
        new_pseudotime_count[~fil] +=1
        new_gradient[~fil, :] += gra[~fil, :]
        gradient_count[~fil, :] += 1

    new_pseudotime[new_pseudotime_count != 0] /= new_pseudotime_count[new_pseudotime_count != 0]
    new_gradient[gradient_count != 0] /= gradient_count[gradient_count != 0]
    new_mass_filter = (gradient_count.sum(axis=1) == 0)

    return new_pseudotime, new_gradient, new_mass_filter


def normalize_gradient(gradient, method="sqrt"):
    """
    Normalize length of 2D vector
    """

    if method == "sqrt":

        size = np.sqrt(np.power(gradient, 2).sum(axis=1))
        size_sq = np.sqrt(size)
        size_sq[size_sq == 0] = 1
        factor = np.repeat(np.expand_dims(size_sq, axis=1), 2, axis=1)

    return gradient / factor

from scipy.stats import norm as normal
from sklearn.neighbors import NearestNeighbors


def calculate_p_mass(embedding, smooth=0.5, steps=(40, 40),
                          n_neighbors=100, n_jobs=4, xylim=((None, None), (None, None))):
    """Calculate the velocity using a points on a regular grid and a gaussian kernel

    Note: the function should work also for n-dimensional grid

    Arguments
    ---------
    embedding:

    smooth: float, smooth=0.5
        Higher value correspond to taking in consideration further points
        the standard deviation of the gaussian kernel is smooth * stepsize
    steps: tuple, default
        the number of steps in the grid for each axis
    n_neighbors:
        number of neighbors to use in the calculation, bigger number should not change too much the results..
        ...as soon as smooth is small
        Higher value correspond to slower execution time
    n_jobs:
        number of processes for parallel computing
    xymin:
        ((xmin, xmax), (ymin, ymax))

    Returns
    -------
    total_p_mass: np.ndarray
        density at each point of the grid

    """

    # Prepare the grid
    grs = []
    for dim_i in range(embedding.shape[1]):
        m, M = np.min(embedding[:, dim_i]), np.max(embedding[:, dim_i])

        if xylim[dim_i][0] is not None:
            m = xylim[dim_i][0]
        if xylim[dim_i][1] is not None:
            M = xylim[dim_i][1]

        m = m - 0.025 * np.abs(M - m)
        M = M + 0.025 * np.abs(M - m)
        gr = np.linspace(m, M, steps[dim_i])
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
    nn.fit(embedding)
    dists, neighs = nn.kneighbors(gridpoints_coordinates)

    std = np.mean([(g[1] - g[0]) for g in grs])
    # isotropic gaussian kernel
    gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
    total_p_mass = gaussian_w.sum(1)
    gridpoints_coordinates

    return total_p_mass, gridpoints_coordinates

def get_gradient(value_on_grid):
    # Gradient calculation
    n = int(np.sqrt(value_on_grid.shape[0]))
    value_on_grid_as_matrix = value_on_grid.reshape(n, n)
    dy, dx = np.gradient(value_on_grid_as_matrix)
    gradient = np.stack([dx.flatten(), dy.flatten()], axis=1)

    return gradient


def visualize_dev_flow(self, scale_for_pseudotime=30, s=10, s_grid=30):

    embedding_whole = self.embedding_whole
    embedding_of_interest= self.embedding
    mass_filter = self.mass_filter
    mass_filter_whole = self.mass_filter_whole
    gridpoints_coordinates=self.gridpoints_coordinates

    pseudotime_raw = self.pseudotime
    pseudotime_on_grid=self.pseudotime_on_grid

    gradient_pseudotime=self.gradient


    fig, ax = plt.subplots(1, 5, figsize=[25,5])

    ##
    ax_ = ax[0]
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.scatter(embedding_of_interest[:, 0], embedding_of_interest[:, 1], c=pseudotime_raw, cmap="rainbow", s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####
    ax_ = ax[1]
    ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[~mass_filter_whole, 0], gridpoints_coordinates[~mass_filter_whole, 1],
     c="lightgray", s=s_grid)
    ax_.scatter(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
     c=pseudotime_on_grid[~mass_filter], cmap="rainbow", s=s_grid)
    ax_.set_title("Pseudotime on grid")
    ax_.axis("off")



    ###
    ax_ = ax[2]
    #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[~mass_filter_whole, 0], gridpoints_coordinates[~mass_filter_whole, 1],
     c="lightgray", s=s_grid)
    ax_.scatter(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
     c=pseudotime_on_grid[~mass_filter], cmap="rainbow", s=s_grid)

    ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    ###
    ax_ = ax[3]
    #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    ####
    ax_ = ax[4]
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.scatter(embedding_of_interest[:, 0], embedding_of_interest[:, 1], c=pseudotime_raw, cmap="rainbow", s=s)

    ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Pseudotime + \nDevelopment flow")
    ax_.axis("off")
