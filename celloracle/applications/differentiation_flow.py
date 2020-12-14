# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys
import math

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
#import h5py
from sklearn.neighbors import KNeighborsRegressor


from ..trajectory.oracle_core import Oracle
from .scatter_to_grid import scatter_value_to_grid_value


from ..utility.hdf5_processing import dump_hdf5, load_hdf5

from ..visualizations.config import CONFIG
from ..visualizations.development_module_visualization import (\
                            plot_cluster_whole,
                            plot_cluster_cells_use,
                            plot_background,
                            plot_pseudotime,
                            plot_reference_flow_on_grid,
                            plot_pseudotime_on_grid)

def load_gradient(file_path):

    """
    Load gradient object saved as hdf5 file.

    Args:
        file_path (str): File path to the hdf5 file.
    """
    obj = load_hdf5(filename=file_path, obj_class=Gradient_calculator, ignore_attrs_if_err=[])


    return obj





class Gradient_calculator():
    def __init__(self, oracle_object=None, adata=None, obsm_key=None, pseudotime_key="Pseudotime", cell_idx_use=None, name=None, gt=None):
        """
        Estimate the direction of differentiation by calculation gradient of pseudotime on the embedding space.
        Please look at web tutorial for example scripts.

        Args:
            adata (anndata): scRNA-seq data in anndata class
            obsm_key (str): Name of dimensional reduction. You can check the list of dimensional reduction data name with "adata.obsm.keys()"
            pseudotime_key (str): Pseudotime data should be stored in adata.obs[pseudotime_key]. Please set the name of pseudotime data in adata.obs
            cluster_column_name (str): If you set cluster_column_name and cluster, you can subset cells to calculate gradient.
                Please look at web tutorial for example codes.
            cluster (str): See above.

        """
        self.cell_idx_use = None
        self.n_neighbors = None
        self.min_mass = None
        self.smooth = None
        self.n_grid = None

        if oracle_object is not None:
            self.load_oracle_object(oracle_object=oracle_object,
                                    cell_idx_use=cell_idx_use,
                                    name=name,
                                    pseudotime_key=pseudotime_key)

        elif adata is not None:
            self.load_adata(adata=adata, obsm_key=obsm_key,
                            pseudotime_key=pseudotime_key,
                            cell_idx_use=cell_idx_use,
                            name=name)

        elif gt is not None:
            self.embedding = gt.embedding.copy()
            self.mass_filter = gt.mass_filter_whole.copy()
            self.mass_filter_whole = gt.mass_filter_whole.copy()
            self.gridpoints_coordinates = gt.gridpoints_coordinates.copy()

            self.n_neighbors = gt.n_neighbors
            self.min_mass = gt.min_mass
            self.smooth = gt.smooth
            self.n_grid = gt.n_grid


    def copy(self):
        """
        Deepcopy itself.
        """
        return deepcopy(self)

    def to_hdf5(self, file_path):
        """
        Save object as hdf5.

        Args:
            file_path (str): file path to save file. Filename needs to end with '.celloracle.oracle'
        """
        if file_path.endswith(".celloracle.gradient"):
            pass
        else:
            raise ValueError("Filename needs to end with '.celloracle.gradient'")

        compression_opts = 5
        dump_hdf5(obj=self, filename=file_path,
                  data_compression=compression_opts,  chunks=(2048, 2048),
                  noarray_compression=compression_opts, pickle_protocol=4)

    def load_adata(self, adata, obsm_key, cell_idx_use=None, name=None, pseudotime_key="Pseudotime"):

        self.name = name
        self.embedding = adata.obsm[obsm_key]
        self.pseudotime = adata.obs[pseudotime_key].values

        if cell_idx_use is not None:
            self.cell_idx_use = np.array(cell_idx_use)

    def load_oracle_object(self, oracle_object, cell_idx_use=None, name=None, pseudotime_key="Pseudotime"):
        self.load_adata(adata=oracle_object.adata,
                        obsm_key=oracle_object.embedding_name,
                        cell_idx_use=cell_idx_use,
                        name=name,
                        pseudotime_key=pseudotime_key)

    def calculate_p_mass(self, smooth=0.8, n_grid=40, n_neighbors=200, n_jobs=-1):

        x_min, y_min = self.embedding.min(axis=0)
        x_max, y_max = self.embedding.max(axis=0)
        xylim = ((x_min, x_max), (y_min, y_max))
        steps = (n_grid, n_grid)

        if self.cell_idx_use is None:

            total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)
            total_p_mass_whole = total_p_mass.copy()

        else:
            total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding[self.cell_idx_use, :], smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

            total_p_mass_whole, _ = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        self.n_neighbors = n_neighbors
        self.smooth = smooth
        self.n_grid = n_grid

        self.total_p_mass = total_p_mass
        self.total_p_mass_whole = total_p_mass_whole
        self.gridpoints_coordinates = gridpoints_coordinates


    def suggest_mass_thresholds(self, n_suggestion=12, s=1, n_col=4):

        min_ = self.total_p_mass.min()
        max_ = self.total_p_mass.max()
        suggestions = np.linspace(min_, max_/2, n_suggestion)

        n_rows = math.ceil(n_suggestion / n_col)

        fig, ax = plt.subplots(n_rows, n_col, figsize=[5*n_col, 5*n_rows])
        if n_rows == 1:
            ax = ax.reshape(1, -1)

        row = 0
        col = 0
        for i in range(n_suggestion):

            ax_ = ax[row, col]

            col += 1
            if col == n_col:
                col = 0
                row += 1

            idx = self.total_p_mass > suggestions[i]

                #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax_.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=s)
            ax_.scatter(self.gridpoints_coordinates[idx, 0],
                       self.gridpoints_coordinates[idx, 1],
                       c="black", s=s)
            ax_.set_title(f"min_mass: {suggestions[i]: .2g}")
            ax_.axis("off")


    def calculate_mass_filter(self, min_mass=0.01, plot=False):

        self.min_mass = min_mass
        self.mass_filter = (self.total_p_mass < min_mass)
        self.mass_filter_whole = (self.total_p_mass_whole < min_mass)

        if plot:
            fig, ax = plt.subplots(figsize=[5,5])

            #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=10)
            ax.scatter(self.gridpoints_coordinates[~self.mass_filter, 0],
                       self.gridpoints_coordinates[~self.mass_filter, 1],
                       c="black", s=0.5)
            ax.set_title("Grid points selected")
            ax.axis("off")

    def transfer_data_into_grid(self, args={}, plot=False):

        if not args:
            args = {"method": "knn",
                    "n_knn": 30}

        if self.cell_idx_use is None:
            self.pseudotime_on_grid = scatter_value_to_grid_value(embedding=self.embedding,
                                                                  grid=self.gridpoints_coordinates,
                                                                  value=self.pseudotime,
                                                                  **args)
        else:
            self.pseudotime_on_grid = scatter_value_to_grid_value(embedding=self.embedding[self.cell_idx_use, :],
                                                                  grid=self.gridpoints_coordinates,
                                                                  value=self.pseudotime[self.cell_idx_use],
                                                                  **args)

        if plot:
            fig, ax = plt.subplots(1, 2, figsize=[10,5])

            s = 10
            s_grid = 20
            show_background = True
            ##
            ax_ = ax[0]
            plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)
            ax_.set_title("Pseudotime")


            ####
            ax_ = ax[1]
            plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background)
            ax_.set_title("Pseudotime on grid")




    def calculate_gradient(self, scale_factor="l2_norm_mean", normalization="sqrt"):

        # Gradient calculation
        gradient = get_gradient(value_on_grid=self.pseudotime_on_grid.copy())

        if normalization == "sqrt":
            gradient = normalize_gradient(gradient, method="sqrt")

        if scale_factor == "l2_norm_mean":
            # divide gradient by the mean of l2 norm.
            l2_norm = np.linalg.norm(gradient, ord=2, axis=1)
            scale_factor = 1 / l2_norm.mean()

        self.ref_flow = gradient * scale_factor


    def plot_dev_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"], show_background=True, s=CONFIG["s_scatter"], args={}):
        plot_reference_flow_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_reference_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"], show_background=True, s=CONFIG["s_scatter"], args={}):
        plot_reference_flow_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def visualize_results(self, scale=30, s=1, s_grid=30, show_background=True):


        fig, ax = plt.subplots(1, 5, figsize=[25,5])

        ##
        ax_ = ax[0]
        plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)
        ax_.set_title("Pseudotime")

        ####
        ax_ = ax[1]
        plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background)
        ax_.set_title("Pseudotime on grid")

        ###
        ax_ = ax[2]
        plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background)
        plot_reference_flow_on_grid(self, ax=ax_, scale=scale, show_background=False, s=s)
        ax_.set_title("Gradient of pseudotime \n(Development flow)")

        ###
        ax_ = ax[3]
        plot_reference_flow_on_grid(self, ax=ax_, scale=scale, show_background=show_background, s=s)
        ax_.set_title("Gradient of pseudotime \n(=Development flow)")

        ####
        ax_ = ax[4]
        plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)
        plot_reference_flow_on_grid(self, ax=ax_, scale=scale, show_background=False, s=s)
        ax_.set_title("Pseudotime + \nDevelopment flow")




    def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"],show_background=True, args=CONFIG["default_args"], cmap="rainbow"):
        plot_pseudotime(self, ax=None, s=s, show_background=show_background, args=args)


def aggregate_Gradient_objects(gradient_object_list, base_gt=None, fill_na=True):

    pseudotime_stack = [i.pseudotime_on_grid for i in gradient_object_list]
    gradient_stack = [i.ref_flow for i in gradient_object_list]
    mass_filter_stack = [i.mass_filter for i in gradient_object_list]

    new_pseudotime, new_gradient, new_mass_filter = _aggregate_gradients(pseudotime_stack=pseudotime_stack,
                                                    gradient_stack=gradient_stack,
                                                    mass_filter_stack=mass_filter_stack)

    if base_gt is None:
        gt = Gradient_calculator(gt=gradient_object_list[0])
        gt.pseudotime_on_grid = new_pseudotime
        gt.ref_flow = new_gradient

    else:
        gt = base_gt
        gt.pseudotime_on_grid[~new_mass_filter] = new_pseudotime[~new_mass_filter]
        gt.ref_flow[~new_mass_filter, :] = new_gradient[~new_mass_filter, :]

    if fill_na:
        y_with_nan= gt.pseudotime_on_grid.copy()
        y_with_nan[(new_mass_filter & (~gt.mass_filter))] = np.nan
        y_with_nan = pd.Series(y_with_nan)
        # Replace na and inf with NN's value
        y_filled =  _fill_inf_and_na(X=gt.gridpoints_coordinates[~gt.mass_filter, :],
                            y_with_inf=y_with_nan[~gt.mass_filter])

        gt.pseudotime_on_grid[~gt.mass_filter] = y_filled

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
