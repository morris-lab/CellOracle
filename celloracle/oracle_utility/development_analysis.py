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
    inner_product_stats = pd.DataFrame({"score": oracle_object.inner_product,
                                        "pseudotime": oracle_object.new_pseudotime})

    bins = _get_bins(inner_product_stats.pseudotime, n_bins)
    inner_product_stats["pseudotime_id"] = np.digitize(inner_product_stats.pseudotime, bins) - 1

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
        self.oracle_dev.mass_filter = mass_filter = oracle_object.total_p_mass < min_mass

        ## 2. Extract pseudotime data
        self.oracle_dev.pseudotime = oracle_object.adata.obs["pseudotime"].values



    def transfer_data_into_grid(self, args={}):

        if not args:
            args = {"method": "knn",
                    "n_knn": 30}

        self.oracle_dev.new_pseudotime = scatter_value_to_grid_value(embedding=self.oracle_dev.embedding,
                                                          grid=self.oracle_dev.flow_grid,
                                                          value=self.oracle_dev.pseudotime,
                                                          **args)

    def calculate_gradient_and_inner_product(self, scale_factor=60):

        # Gradient calculation
        new_pseudotime = self.oracle_dev.new_pseudotime
        n = int(np.sqrt(new_pseudotime.shape[0]))
        new_pseudotime_as_grid = new_pseudotime.reshape(n, n)
        dy, dx = np.gradient(new_pseudotime_as_grid)
        self.oracle_dev.gradient = np.stack([dx.flatten(), dy.flatten()], axis=1) * scale_factor

        # Calculate inner product between the pseudotime-gradient and the perturb-gradient
        self.oracle_dev.inner_product = np.array([np.dot(i, j) for i, j in zip(self.oracle_dev.flow, self.oracle_dev.gradient)])


    def calculate_stats(self, n_bins=10):
        inner_product_stats, inner_product_stats_grouped, = \
            get_stat_for_inner_product(self.oracle_dev, n_bins)

        self.oracle_dev.inner_product_stats = inner_product_stats
        self.oracle_dev.inner_product_stats_grouped = inner_product_stats_grouped
