# -*- coding: utf-8 -*-


import os, sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from ..trajectory.oracle_utility import _adata_to_color_dict


def visualize_developmental_analysis_ver1(self, scale_for_pseudotime=30, scale_for_simulated=30, s=10, s_grid=30, vmin=-1, vmax=1):

    original_embedding=self.oracle_dev.embedding
    original_value=self.oracle_dev.pseudotime
    mass_filter=self.oracle_dev.mass_filter
    grid=self.oracle_dev.flow_grid
    value_on_grid=self.oracle_dev.new_pseudotime
    gradient_pseudotime=self.oracle_dev.gradient
    gradient_simulated=self.oracle_dev.flow
    inner_product=self.oracle_dev.inner_product


    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    ax[0, 0].scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax[0, 0].set_title("Pseudotime")
    ax[0,0].axis("off")

    ax[0, 1].scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=value_on_grid[~mass_filter], cmap="rainbow", s=s_grid)
    ax[0, 1].set_title("Pseudotime on grid")
    ax[0, 1].axis("off")

    ax[0, 2].quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax[0, 2].set_title("Gradient of pseudotime \n(=Development flow)")
    ax[0, 2].axis("off")



    ax[0, 3].scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax[0, 3].quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax[0, 3].set_title("Pseudotime + \nDevelopment flow")
    ax[0, 3].axis("off")

    ax[1, 0].quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax[1, 0].set_title("Perturb simulation")
    ax[1, 0].axis("off")


    ax[1, 1].scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter], cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)
    ax[1, 1].set_title("Inner product of \n Perturb simulation * Development flow")
    ax[1, 1].axis("off")


    ax[1, 2].scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter], cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)
    ax[1, 2].quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax[1, 2].axis("off")
    ax[1, 2].set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")

    ax_ = ax[1, 3]
    pcm = ax_.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)
    ax_.set_ylim([vmin*1.1, vmax*1.1])

    ax[1, 3].axhline(0, color="lightgray")
    pp = fig.colorbar(pcm, ax=ax_, orientation="vertical")
    sns.despine()

    ax[1, 3].set_xlabel("pseudotime")
    ax[1, 3].set_ylabel("inner product score")




def visualize_developmental_analysis_ver2(self, scale_for_pseudotime=30, scale_for_simulated=30, s=10, s_grid=30, vmin=-1, vmax=1):

    whole_embedding = self.oracle.embedding
    original_embedding=self.oracle_dev.embedding
    original_value=self.oracle_dev.pseudotime
    mass_filter=self.oracle_dev.mass_filter
    grid=self.oracle_dev.flow_grid
    value_on_grid=self.oracle_dev.new_pseudotime
    gradient_pseudotime=self.oracle_dev.gradient
    gradient_simulated=self.oracle_dev.flow
    inner_product=self.oracle_dev.inner_product

    cluster = self.oracle_dev.cluster_loaded
    cluster_column_name = self.oracle_dev.cluster_column_name_loaded

    cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]


    inner_product_stats = self.oracle_dev.inner_product_stats
    inner_product_stats_grouped = self.oracle_dev.inner_product_stats_grouped

    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    ax[0, 0].scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax[0, 0].scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
    ax[0, 0].set_title(f"Cluster of interest: {cluster}")
    ax[0, 0].axis("off")

    ax[0, 1].scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax[0, 1].set_title("Pseudotime")
    ax[0, 1].axis("off")



    ax_ = ax[0, 2]
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    ax_ = ax[0, 3]
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation")
    ax_.axis("off")


    ax[1, 0].scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter], cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)
    ax[1, 0].set_title("Inner product of \n Perturb simulation * Development flow")
    ax[1, 0].axis("off")


    ax[1, 1].scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter], cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)
    ax[1, 1].quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax[1, 1].axis("off")
    ax[1, 1].set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")

    ax_ = ax[1, 2]
    pcm = ax_.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)
    ax_.set_ylim([vmin*1.1, vmax*1.1])

    ax_.axhline(0, color="lightgray")
    pp = fig.colorbar(pcm, ax=ax_, orientation="vertical")
    sns.despine()

    ax_.set_xlabel("pseudotime")
    ax_.set_ylabel("inner product score")


    ax_ = ax[1, 3]
    sns.boxplot(data=inner_product_stats, x="pseudotime_id", y="score", color="white", ax=ax_)
    ax_.set_xlabel("Digitized_pseudotime")
    ax_.set_ylabel("inner product score")
    ax_.axhline(0, color="gray")
    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.tick_params(
                labelleft=False)
