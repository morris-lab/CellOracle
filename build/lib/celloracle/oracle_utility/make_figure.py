# -*- coding: utf-8 -*-


import os, sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from ..trajectory.oracle_utility import _adata_to_color_dict


def _get_ix_for_a_cluster(oracle, cluster_column_name, cluster):
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix

def _plot_quiver_for_a_cluster(oracle, cluster_column_name, cluster, quiver_scale, color=None, plot_whole_cells=True, args={}):

    if cluster == "whole":
        ix_choice = ix = np.arange(oracle.adata.shape[0])
    else:
        ix_choice = _get_ix_for_a_cluster(oracle, cluster_column_name, cluster)


    if plot_whole_cells:

        plt.scatter(oracle.embedding[:, 0], oracle.embedding[:, 1],
                    c="lightgray", alpha=1, lw=0.3, rasterized=True,  **args)

    plt.scatter(oracle.embedding[ix_choice, 0], oracle.embedding[ix_choice, 1],
                c="lightgray", alpha=0.2, edgecolor=(0,0,0,1), lw=0.3, rasterized=True, **args)



    if color is None:
        color=oracle.colorandum[ix_choice]

    quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                       linewidths=0.25, width=0.0045,edgecolors="k",
                       color=color, alpha=1)

    plt.quiver(oracle.embedding[ix_choice, 0], oracle.embedding[ix_choice, 1],
               oracle.delta_embedding[ix_choice, 0],
               oracle.delta_embedding[ix_choice, 1],
               scale=quiver_scale, **quiver_kwargs)

    plt.axis("off")

def plot_scatter_with_anndata(adata, obsm_key, cluster_column_name, args={}):

    embedding = adata.obsm[obsm_key]
    colors = _adata_to_color_dict(adata=adata, cluster_use=cluster_column_name)

    for cluster, color in colors.items():
        idx = np.where(adata.obs[cluster_column_name] == cluster)[0]
        plt.scatter(embedding[idx, 0], embedding[idx, 1], c=color, label=cluster, **args)





def figures_for_trajectories301(self, save_folder, scale_for_pseudotime=30, scale_for_simulated=30, quiver_scale=30, s=10, s_grid=30, vmin=-1, vmax=1, figsize=[5, 5], fontsize=15):

    whole_embedding = self.oracle.embedding
    original_embedding=self.oracle_dev.embedding
    original_value=self.oracle_dev.pseudotime
    mass_filter=self.oracle_dev.mass_filter
    grid=self.oracle_dev.flow_grid
    value_on_grid=self.oracle_dev.new_pseudotime
    gradient_pseudotime=self.oracle_dev.gradient
    gradient_simulated=self.oracle_dev.flow
    inner_product=self.oracle_dev.inner_product

    inner_product_stats = self.oracle_dev.inner_product_stats
    inner_product_stats_grouped = self.oracle_dev.inner_product_stats_grouped

    alpha = 1

    cluster = self.oracle_dev.cluster_loaded
    #if cluster == "True":
    #    cluster = True

    cluster_column_name = self.oracle_dev.cluster_column_name_loaded



    fig = plt.figure(figsize=figsize)
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key="X_umap",
                                  cluster_column_name=cluster_column_name,
                                  args={"s": s})

    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]

        plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"scatter_{cluster_column_name}_{cluster}.png"), transparent=True)


    ##
    fig = plt.figure(figsize=figsize)
    plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    #plt.title("Pseudotime")
    plt.axis("off")

    plt.savefig(os.path.join(save_folder, f"pseudotime_{cluster_column_name}_{cluster}.png"), transparent=True)


    ####
    ###
    fig = plt.figure(figsize=figsize)
    plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    #plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    plt.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    #plt.title("Differentiation")
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"differentiation_{cluster_column_name}_{cluster}.png"), transparent=True)



    fig = plt.figure(figsize=figsize)
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key="X_umap",
                                  cluster_column_name=cluster_column_name,
                                  args={"s": s})
    else:
        plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    #plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    plt.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    #plt.title("Gradient of pseudotime \n(=Development flow)")
    plt.axis("off")

    plt.savefig(os.path.join(save_folder, f"differentiation_with_cells_of_interest_{cluster_column_name}_{cluster}.png"), transparent=True)



def figures_for_perturb_analysis_301(self, save_folder, scale_for_pseudotime=30, scale_for_simulated=30, quiver_scale=30, s=10, s_grid=30, vmin=-1, vmax=1, figsize=[5, 5], fontsize=15):

    whole_embedding = self.oracle.embedding
    original_embedding=self.oracle_dev.embedding
    original_value=self.oracle_dev.pseudotime
    mass_filter=self.oracle_dev.mass_filter
    grid=self.oracle_dev.flow_grid
    value_on_grid=self.oracle_dev.new_pseudotime
    gradient_pseudotime=self.oracle_dev.gradient
    gradient_simulated=self.oracle_dev.flow
    inner_product=self.oracle_dev.inner_product

    inner_product_stats = self.oracle_dev.inner_product_stats
    inner_product_stats_grouped = self.oracle_dev.inner_product_stats_grouped

    alpha = 1

    cluster = self.oracle_dev.cluster_loaded
    #if cluster == "True":
    #    cluster = True

    cluster_column_name = self.oracle_dev.cluster_column_name_loaded

    if cluster == "whole":
        pass
    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]



    #####
    fig = plt.figure(figsize=figsize)
    #plt.title(f"Perturb simulation \n color: {cluster_column_name}")
    if cluster == "whole":
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=None,
                                   cluster=cluster, quiver_scale=30, args={"s": s})
    else:
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=cluster_color,
                                   cluster=cluster, quiver_scale=30, args={"s": s})
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"quiver_full_on_cellls_of_interest_{cluster_column_name}_{cluster}.png"), transparent=True)


    #######
    fig = plt.figure(figsize=figsize)
    #plt.title("Perturb simulation \n color: cluster")
    _plot_quiver_for_a_cluster(oracle=self.oracle,
                               cluster_column_name=cluster_column_name,
                               cluster=cluster, quiver_scale=30, args={"s": s})
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"quiver_ful_on_cluster_{cluster_column_name}_{cluster}.png"), transparent=True)



    ########
    fig = plt.figure(figsize=figsize)
    plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    plt.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    #plt.title("Perturb simulation result on grid")
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"quiver_grid_{cluster_column_name}_{cluster}.png"), transparent=True)


    ##########
    fig = plt.figure(figsize=figsize)
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key="X_umap",
                                  cluster_column_name=cluster_column_name,
                                  args={"s": s})
    else:
        plt.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        plt.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    plt.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    #plt.title("Perturb simulation result on grid")
    plt.axis("off")
    plt.savefig(os.path.join(save_folder, f"quiver_grid_on_cells_of_interest_{cluster_column_name}_{cluster}.png"), transparent=True)



    #########

    fig = plt.figure(figsize=figsize)
    plt.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    plt.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    plt.axis("off")
    #plt.title("Inner product of \n Perturb simulation * Development flow")
    plt.savefig(os.path.join(save_folder, f"inner_product_score_{cluster_column_name}_{cluster}.png"), transparent=True)



    fig = plt.figure(figsize=figsize)
    plt.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    plt.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    plt.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    plt.axis("off")
    #plt.title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")
    plt.savefig(os.path.join(save_folder, f"inner_product_score_and_perturb_quiver_grid_{cluster_column_name}_{cluster}.png"), transparent=True)


    #####
    #fig = plt.figure(figsize=figsize)
    fig, ax = plt.subplots(figsize=figsize)
    pcm = ax.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)

    plt.ylim([vmin*1.1, vmax*1.1])
    plt.axhline(0, color="lightgray")
    pp = fig.colorbar(pcm, ax=ax, orientation="vertical")
    sns.despine()
    plt.xlabel("pseudotime")
    plt.ylabel("inner product score")
    plt.savefig(os.path.join(save_folder, f"inner_product_score_distribution_{cluster_column_name}_{cluster}.png"), transparent=True)


    fig = plt.figure(figsize=figsize)
    #fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(data=inner_product_stats, x="pseudotime_id", y="score", color="white")
    plt.xlabel("Digitized_pseudotime")
    plt.ylabel("inner product score")
    plt.axhline(0, color="gray")
    plt.ylim([vmin*1.1, vmax*1.1])
    plt.tick_params(
                labelleft=False)
    plt.show()
    plt.savefig(os.path.join(save_folder, f"inner_product_score_distribution_box_plot_{cluster_column_name}_{cluster}.png"), transparent=True)
