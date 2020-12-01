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


def plot_scatter_with_anndata(adata, obsm_key, cluster_column_name, ax, args={}):

    embedding = adata.obsm[obsm_key]
    colors = _adata_to_color_dict(adata=adata, cluster_use=cluster_column_name)

    for cluster, color in colors.items():
        idx = np.where(adata.obs[cluster_column_name] == cluster)[0]
        ax.scatter(embedding[idx, 0], embedding[idx, 1], c=color, label=cluster, **args)



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

    inner_product_stats = self.oracle_dev.inner_product_stats
    inner_product_stats_grouped = self.oracle_dev.inner_product_stats_grouped


    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax[0, 0], args={"s": s})
        ax[0, 0].set_title(f"Cluster of interest: all clusters")
        ax[0, 0].axis("off")
    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]

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




def visualize_developmental_analysis_ver101(self, scale_for_pseudotime=30, scale_for_simulated=30, s=10, s_grid=30, vmin=-1, vmax=1):

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
    #if cluster == "True":
    #    cluster = True

    cluster_column_name = self.oracle_dev.cluster_column_name_loaded


    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    ax_ = ax[0, 0]

    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax[0, 0], args={"s": s})
        ax_.set_title(f"Cluster of interest: all clusters")

    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]


        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
        if cluster == "True" :
            ax_.set_title(f"Cells of interest: \n{cluster_column_name}")

    ax_.axis("off")

    ##
    ax_ = ax[0, 1]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####
    ax_ = ax[0, 2]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=value_on_grid[~mass_filter], cmap="rainbow", s=s_grid)
    ax_.set_title("Pseudotime on grid")
    ax_.axis("off")


    ###
    ax_ = ax[0, 3]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")


    ####
    ax_ = ax[1, 0]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Pseudotime + \nDevelopment flow")
    ax_.axis("off")

    ####
    ax_ = ax[1, 1]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation")
    ax_.axis("off")

    ax_ = ax[1, 2]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")

    #####
    ax_ = ax[1, 3]
    pcm = ax_.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)

    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.axhline(0, color="lightgray")
    pp = fig.colorbar(pcm, ax=ax_, orientation="vertical")
    sns.despine()
    ax_.set_xlabel("pseudotime")
    ax_.set_ylabel("inner product score")



def visualize_developmental_analysis_ver201(self, scale_for_pseudotime=30, scale_for_simulated=30, s=10, s_grid=30, vmin=-1, vmax=1):

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


    cluster = self.oracle_dev.cluster_loaded
    #if cluster == "True":
    #    cluster = True

    cluster_column_name = self.oracle_dev.cluster_column_name_loaded


    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    ax_ = ax[0, 0]

    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
        ax_.set_title(f"Cluster of interest: all clusters")

    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]

        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
        if cluster == True :
            ax_.set_title(f"Cells of interest: \n{cluster_column_name}")

    ax_.axis("off")

    ##
    ax_ = ax[0, 1]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####

    ###
    ax_ = ax[0, 2]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    #####

    ax_ = ax[0, 3]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation")
    ax_.axis("off")


    ####
    ax_ = ax[1, 0]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[1, 1]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")


    #####
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



#####
def _get_ix_for_a_cluster(oracle, cluster_column_name, cluster):
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix

def _plot_quiver_for_a_cluster(oracle, cluster_column_name, cluster, quiver_scale, ax, color=None, plot_whole_cells=True, args={}):

    if cluster == "whole":
        ix_choice = ix = np.arange(oracle.adata.shape[0])
    else:
        ix_choice = _get_ix_for_a_cluster(oracle, cluster_column_name, cluster)


    if plot_whole_cells:

        ax.scatter(oracle.embedding[:, 0], oracle.embedding[:, 1],
                    c="lightgray", alpha=1, lw=0.3, rasterized=True,  **args)

    ax.scatter(oracle.embedding[ix_choice, 0], oracle.embedding[ix_choice, 1],
                c="lightgray", alpha=0.2, edgecolor=(0,0,0,1), lw=0.3, rasterized=True, **args)



    if color is None:
        color=oracle.colorandum[ix_choice]

    quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                       linewidths=0.25, width=0.0045,edgecolors="k",
                       color=color, alpha=1)

    ax.quiver(oracle.embedding[ix_choice, 0], oracle.embedding[ix_choice, 1],
               oracle.delta_embedding[ix_choice, 0],
               oracle.delta_embedding[ix_choice, 1],
               scale=quiver_scale, **quiver_kwargs)

    plt.axis("off")




def visualize_developmental_analysis_ver301(self, scale_for_pseudotime=30, scale_for_simulated=30, quiver_scale=30, s=10, s_grid=30, vmin=-1, vmax=1):

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


    ''' fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ax_ = ax[0]
    ax_.set_title(f"Clustering results")
    plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,,
                                  cluster_column_name=self.oracle.cluster_column_name,
                                  ax=ax_, args={"s": s})
    ax_.axis("off")'''


    fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ax_ = ax[0]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
        ax_.set_title(f"Cluster of interest: all clusters")

    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]

        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
        if cluster == "True" :
            ax_.set_title(f"Cells of interest: \n{cluster_column_name}")
    ax_.axis("off")

    ##
    ax_ = ax[1]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####
    ###
    ax_ = ax[2]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    #ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")


    ax_ = ax[3]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
    else:
        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    #ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    plt.show()



    fig, ax = plt.subplots(1, 4, figsize=[20,5])


    #####
    ax_ = ax[0]
    ax_.set_title(f"Perturb simulation \n color: {cluster_column_name}")
    if cluster == "whole":
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=None,
                                   cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    else:
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=cluster_color,
                                   cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    ax_.axis("off")

    ax_ = ax[1]
    ax_.set_title("Perturb simulation \n color: cluster")
    _plot_quiver_for_a_cluster(oracle=self.oracle,
                               cluster_column_name=cluster_column_name,
                               cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    ax_.axis("off")

    ax_ = ax[2]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation result on grid")
    ax_.axis("off")


    ax_ = ax[3]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
    else:
        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation result on grid")
    ax_.axis("off")

    plt.show()

    fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ####
    ax_ = ax[0]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[1]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")


    #####
    ax_ = ax[2]
    pcm = ax_.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)

    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.axhline(0, color="lightgray")
    pp = fig.colorbar(pcm, ax=ax_, orientation="vertical")
    sns.despine()
    ax_.set_xlabel("pseudotime")
    ax_.set_ylabel("inner product score")

    ax_ = ax[3]
    sns.boxplot(data=inner_product_stats, x="pseudotime_id", y="score", color="white", ax=ax_)
    ax_.set_xlabel("Digitized_pseudotime")
    ax_.set_ylabel("inner product score")
    ax_.axhline(0, color="gray")
    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.tick_params(
                labelleft=False)
    plt.show()



def visualize_developmental_analysis_ver401(self, scale_for_pseudotime=30, scale_for_simulated=30, quiver_scale=30, s=10, s_grid=30, vmin=-1, vmax=1):

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


    ''' fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ax_ = ax[0]
    ax_.set_title(f"Clustering results")
    plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,,
                                  cluster_column_name=self.oracle.cluster_column_name,
                                  ax=ax_, args={"s": s})
    ax_.axis("off")'''


    fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ax_ = ax[0]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
        ax_.set_title(f"Cluster of interest: all clusters")

    else:
        cluster_color = _adata_to_color_dict(self.oracle.adata, cluster_column_name)[cluster]

        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s)
        if cluster == "True" :
            ax_.set_title(f"Cells of interest: \n{cluster_column_name}")
    ax_.axis("off")

    ##
    ax_ = ax[1]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####
    ###
    ax_ = ax[2]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    #ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")


    ax_ = ax[3]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
    else:
        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    #ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=original_value, cmap="rainbow", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    plt.show()



    fig, ax = plt.subplots(1, 4, figsize=[20,5])


    #####
    ax_ = ax[0]
    ax_.set_title(f"Perturb simulation \n color: {cluster_column_name}")
    if cluster == "whole":
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=None,
                                   cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    else:
        _plot_quiver_for_a_cluster(oracle=self.oracle,
                                   cluster_column_name=cluster_column_name,
                                   color=cluster_color,
                                   cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    ax_.axis("off")

    ax_ = ax[1]
    ax_.set_title("Perturb simulation \n color: cluster")
    _plot_quiver_for_a_cluster(oracle=self.oracle,
                               cluster_column_name=cluster_column_name,
                               cluster=cluster, quiver_scale=30, ax=ax_, args={"s": s})
    ax_.axis("off")

    ax_ = ax[2]
    ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation result on grid")
    ax_.axis("off")


    ax_ = ax[3]
    if cluster == "whole":
        plot_scatter_with_anndata(adata=self.oracle.adata, obsm_key=self.obsm_key,
                                  cluster_column_name=cluster_column_name,
                                  ax=ax_, args={"s": s})
    else:
        ax_.scatter(whole_embedding[:, 0], whole_embedding[:, 1], c="lightgray", s=s)
        ax_.scatter(original_embedding[:, 0], original_embedding[:, 1], c=cluster_color, s=s, alpha=alpha)
    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.set_title("Perturb simulation result on grid")
    ax_.axis("off")

    plt.show()

    fig, ax = plt.subplots(1, 4, figsize=[20,5])

    ####
    ax_ = ax[0]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[1]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    ax_.scatter(grid[~mass_filter, 0], grid[~mass_filter, 1], c=inner_product[~mass_filter],
                cmap="coolwarm", s=s_grid, vmin=vmin, vmax=vmax)

    ax_.quiver(grid[~mass_filter, 0], grid[~mass_filter, 1],
                    gradient_simulated[~mass_filter, 0], gradient_simulated[~mass_filter, 1],
                    scale=scale_for_simulated, zorder=20000)
    ax_.axis("off")
    ax_.set_title("Inner product of \n Perturb simulation * Development flow \n + Perturb simulation")


    #####
    ax_ = ax[2]
    pcm = ax_.scatter(value_on_grid[~mass_filter], inner_product[~mass_filter],
                     c=inner_product[~mass_filter], cmap="coolwarm",
                      vmin=vmin, vmax=vmax, s=s_grid)

    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.axhline(0, color="lightgray")
    #pp = fig.colorbar(pcm, ax=ax_, orientation="vertical")
    sns.despine()
    ax_.set_xlabel("pseudotime")
    ax_.set_ylabel("inner product score")

    ax_ = ax[3]
    sns.boxplot(data=inner_product_stats, x="pseudotime_id", y="score", color="white", ax=ax_)
    ax_.set_xlabel("Digitized_pseudotime")
    ax_.set_ylabel("inner product score")
    ax_.axhline(0, color="gray")
    ax_.set_ylim([vmin*1.1, vmax*1.1])
    ax_.tick_params(
                labelleft=False)
    plt.show()



    fig, ax = plt.subplots(1, 4, figsize=[20,5])

    stage_grid = self.oracle_dev.stage_grid
    colors_dict = _adata_to_color_dict(adata=self.oracle.adata, cluster_use="Stage")
    stage_order = return_order(stage_list=self.oracle_dev.inner_product_stats.stage)


    ####
    ax_ = ax[0]
    ax_.scatter(grid[mass_filter, 0], grid[mass_filter, 1], s=0)
    plot_grid_with_categprocal_color(grid=grid, mass_filter=mass_filter, color_array=stage_grid,
                                 colors_dict=colors_dict,
                                 ax=ax_, args={})

    ax_.axis("off")
    ax_.set_title("Stage (hpf) on grid")


    ax_ = ax[1]
    plot_legend(labels=stage_order, palette=colors_dict, ax_=ax_)

    #####
    ax_ = ax[2]
    sns.violinplot(data=self.oracle_dev.inner_product_stats,
                   x="pseudotime", y="stage",
                   palette=colors_dict, order=stage_order[::-1],
                   ax=ax_)
    sns.despine()
    #ax_.set_xlabel("pseudotime")
    ax_.set_ylabel("stage (hpf)")

    ax_ = ax[3]
    plot_stackedvar(pd.crosstab(self.oracle_dev.inner_product_stats['stage'],
                                self.oracle_dev.inner_product_stats['pseudotime_id']).loc[stage_order],
                                ax=ax_, palette=colors_dict)
    ax_.set_xlabel("Digitized_pseudotime")
    ax_.set_ylabel("Grid point count")
    sns.despine()
    plt.show()

def plot_grid_with_categprocal_color(grid, mass_filter, color_array, colors_dict, ax, args={}):
    x = grid[~mass_filter, 0]
    y = grid[~mass_filter, 1]
    color_array_filtered = color_array[~mass_filter]

    for cluster, color in colors_dict.items():
        idx = np.where(color_array_filtered == cluster)[0]
        ax.scatter(x[idx], y[idx], c=color, label=cluster, **args)

def return_order(stage_list):
    stage_unique = np.unique(stage_list)
    hpf_unique = [float(i[:4]) for i in stage_unique]
    stage_order = list(stage_unique[np.argsort(hpf_unique)])
    return stage_order

def plot_legend(labels, palette, ax_):

    for i, label in enumerate(labels):
        ax_.scatter([0], [i], s=100, c=palette[label])
        ax_.text(1, i-len(labels)*0.015, s=label)
    ax_.set_ylim([-1, len(labels)])
    ax_.set_xlim([-1, 10])
    ax_.axis("off")

def plot_stackedvar(df, ax, palette=None):

    bottom_feats=[]
    if palette is None:
        for i, j in enumerate(df.index.values):
            if i==0:
                ax.bar(df.columns.values, df.loc[j].values, edgecolor='white', label=j)
            else:
                ax.bar(df.columns.values, df.loc[j].values, label=j,
                        bottom=df.loc[bottom_feats].sum(axis=0).values,
                        edgecolor='white')
            bottom_feats.append(j)
    else:
        for i, j in enumerate(df.index.values):
            if i==0:
                ax.bar(df.columns.values, df.loc[j].values,
                    edgecolor='white', color=palette[j], label=j)
            else:
                ax.bar(df.columns.values, df.loc[j].values, label=j, color=palette[j],
                        bottom=df.loc[bottom_feats].sum(axis=0).values,
                        edgecolor='white')
            bottom_feats.append(j)
    #plt.legend()
    ax.set_xticks(df.columns)
