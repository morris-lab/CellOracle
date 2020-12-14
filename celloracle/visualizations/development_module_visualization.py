# -*- coding: utf-8 -*-


import os, sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from .config import CONFIG

def plot_cluster_whole(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):

    if ax is None:
        ax = plt

    ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.colorandum, s=s, **args)
    ax.axis("off")

def plot_cluster_cells_use(self, ax=None, s=CONFIG["s_scatter"], color=None, show_background=True, args=CONFIG["default_args"]):

    if ax is None:
        ax = plt

    if s == 0:
        color = "white"

    if show_background:
        plot_background(self=self, ax=ax, s=s, args=args)

    if not hasattr(self, "cell_idx_use"):
        self.cell_idx_use = None

    if self.cell_idx_use is None:
        if color is None:
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.colorandum, s=s, **args)
        else:
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=color, s=s, **args)

    else:
        if color is None:
            ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                       c=self.colorandum[self.cell_idx_use, :],s=s, **args)
        else:
            ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                       c=color, s=s, **args)


    ax.axis("off")


def plot_background(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):

    if ax is None:
        ax = plt

    ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=s, **args)

    #ax.set_title("Pseudotime")
    ax.axis("off")

def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"], show_background=True, cmap="rainbow", args=CONFIG["default_args"]):

    if ax is None:
        ax = plt

    if show_background:
        plot_background(self=self, ax=ax, s=s, args=args)

    if self.cell_idx_use is None:
        ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.pseudotime, cmap=cmap, s=s, **args)
    else:
        ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                    c=self.pseudotime[self.cell_idx_use], cmap=cmap, s=s, **args)

    ax.axis("off")


def plot_background_on_grid(self, ax=None, s=CONFIG["s_grid"], args={}):


    if ax is None:
        ax = plt

    if hasattr(self, "mass_filter_whole_reference"):
        mass_filter = self.mass_filter_whole_reference
    elif hasattr(self, "mass_filter_whole"):
        mass_filter = self.mass_filter_whole

    ax.scatter(self.gridpoints_coordinates[:, 0],
               self.gridpoints_coordinates[:, 1], s=0)

    if "c" not in args.keys():
        ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               c="lightgray", s=s, **args)
    else:
        ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               s=s, **args)
    ax.axis("off")


def plot_pseudotime_on_grid(self, ax=None, s=CONFIG["s_grid"], show_background=True, args={}):

    if ax is None:
        ax = plt

    if hasattr(self, "mass_filter_simulation"):
        mass_filter = self.mass_filter_simulation
    elif hasattr(self, "mass_filter"):
        mass_filter = self.mass_filter

    if show_background:
        plot_background_on_grid(self=self, ax=ax, s=s, args=args)
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color="white", show_background=False, args={})


    ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               c=self.pseudotime_on_grid[~mass_filter],
               cmap="rainbow", s=s, **args)

    ax.axis("off")



def plot_reference_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):

    if ax is None:
        ax = plt

    if hasattr(self, "mass_filter_simulation"):
        mass_filter = self.mass_filter_simulation
    elif hasattr(self, "mass_filter"):
        mass_filter = self.mass_filter

    if show_background:
        plot_background(self=self, ax=ax, s=s, args=CONFIG["default_args"])
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color="white", show_background=False, args={})

    ax.quiver(self.gridpoints_coordinates[~mass_filter, 0],
              self.gridpoints_coordinates[~mass_filter, 1],
              self.ref_flow[~mass_filter, 0],
              self.ref_flow[~mass_filter, 1],
              scale=scale, **args)

    ax.axis("off")
def plot_simulation_flow_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
    _plot_simulation_flow_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, data_random=False, args=args)

def plot_simulation_flow_random_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
    _plot_simulation_flow_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, data_random=True, args=args)

def _plot_simulation_flow_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], data_random=False, args=CONFIG["default_args_quiver"]):

    if ax is None:
        ax = plt

    if show_background:
        plot_background(self=self, ax=ax, s=s, args=CONFIG["default_args"])
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color=None, show_background=False, args={})

    # mass filter selection
    if hasattr(self, "mass_filter_simulation"):
        mass_filter = self.mass_filter_simulation
    elif hasattr(self, "mass_filter"):
        mass_filter = self.mass_filter

    # Gridpoint cordinate selection
    if hasattr(self, "gridpoints_coordinates"):
        gridpoints_coordinates = self.gridpoints_coordinates
    elif hasattr(self, "mass_filter"):
        gridpoints_coordinates = self.flow_grid

    # Arrow selection
    if data_random:
        flow = self.flow_rndm
    else:
        flow = self.flow

    ax.quiver(gridpoints_coordinates[~mass_filter, 0],
              gridpoints_coordinates[~mass_filter, 1],
              flow[~mass_filter, 0],
              flow[~mass_filter, 1], #zorder=20000,
              scale=scale, **args)

    ax.axis("off")

def plot_inner_product_on_grid(self, ax=None, vm=1,s=CONFIG["s_grid"], show_background=True, args={}):

    if ax is None:
        ax = plt

    if show_background:
        plot_background_on_grid(self=self, ax=ax, s=s,
                                args={"facecolor": "None",
                                      "c": "None",
                                      "edgecolors":'black',
                                      "linewidths": 0.05})
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color=None, show_background=False, args={})

    ax.scatter(self.gridpoints_coordinates[~self.mass_filter_simulation, 0],
               self.gridpoints_coordinates[~self.mass_filter_simulation, 1],
               c=self.inner_product[~self.mass_filter_simulation],
               cmap="coolwarm", vmin=-vm, vmax=vm, s=s, **args)

    ax.axis("off")



def plot_inner_product_on_pseudotime(self, ax=None, vm=1, s=CONFIG["s_grid"], args={}):

    if ax is None:
        fig, ax = plt.subplots()

    pcm = ax.scatter(self.pseudotime_on_grid[~self.mass_filter_simulation],
                     self.inner_product[~self.mass_filter_simulation],
                     c=self.inner_product[~self.mass_filter_simulation],
                     cmap="coolwarm",
                     vmin=-vm, vmax=vm, s=s, **args)

    ax.set_ylim([-vm*1.1, vm*1.1])
    ax.axhline(0, color="lightgray")
    pp = plt.colorbar(pcm, ax=ax, orientation="vertical")
    sns.despine()
    ax.set_xlabel("pseudotime")
    ax.set_ylabel("inner product score")

def plot_inner_product_as_box(self, ax=None, vm=1, args={}):

    if ax is None:
        fig, ax = plt.subplots()

    sns.boxplot(data=self.inner_product_df, x="pseudotime_id", y="score", color="white", ax=ax)
    ax.set_xlabel("Digitized_pseudotime")
    ax.set_ylabel("inner product score")
    ax.axhline(0, color="gray")
    ax.set_ylim([-vm*1.1, vm*1.1])
    ax.tick_params(
                labelleft=False)
    sns.despine()

def plot_quiver(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
    _plot_quiver(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args, data_random=False)

def plot_quiver_random(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
    _plot_quiver(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args, data_random=True)

def _plot_quiver(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"], data_random=False):

    if ax is None:
        ax = plt

    if not hasattr(self, "cell_idx_use"):
        self.cell_idx_use = None

    if self.cell_idx_use is None:
        ix_choice = np.arange(self.embedding.shape[0])
    else:
        ix_choice = self.cell_idx_use

    # Plot whole cell with lightgray
    if show_background:
        ax.scatter(self.embedding[:, 0], self.embedding[:, 1],
                   c="lightgray", alpha=1, s=s, **args)


    ax.scatter(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
               c="lightgray", alpha=0.2, edgecolor=(0,0,0,1), s=s, **args)



    if color is None:
        color=self.colorandum[ix_choice]

    quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,
                       linewidths=0.25, width=0.0045,edgecolors="k",
                       color=color, alpha=1)

    if data_random:
        quiver = self.delta_embedding_random
    else:
        quiver = self.delta_embedding

    ax.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
              quiver[ix_choice, 0],
              quiver[ix_choice, 1],
              scale=scale, **quiver_kwargs)

    ax.axis("off")




def visualize_development_module_layout_2(self, scale_for_pseudotime=CONFIG["scale_dev"],
    scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):

    if self.name is None:
        name = "Selected lineage"
    else:
        name = self.name

    fig, ax = plt.subplots(3, 4, figsize=[20,15])

    ax_ = ax[0, 0]
    plot_cluster_whole(self, ax=ax_, s=s)
    ax_.set_title("Whole population")


    ##
    ax_ = ax[0, 1]
    plot_cluster_cells_use(self, ax=ax_, s=s, color="#EC7063", show_background=show_background)
    ax_.set_title(f"{name}")

    ####
    ax_ = ax[0, 2]
    plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)
    ax_.set_title("Pseudotime")

    ###
    ax_ = ax[0, 3]
    plot_reference_flow_on_grid(self, ax=ax_, scale=scale_for_pseudotime, show_background=show_background, s=s)
    ax_.set_title("Development flow")

    ####
    ax_ = ax[1, 0]
    plot_quiver(self, ax=ax_, scale=scale_for_simulation, color="#EC7063", s=s, show_background=show_background)
    ax_.set_title(f"Perturb simulation \n color: {name}")

    ####
    ax_ = ax[1, 1]
    plot_quiver(self, ax=ax_, scale=scale_for_simulation, color=None, s=s, show_background=show_background)
    ax_.set_title("Perturb simulation \n color: cluster")


    ax_ = ax[1, 2]
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=show_background, s=s)
    ax_.set_title("Perturb simulation")


    #####
    ax_ = ax[1, 3]
    plot_cluster_cells_use(self, ax=ax_, s=s, color="#EC7063", show_background=show_background)
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=False, s=s)
    ax_.set_title("Perturb simulation")


    ax_ = ax[2, 0]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")

    ax_ = ax[2, 1]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=False, s=s)
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[2, 2]
    plot_inner_product_on_pseudotime(self, ax=ax_, vm=vm, s=s_grid)

    ax_ = ax[2, 3]
    plot_inner_product_as_box(self, ax=ax_, vm=vm)


def visualize_development_module_layout_1(self, scale_for_pseudotime=CONFIG["scale_dev"],
    scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):


    fig, ax = plt.subplots(2, 4, figsize=[20,10])

    ax_ = ax[0, 0]
    plot_cluster_cells_use(self, ax=ax_, s=s, show_background=show_background)

    ##
    ax_ = ax[0, 1]
    plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)

    ####
    ax_ = ax[0, 2]
    plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background)
    ax_.set_title("Pseudotime on grid")

    ###
    ax_ = ax[0, 3]
    plot_reference_flow_on_grid(self, ax=ax_, scale=scale_for_pseudotime, show_background=show_background, s=s)
    ax_.set_title("Development flow")


    ####
    ax_ = ax[1, 0]
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=show_background, s=s)
    ax_.set_title("Perturb simulation")

    ####
    ax_ = ax[1, 1]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[1, 2]
    plot_inner_product_on_pseudotime(self, ax=ax_, vm=vm, s=s_grid)

    #####
    ax_ = ax[1, 3]
    plot_inner_product_as_box(self, ax=ax_, vm=vm)

def visualize_development_module_layout_0(self, scale_for_pseudotime=CONFIG["scale_dev"],
    scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):


    fig, ax = plt.subplots(2, 3, figsize=[20, 13.5])

    ax_ = ax[0, 0]
    plot_cluster_cells_use(self, ax=ax_, s=s, show_background=show_background)
    ax_.set_title("Cluster")

    ax_ = ax[0, 1]
    plot_reference_flow_on_grid(self, ax=ax_, scale=scale_for_pseudotime, show_background=show_background, s=s)
    ax_.set_title("Development flow")

    ##
    ax_ = ax[0, 2]
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=show_background, s=s)
    ax_.set_title("Perturb simulation")

    ####
    ax_ = ax[1, 0]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    ax_.set_title("Inner product of \n Perturb simulation * Development flow")


    ax_ = ax[1, 1]
    plot_inner_product_on_pseudotime(self, ax=ax_, vm=vm, s=s_grid)

    #####
    ax_ = ax[1, 2]
    plot_inner_product_as_box(self, ax=ax_, vm=vm)





#####

'''

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

'''
