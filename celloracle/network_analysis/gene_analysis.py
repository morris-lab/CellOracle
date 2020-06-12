# -*- coding: utf-8 -*-
'''
This is a series of custom functions for the inferring of GRN from single cell RNA-seq data.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tqdm.notebook import tqdm

from .cartography import plot_cartography_kde

#import seaborn as sns

settings = {"save_figure_as": "png"}

def plot_scores_as_rank(links, cluster, n_gene=50, save=None):
    """
    Pick up top n-th genes wich high-network scores and make plots.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        cluster (str): Cluster nome to analyze.
        n_gene (int): Number of genes to plot. Default is 50.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """

    col = ['degree_centrality_all',
   'degree_centrality_in', 'degree_centrality_out',
   'betweenness_centrality',  'eigenvector_centrality', "participation"]
    for value in col:

        res = links.merged_score[links.merged_score.cluster == cluster]
        res = res[value].sort_values(ascending=False)
        res = res[:n_gene]

        fig = plt.figure()

        plt.scatter(res.values, range(len(res)))
        plt.yticks(range(len(res)), res.index.values)#, rotation=90)
        plt.xlabel(value)
        plt.title(f" {value} \n top {n_gene} in {cluster}")
        plt.gca().invert_yaxis()
        plt.subplots_adjust(left=0.5, right=0.85)

        if not save is None:
            os.makedirs(save, exist_ok=True)
            path = os.path.join(save, f"ranked_values_in_{links.name}_{value}_{links.thread_number}_in_{cluster}.{settings['save_figure_as']}")
            fig.savefig(path, transparent=True)
        plt.show()


def _plot_goi(x, y, goi, args_annot, scatter=False, x_shift=0.1, y_shift=0.1):
    """
    Plot an annoation to highlight one point in scatter plot.

    Args:
        x (float): Cordinate-x.
        y (float): Cordinate-y.
        args_annot (dictionary): arguments for matplotlib.pyplot.annotate().
        scatter (bool): Whether to plot dot for the point of interest.
        x_shift (float): distance between the annotation and the point of interest in the x-axis.
        y_shift (float): distance between the annotation and the point of interest in the y-axis.
    """

    default = {"size": 10}
    default.update(args_annot)
    args_annot = default.copy()

    arrow_dict = {"width": 0.5, "headwidth": 0.5, "headlength": 1, "color": "black"}
    if scatter:
        plt.scatter(x, y, c="none", edgecolor="black")
    plt.annotate(f"{goi}", xy=(x, y), xytext=(x+x_shift, y+y_shift),
                 color="black", arrowprops=arrow_dict, **args_annot)



def plot_score_comparison_2D(links, value, cluster1, cluster2, percentile=99, annot_shifts=None, save=None, fillna_with_zero=True):
    """
    Make a scatter plot that shows the relationship of a specific network score in two groups.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        value (srt): The network score to be shown.
        cluster1 (str): Cluster nome to analyze. Network scores in the cluste1 are shown as x-axis.
        cluster2 (str): Cluster nome to analyze. Network scores in the cluste2 are shown as y-axis.
        percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 99.
        annot_shifts ((float, float)): Shift x and y cordinate for annotations.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    res = links.merged_score[links.merged_score.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="index")
    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))

    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index

    gois = np.union1d(goi1, goi2)

    x, y = piv[cluster1], piv[cluster2]
    plt.scatter(x, y, c="none", edgecolor="black")

    if annot_shifts is None:
        x_shift, y_shift = (x.max() - x.min())*0.03, (y.max() - y.min())*0.03
    else:
        x_shift, y_shift = annot_shifts
    for goi in gois:
        x, y = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        _plot_goi(x, y, goi, {}, scatter=False, x_shift=x_shift, y_shift=y_shift)

    plt.xlabel(cluster1)
    plt.ylabel(cluster2)
    plt.title(f"{value}")
    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save, f"values_comparison_in_{links.name}_{value}_{links.thread_number}_{cluster1}_vs_{cluster2}.{settings['save_figure_as']}")
        plt.savefig(path, transparent=True)
    plt.show()




def _test_ver_plot_score_comparison_2D(links, value, cluster1, cluster2, percentile=99,
                             annot_shifts=None, save=None, fillna_with_zero=True,
                             plot_line=True, threshold = 0.25, figsize=(4, 4)):
    from adjustText import adjust_text

    """
    Make a scatter plot that shows the relationship of a specific network score in two groups.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        value (srt): The network score to be shown.
        cluster1 (str): Cluster nome to analyze. Network scores in the cluste1 are shown as x-axis.
        cluster2 (str): Cluster nome to analyze. Network scores in the cluste2 are shown as y-axis.
        percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 99.
        annot_shifts ((float, float)): Shift x and y cordinate for annotations.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    def rotation_n45(s):
        k = 1/np.sqrt(2)
        x = k*(s[0] + s[1])
        y = k*(-s[0] + s[1])
        return x, y

    def distance_(s):
        x, y = rotation_n45(s)
        return abs(y)

    def rotation_p45(s):
        k = 1/np.sqrt(2)
        x = k*(s[0] - s[1])
        y = k*(s[0] + s[1])
        return x, y

    def distance_1(s):
        x, y = rotation_p45(s)
        return abs(y)

    fig, ax = plt.subplots(figsize=figsize)
    plt.subplots_adjust(left=0.15)
    plt.subplots_adjust(bottom=0.15)


    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    res = links.merged_score[links.merged_score.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="index")
    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))

    # pick up genes with high score
    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index
    gois = np.union1d(goi1, goi2)


    # plot lines
    x, y = piv[cluster1], piv[cluster2]
    if plot_line:
        # y = x line
        max_ = max(x.max(), y.max())
        plt.plot((0, max_),(0, max_), linestyle=(0, (5, 10)), color="gray", linewidth=0.5)
        # other lines
        plt.plot((0, max_-threshold),(threshold, max_),
                 linestyle="dashed", color="gray", linewidth=0.5, alpha=1)
        plt.plot((threshold, max_), (0, max_-threshold),
                 linestyle="dashed", color="gray", linewidth=0.5, alpha=1)
        plt.plot((0, threshold),(threshold, 0),
                    linestyle="dashed", color="gray", linewidth=0.5, alpha=1)

    # plot all gene
    plt.scatter(x, y, c="lightgray")

    # plot genes that have high scores
    piv["dis_"] = piv.apply(distance_1, axis=1).values
    piv_ = piv[piv.dis_ >= threshold/np.sqrt(2)]
    x, y = piv_[cluster1], piv_[cluster2]
    plt.scatter(x, y, c="dimgray")

    # plot genes that are differentially connected
    piv_["dis"] = piv_.apply(distance_, axis=1).values
    x, y = piv_[piv_.dis>=threshold/np.sqrt(2)][cluster1], piv_[piv_.dis>=threshold/np.sqrt(2)][cluster2]
    gois = np.concatenate([gois, x.index])
    plt.scatter(x, y, c="indianred")

    # plot annotations
    texts = []
    gois = np.unique(gois)
    for goi in gois:
        x, y = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        plt_text = ax.annotate(goi, (x, y), fontsize=8, color='black')
        texts.append(plt_text)

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='dimgray', lw=0.5))

    plt.xlabel(cluster1)
    plt.ylabel(cluster2)
    plt.title(f"{value}")
    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save, f"values_comparison_in_{links.name}_{value}_{links.thread_number}_{cluster1}_vs_{cluster2}.{settings['save_figure_as']}")
        plt.savefig(path, transparent=True)
    plt.show()


######################
### score dynamics ###
######################
def plot_score_per_cluster(links, goi, save=None):
    """
    Plot network score for a specific gene.
    This function can be used to compare network score of a specific gene between clusters
    and get insight about the dynamics of the gene.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        goi (srt): Gene name.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    print(goi)
    res = links.merged_score[links.merged_score.index==goi]
    res = res.rename(
        columns={"degree_centrality_all": "degree\ncentrality",
                 "betweenness_centrality": "betweenness\ncentrality",
                 "eigenvector_centrality": "eigenvector\ncentrality"})
    # make plots
    values = [ "degree\ncentrality",  "betweenness\ncentrality",
              "eigenvector\ncentrality"]
    for i, value in zip([1, 2, 3], values):
        plt.subplot(1, 3,  i)
        ax = sns.stripplot(data=res, y="cluster", x=value,
                      size=10, orient="h",linewidth=1, edgecolor="w",
                      order=links.palette.index.values,
                      palette=links.palette.palette.values)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.yaxis.grid(True)
        ax.tick_params(bottom=False,
                        left=False,
                        right=False,
                        top=False)
        if i > 1:
            plt.ylabel(None)
            ax.tick_params(labelleft=False)

    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save,
                           f"score_dynamics_in_{links.name}_{links.thread_number}_{goi}.{settings['save_figure_as']}")
        plt.savefig(path, transparent=True)
    plt.show()


###################
### cartography ###
###################

def plot_cartography_scatter_per_cluster(links, gois=None, clusters=None,
                                         scatter=False, kde=True,
                                         auto_gene_annot=False, percentile=98,
                                         args_dot={}, args_line={},
                                         args_annot={}, save=None):
    """
    Plot gene network cartography in a scatter or kde plot.
    Please read the original paper of gene network cartography for detail.
    https://www.nature.com/articles/nature03288

    Args:
        links (Links object): See network_analisis.Links class for detail.
        gois (list of srt): List of Gene name to highlight.
        clusters (list of str): List of cluster name to analyze. If None, all clusters in Links object will be analyzed.
        scatter (bool): Whether to make scatter plot.
        auto_gene_annot (bool): Whether to pick up genes to make annotation.
        percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 98.
        args_dot (dictionary): Arguments for scatter plot.
        args_line (dictionary): Arguments for lines in cartography plot.
        args_annot (dictionary): Arguments for annoation in plots.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    if clusters is None:
        clusters = links.cluster

    for cluster in clusters:
        print(cluster)
        data = links.merged_score[links.merged_score.cluster == cluster]
        data = data[["connectivity", "participation"]]
        #data = data[["within_module_degree", "participation_coefficient"]]

        if auto_gene_annot:
            goi1 = data[data["connectivity"] > np.percentile(data["connectivity"].values, percentile)].index
            goi2 = data[data["participation"] > np.percentile(data["participation"].values, percentile)].index

            if gois is None:
                gois_ = np.union1d(goi1, goi2)
            else:
                gois_ = np.union1d(gois, goi1)
                gois_ = np.union1d(gois_, goi2)
        else:
            gois_ = gois

        fig = plt.figure()

        plot_cartography_kde(data, gois_, scatter, kde,
                             args_dot, args_line, args_annot)
        plt.title(f"cartography in {cluster}")
        plt.subplots_adjust(left=0.2, bottom=0.25)

        if not save is None:
            os.makedirs(save, exist_ok=True)
            path = os.path.join(save, f"cartography_in_{links.name}_{links.thread_number}_{cluster}.{settings['save_figure_as']}")
            fig.savefig(path, transparent=True)
        plt.show()


def plot_cartography_term(links, goi, save=None):
    """
    Plot the summary of gene network cartography like a heatmap.
    Please read the original paper of gene network cartography for detail.
    https://www.nature.com/articles/nature03288

    Args:
        links (Links object): See network_analisis.Links class for detail.
        gois (list of srt): List of Gene name to highlight.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    print(goi)
    tt = pd.get_dummies(links.merged_score[["cluster", "role"]],columns=["role"])
    tt = tt.loc[goi].set_index("cluster")
    tt.columns = [i.replace("role_", "") for i in tt.columns]

    order = ["Ultra peripheral", "Peripheral", "Connector","Kinless","Provincical Hub","Connector Hub", "Kinless Hub"]

    tt = tt.reindex(index=links.palette.index.values, columns=order).fillna(0)

    sns.heatmap(data=tt, cmap="Blues", cbar=False)
    if not save is None:
        os.makedirs(save, exist_ok=True)
        path = os.path.join(save,
                           f"cartography_role_in_{links.name}_{links.thread_number}_{goi}.{settings['save_figure_as']}")
        plt.savefig(path, transparent=True)
    plt.show()
