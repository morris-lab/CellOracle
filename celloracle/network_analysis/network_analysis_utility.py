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

from scipy import stats

from tqdm.notebook import tqdm
import networkx as nx


#import seaborn as sns

settings = {"save_figure_as": "png"}


def transfer_scores_from_links_to_adata(adata, links, method="median"):
    """
    Transfer the summary of network scores (median or mean) per group from Links object into adata.

    Args:
        adata (anndata): anndata

        links (Links): Likns object

        method (str): The method to summarize data.
    """
    cluster_name = links.name
    index_by_cluster = adata.obs[cluster_name].values

    if links.entropy is None:
        links.get_network_entropy()

    if method == "median":
        grouped = links.merged_score.groupby(by="cluster").median()
        grouped_entropy = links.entropy.groupby("cluster").median()
    elif method == "mean":
        grouped = links.merged_score.groupby(by="cluster").mean()
        grouped_entropy = links.entropy.groupby("cluster").mean()


    grouped = grouped.loc[index_by_cluster]
    grouped_entropy= grouped_entropy.loc[index_by_cluster]
    grouped = pd.concat([grouped, grouped_entropy], axis=1)
    grouped.index = adata.obs.index

    adata.obs = pd.concat([adata.obs, grouped], axis=1)
#####################
### Use Network X ###
#####################

def linkList_to_networkgraph(filteredlinkList):
    """
    Convert linkList into Graph object in NetworkX.

    Args:
       filteredlinkList (pandas.DataFrame): GRN saved as linkList.

    Returns:
        Graph object: Network X graph objenct.
    """
    G=nx.DiGraph()
    G_ = nx.from_pandas_edgelist(filteredlinkList, edge_attr=True)
    G.add_edges_from(G_.edges())

    return G


def draw_network(linkList, return_graph=False):
    """
    Plot network graph.

    Args:
       linkList (pandas.DataFrame): GRN saved as linkList.
       return_graph (bool): Whether to return graph object.

    Returns:
        Graph object: Network X graph objenct.
    """
    G = linkList_to_networkgraph(linkList)
    # レイアウトの取得
    pos = nx.spring_layout(G)
    # 可視化
    #nx.draw_networkx_edges(G, pos)
    nx.draw(G, pos,edge_color='black',width=1,linewidths=1,\
node_size=500,node_color='pink',alpha=0.9,arrowstyle='->',
            labels={node:node for node in G.nodes()})
    plt.axis('off')
    plt.show()

    if return_graph:
        return G
