# -*- coding: utf-8 -*-
'''
This is a series of custom functions for the inferring of GRN from single cell RNA-seq data.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

import sys
import os
import glob
from copy import deepcopy
import pickle
from tqdm.auto import tqdm
import seaborn as sns


from joblib import dump, load


##########################
### internal functions ###
##########################

def _plot_base_cartography(z_min, z_max, args_line={}):


    default = {"linestyle": "dashed",
               "alpha": 0.5,
               "c": "gray"}
    default.update(args_line)
    args_line = default.copy()

    plt.plot([-0.05, 1.05], [2.5, 2.5], **args_line)
    plt.plot([0.05, 0.05], [z_min, 2.5], **args_line)
    plt.plot([0.62, 0.62], [z_min, 2.5], **args_line)
    plt.plot([0.8, 0.8], [z_min, 2.5], **args_line)
    plt.plot([0.3, 0.3], [2.5, z_max + 0.5], **args_line)
    plt.plot([0.75, 0.75], [2.5, z_max + 0.5], **args_line)

    plt.xlabel("Participation coefficient (P)")
    plt.ylabel("Whithin-module degree (z)")

    plt.xlim([-0.1, 1.1])


def __plot_goi(data, goi, args_annot, scatter=False):


    x = data.loc[goi, "p"]
    y = data.loc[goi, "z"]
    if scatter:
        plt.scatter(x, y, c="none", edgecolor="black")

    __plot_annotations(x, y, goi, args=args_annot)


@np.vectorize
def __plot_annotations(x, y, label, x_shift=0.05, y_shift=0.05,
                       args={}):

    args_annot = {"size": 10, "color": "black"}
    args_annot.update(args)

    arrow_dict = {"width": 0.5, "headwidth": 0.5, "headlength": 1, "color": "black"}

    plt.annotate(label, xy=(x, y), xytext=(x+x_shift, y+y_shift),
                 arrowprops=arrow_dict, **args_annot)



def add_label_cartography(data):

    # categorize nodes
    R1_nodes = data[(data.p <= 0.05) & (data.z < 2.5)].index
    R2_nodes = data[(data.p > 0.05) & (data.p <= 0.62) & (data.z < 2.5)].index
    R3_nodes = data[(data.p > 0.62) & (data.p <= 0.80) & (data.z < 2.5)].index
    R4_nodes = data[(data.p > 0.80) & (data.z < 2.5)].index
    R5_nodes = data[(data.p <= 0.30) & (data.z >= 2.5)].index
    R6_nodes = data[(data.p > 0.30) & (data.p <= 0.75) & (data.z >= 2.5)].index
    R7_nodes = data[(data.p > 0.75) & (data.z >= 2.5)].index

    node_kind = ["R1: Ultra peripheral", "R2: Peripheral", "R3: Non-hub",
                 "R4: Non-hub kinless", "R5: Provincial",
                 "R6: Connector hubs", "R7: Kinless hubs"]

    data["label_cartography"] = "0"

    for i in range(1, 8):
        data.loc[eval(f"R{i}_nodes"), "label_cartography"] = node_kind[(i-1)]

    return data


##################################
### functions for making plots ###
##################################

### scatter plots
def plot_cartography_color_scatter(data_, gois=None, args_dot={}, args_line={},
                     args_annot={}):

    data = data_.copy()
    data.columns = ["z", "p"]
    data = add_label_cartography(data)

    _plot_base_cartography(z_min=0, z_max=data.z.max(),
                           args_line=args_line)

    labels = list(data.label_cartography.unique())
    labels.sort()
    for i in labels:
        tmp_data = data[data.label_cartography == i]
        plt.scatter(tmp_data.p,
                    tmp_data.z,
                    label=i, **args_dot)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

    if not gois is None:
        for i in gois:
            __plot_goi(data, i, args_annot)

def plot_cartography_kde(data_, gois=None, scatter=False, kde=True,
                             args_kde={},
                             args_line={},
                             args_annot={}):
    """

    kind: str
        "kde", "hex"
    """
    default = {"cmap": "RdPu",
               "shade":True,
               "shade_lowest":False,
               "n_levels": 12}#"RdPu"
    default.update(args_kde)
    args_kde = default.copy()

    default = {"alpha": 0.5, "c":"white"}
    default.update(args_line)
    args_line = default.copy()

    if not args_kde["shade"]:
        args_line.update({"c": "gray"})


    data = data_.copy()
    data = data.dropna(axis=0)

    data.columns = ["z", "p"]
    #data = data[(data.p>=0) & (data.z >= 0)]

    data = add_label_cartography(data)

    _plot_base_cartography(z_min=-1, z_max=data.z.max(), args_line=args_line)

    labels = list(data.label_cartography.unique())
    labels.sort()

    x, y = data.p, data.z,
    if kde:
        sns.kdeplot(x, y, **args_kde)

    if scatter:
        plt.scatter(x, y,  marker='o', edgecolor="lightgray",
                    c="none", alpha=1)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

    if not gois is None:
        gois = np.intersect1d(data.index.values, gois)
        for i in gois:
            __plot_goi(data, i, args_annot, True)
    plt.xlabel("Participation coefficient (P)")
    plt.ylabel("Whithin-module\ndegree (z)")
