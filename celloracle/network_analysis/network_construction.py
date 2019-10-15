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

import numpy as np
import pandas as pd

from scipy import stats

from tqdm import tqdm_notebook as tqdm

from ..network.net_core import Net
from ..utility import standard
from .links_object import Links
from ..trajectory.oracle_utility import _adata_to_df

RIDGE_SOLVER = "auto"


########################################
### Construct cluster specific GRNs  ###
########################################

def get_links(oracle_object, cluster_name_for_GRN_unit=None, alpha=10, bagging_number=20, verbose_level=1, test_mode=False):
    """
    Make GRN for each cluster and returns results as a Links object.
    Several preprocessing should be done before using this function.

    Args:
        oracle_object (Oracle): See Oracle module for detail.

        cluster_name_for_GRN_unit (str): Cluster name for GRN calculation. The cluster information should be stored in Oracle.adata.obs.

        alpha (float or int): the strength of regularization.
            If you set a lower value, the sensitivity increase, and you can detect a weak network connection, but it might get more noize.
            With a higher value of alpha may reduce the chance of overfitting.

        bagging_number (int): The number for bagging calculation.


        verbose_level (int): if [verbose_level>1], most detailed progress information will be shown.
            if [verbose_level > 0], one progress bar will be shown.
            if [verbose_level == 0], no progress bar will be shown.

        test_mode (bool): If test_mode is True, GRN calculation will be done for only one cluster rather than all clusters.

    """
    if cluster_name_for_GRN_unit is None:
        cluster_name_for_GRN_unit = oracle_object.cluster_column_name

    # calculate GRN for each cluster
    linkLists = _fit_GRN_for_network_analysis(oracle_object, cluster_name_for_GRN_unit=cluster_name_for_GRN_unit,
                                  alpha=alpha, bagging_number=bagging_number,  verbose_level=verbose_level, test_mode=test_mode)

    # initiate links object
    links = Links(name=cluster_name_for_GRN_unit,
                 links_dict=linkLists)
    if not test_mode:
        # extract color infomation
        c = [i.upper() for i in oracle_object.adata.uns[f"{cluster_name_for_GRN_unit}_colors"]]
        cnames = list(oracle_object.adata.obs[cluster_name_for_GRN_unit].unique())
        cnames.sort()
        links.add_palette(cnames, c)

    #links.merge_links()
    links.ALPHA_used = alpha

    return links


def _fit_GRN_for_network_analysis(oracle_object, cluster_name_for_GRN_unit, alpha=10, bagging_number=20,
                                  verbose_level=1, test_mode=False):

    # extract information from oracle_object
    gem_imputed = _adata_to_df(oracle_object.adata, "imputed_count")
    gem_imputed_std = standard(gem_imputed)
    cluster_info = oracle_object.adata.obs[cluster_name_for_GRN_unit]
    linkLists = {}

    # setting about verbose
    if verbose_level < 0:
        raise ValueError("varbose_level should be positive number.")
    elif verbose_level == 0:
        loop = np.unique(cluster_info)
        verbose = False
    else:
        loop = tqdm(np.unique(cluster_info))
        if verbose_level <= 1:
            verbose = False
        else:
            verbose = True

    First = True
    for cluster in loop:
        if (not test_mode) | First:
            First = False
            if verbose:
                print(f"inferring GRN for {cluster}...")

            cells_in_the_cluster_bool = (cluster_info == cluster)
            gem_ = gem_imputed[cells_in_the_cluster_bool]
            gem_std = gem_imputed_std[cells_in_the_cluster_bool]


            tn_ = Net(gene_expression_matrix=gem_,
                         gem_standerdized=gem_std,
                         TFinfo_dic=oracle_object.TFdict,
                         verbose=False)
            tn_.fit_All_genes(bagging_number=bagging_number,
                              alpha=alpha, verbose=verbose)


            #oracle_object.linkMat[cluster] = tn_.returnResultAs_TGxTFs("coef_abs")
            tn_.updateLinkList(verbose=False)
            linkLists[cluster] = tn_.linkList.copy()

    return linkLists
