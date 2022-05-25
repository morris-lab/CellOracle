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

from tqdm.auto import tqdm

from ..network.net_core import Net
from ..utility import standard
from .links_object import Links
from ..trajectory.oracle_utility import _adata_to_df, _get_clustercolor_from_anndata, _check_color_information_and_create_if_not_found


RIDGE_SOLVER = "auto"


########################################
### Construct cluster specific GRNs  ###
########################################

def get_links(oracle_object, cluster_name_for_GRN_unit=None, alpha=10, bagging_number=20, verbose_level=1, test_mode=False, model_method="bagging_ridge", n_jobs=-1):
    """
    Make GRN for each cluster and returns results as a Links object.
    Several preprocessing should be done before using this function.

    Args:
        oracle_object (Oracle): See Oracle module for detail.

        cluster_name_for_GRN_unit (str): Cluster name for GRN calculation. The cluster information should be stored in Oracle.adata.obs.

        alpha (float or int): The strength of regularization.
            If you set a lower value, the sensitivity increases, and you can detect weaker network connections. However, there may be more noise.
            If you select a higher value, it will reduce the chance of overfitting.

        bagging_number (int): The number used in bagging calculation.


        verbose_level (int): if [verbose_level>1], most detailed progress information will be shown.
            if [1 >= verbose_level > 0], one progress bar will be shown.
            if [verbose_level == 0], no progress bar will be shown.

        test_mode (bool): If test_mode is True, GRN calculation will be done for only one cluster rather than all clusters.

        model_method (str): Chose modeling algorithm. "bagging_ridge" or "bayesian_ridge"

        n_jobs (int): Number of cpu cores for parallel calculation.  -1 means using all available cores. Default is -1.

    """
    if model_method not in ["bagging_ridge", "bayesian_ridge"]:
        raise ValueError("model_mothod error. Please set 'bagging_ridge' or 'bayesian_ridge'.")

    if cluster_name_for_GRN_unit is None:
        cluster_name_for_GRN_unit = oracle_object.cluster_column_name

    # calculate GRN for each cluster
    linkLists = _fit_GRN_for_network_analysis(oracle_object, cluster_name_for_GRN_unit=cluster_name_for_GRN_unit,
                                  alpha=alpha, bagging_number=bagging_number,  verbose_level=verbose_level, test_mode=test_mode,
                                  model_method=model_method, n_jobs=n_jobs)

    # initiate links object
    links = Links(name=cluster_name_for_GRN_unit,
                 links_dict=linkLists)

    # extract color infomation
    # update color information
    _check_color_information_and_create_if_not_found(adata=oracle_object.adata,
                                                     cluster_column_name=cluster_name_for_GRN_unit,
                                                     embedding_name=oracle_object.embedding_name)

    links.palette = _get_clustercolor_from_anndata(adata=oracle_object.adata,
                                                   cluster_name=cluster_name_for_GRN_unit,
                                                   return_as="palette")

    #links.merge_links()
    links.ALPHA_used = alpha

    links.model_method = model_method

    return links


def _fit_GRN_for_network_analysis(oracle_object, cluster_name_for_GRN_unit, alpha=10, bagging_number=20,
                                  verbose_level=1, test_mode=False, model_method="bagging_ridge", n_jobs=-1):

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
                print(f"Inferring GRN for {cluster}...")

            cells_in_the_cluster_bool = (cluster_info == cluster)
            gem_ = gem_imputed[cells_in_the_cluster_bool]
            gem_std = gem_imputed_std[cells_in_the_cluster_bool]


            tn_ = Net(gene_expression_matrix=gem_,
                         gem_standerdized=gem_std,
                         TFinfo_dic=oracle_object.TFdict,
                         verbose=False)
            tn_.fit_All_genes(bagging_number=bagging_number,
                              model_method=model_method,
                              alpha=alpha,
                              verbose=verbose,
                              n_jobs=n_jobs)


            #oracle_object.linkMat[cluster] = tn_.returnResultAs_TGxTFs("coef_abs")
            tn_.updateLinkList(verbose=False)
            linkLists[cluster] = tn_.linkList.copy()

    return linkLists
