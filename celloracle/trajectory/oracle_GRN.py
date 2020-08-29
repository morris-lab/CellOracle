# -*- coding: utf-8 -*-

from copy import deepcopy
import warnings
import logging
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import os

from sklearn.linear_model import Ridge
from tqdm.notebook import tqdm

from ..utility import intersect
#from network_fitting import TransNet as tn



def _do_simulation(coef_matrix, simulation_input, gem, n_propagation):
    """
    Silulate signal propagation in GRNs.

    Args:
        coef_matrix (numpy array): 2d matrix that store GRN weights

        simulation_input (pandas.DataFrame): input for simulation

        gem (pandas.DataFrame): input for simulation

        n_propagation (int): number of propagation.

    Returns:
        pandas.DataFrame: simulated gene expression matrix after signal progagation
    """
    delta_input = simulation_input - gem

    delta_simulated = delta_input.copy()
    for i in range(n_propagation):
        delta_simulated = delta_simulated.dot(coef_matrix)
        delta_simulated[delta_input != 0] = delta_input

        # gene expression cannot be negative. adjust delta values to make sure that gene expression are not netavive values.
        gem_tmp = gem + delta_simulated
        gem_tmp[gem_tmp<0] = 0
        delta_simulated = gem_tmp - gem

    gem_simulated = gem + delta_simulated

    return gem_simulated


def _getCoefMatrix(gem, TFdict, alpha=1):
    """
    Calculate GRN and return CoefMatrix (network weights)

    Args:
        gem (pandas.DataFrame): gene expression matrix to calculate GRN

        TFdict (dictionary): python dictionary of potential regulatory gene list

        alpha (float) : strength of regularization in Ridge.

    Returns:
        2d numpy array: numpy array
    """

    genes = gem.columns

    all_genes_in_dict = intersect(gem.columns, list(TFdict.keys()))
    zero_ = pd.Series(np.zeros(len(genes)), index=genes)

    def get_coef(target_gene):
        tmp = zero_.copy()

        # define regGenes
        reggenes = TFdict[target_gene]
        reggenes = intersect(reggenes, genes)

        if target_gene in reggenes:
            reggenes.remove(target_gene)

        if len(reggenes) == 0 :

            tmp[target_gene] = 1
            return(tmp)


        # prepare learning data
        Data = gem[reggenes]
        Label = gem[target_gene]


        # model fitting
        model = Ridge(alpha=alpha, random_state=123)
        model.fit(Data, Label)

        tmp[reggenes] = model.coef_

        return tmp




    li = []
    li_calculated = []

    for i in tqdm(genes):
        if not i in all_genes_in_dict:
            tmp = zero_.copy()
            tmp[i] = 1

        else:
            tmp = get_coef(i)
            li_calculated.append(i)

        li.append(tmp)


    coef_matrix = pd.concat(li, axis=1)

    coef_matrix.columns = genes

    print(f"genes_in_gem: {gem.shape[1]}")
    print(f"models made for {len(li_calculated)} genes")


    return coef_matrix #, li_calculated

def _coef_to_active_gene_list(coef_matrix):
    """
    Args:
        coef_matricx (pd.DataFrame): 2d dataframe (gene x gene). This is a result of GRN calculation.

    Return:
        list: list of active gene that have at least one target gene.
    """
    coef = coef_matrix.copy()
    
    active_TF_list = []
    for i in coef.columns.values:
        coef.loc[i, i] = 0
        if np.any(coef.loc[i, :] != 0):
            active_TF_list.append(i)

    return active_TF_list
