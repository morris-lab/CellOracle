# -*- coding: utf-8 -*-
'''
This is a series of custom functions for the inferring of GRN from single cell RNA-seq data.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing
import numpy as np
import pandas as pd

# 0.2. libraries for Machine learning modeling
from sklearn.ensemble import BaggingRegressor
from sklearn.linear_model import Ridge, BayesianRidge


# 0.3. custom libraries
from ..utility import intersect

###############################
### bayesian_ridge function ###
###############################

def get_bayesian_ridge_coefs(target_gene, gem, gem_scaled, TFdict, cellstate=None,
                             scaling=True):

    #print(target_gene)
    ## 1. Data prep
    if target_gene not in TFdict.keys():
        #print("err")
        return 0, 0, 0

    # define regGenes
    reggenes = TFdict[target_gene]
    allgenes_detected = list(gem.columns)
    reggenes = intersect(reggenes, allgenes_detected)

    if target_gene in reggenes:
        reggenes.remove(target_gene)

    reg_all = reggenes.copy()

    if not reggenes: # if reqgene is empty, return 0
        return 0, 0, 0

    # add cell state data
    if not cellstate is None:
        cellstate_name = list(cellstate.columns)
        reg_all += cellstate_name


    # prepare learning data
    if scaling:
        data = gem_scaled[reg_all]
    else:
        data = gem[reg_all]
    label = gem[target_gene]

    # make model and fitting
    model_br = BayesianRidge()
    model_br.fit(data.values, label)

    coef_names = reggenes
    ind_ = np.where([i in reggenes for i in reg_all])

    coef_mean = model_br.coef_[ind_]
    coef_variance = model_br.sigma_.diagonal()[ind_]

    return coef_mean, coef_variance, coef_names

##############################
### bagging_ridge function ###
##############################

def get_bagging_ridge_coefs(target_gene, gem, gem_scaled, TFdict, cellstate=None,
                 bagging_number=1000, scaling=True, n_jobs=-1, alpha=1, solver="auto"):
    #print(target_gene)
    ## 1. Data prep
    if target_gene not in TFdict.keys():
        #print("err")
        return 0

    # define regGenes
    reggenes = TFdict[target_gene]
    allgenes_detected = list(gem.columns)
    reggenes = intersect(reggenes, allgenes_detected)

    if target_gene in reggenes:
        reggenes.remove(target_gene)

    reg_all = reggenes.copy()

    if not reggenes: # if reqgene is empty, return 0
        return 0

    # add cell state data
    if not cellstate is None:
        cellstate_name = list(cellstate.columns)
        reg_all += cellstate_name


    # prepare learning data
    if scaling:
        data = gem_scaled[reg_all]
    else:
        data = gem[reg_all]
    label = gem[target_gene]

    #print(n_jobs)
    # bagging model
    model = BaggingRegressor(base_estimator=Ridge(alpha=alpha,
                                                  solver=solver,
                                                  random_state=123),
                             n_estimators=bagging_number,
                             bootstrap=True,
                             max_features=0.8,
                             n_jobs=n_jobs,
                             verbose=False,
                             random_state=123)
    model.fit(data, label)

    # get results
    coefs = _get_coef_matrix(model, reg_all)

    # remove cellstate data from coefs
    coefs = coefs[intersect(reggenes, coefs.columns.values)]

    return coefs


# this is a function to extract coef information from sklearn ensemble_model.
def _get_coef_matrix(ensemble_model, feature_names):
    # ensemble_model: trained ensemble model. e.g. BaggingRegressor
    # feature_names: list or numpy array of feature names. e.g. feature_names=X_train.columns
    feature_names = np.array(feature_names)
    n_estimater = len(ensemble_model.estimators_features_)
    coef_list = \
        [pd.Series(ensemble_model.estimators_[i].coef_,
                   index=feature_names[ensemble_model.estimators_features_[i]])\
         for i in range(n_estimater)]

    coef_df = pd.concat(coef_list, axis=1, sort=False).transpose()

    return coef_df
