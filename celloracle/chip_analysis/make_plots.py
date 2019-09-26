# -*- coding: utf-8 -*-
'''
This file contains custom functions for the analysis of ATAC-seq data.
Genomic activity information (peak of ATAC-seq) will be extracted first.
Then the peak DNA sequence will be subjected to TF motif scan.
Finally we will get list of TFs that potentially binds to a specific gene.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc



#################################################
### functions for visualization and filtering ###
#################################################

def get_meta_data(DATA_PATH):
    
    def add_number(list_):
        li = []
        dup = {}
        for i in list_:
            if i in dup.keys():
                i_new = i + "_" + str(dup[i])
                dup[i] += 1
            else:
                dup[i] = 1
                i_new = i
            li.append(i_new)
        return li


    Data = {}

    # file name
    Data["file_names"] = os.listdir(DATA_PATH)
    Data["file_names"] = [i for i in Data["file_names"] if "tn" in i]
    Data["file_names"].sort()

    # file path
    Data["file_path"] = [os.path.join(DATA_PATH, i) for i in Data["file_names"]]

    # sample name
    tissues = [i.split("_")[0] for i in Data["file_names"]]
    li = []
    for i, j in enumerate(tissues):
        if "-10X" in j:
            li.append(j)
        elif "-sorting" in j:
            li.append(j)
        else:
            li.append(j + "-microwell")
    Data["scRNA-seq_data"] = li
    Data["scRNA-seq_data"] = add_number(Data["scRNA-seq_data"])

    # sample type
    Data["sample_type"] = [i.split("-microwell")[0].split("-10X")[0].split("-sorting")[0] for i in Data["scRNA-seq_data"]]
    li = []
    for i in Data["sample_type"]:
        if "Small-" in i:
            li.append(i.split("Small-")[1])
        else:
            li.append(i)
    Data["sample_type"] = li


    # platform
    li = []
    for i in Data["scRNA-seq_data"]:
        if "-10X" in i:
            li.append("10x_chromium")
        elif "-sorting" in i:
            li.append("sorting")
        elif "-microwell" in i:
            li.append("microwell")
    Data["platform"] = li

    data = pd.DataFrame(Data)
    return data
    
    
def pickUpData(meta_df, chip_data_name_list, sample_type_list):
    data = meta_df.copy()
    tmp_data = data[data.sample_type.isin(sample_type_list)]
    li = []
    for i in  tmp_data.index:
        df = pd.read_parquet(data.loc[i, "file_path"])
        df = df[df.sample_id.isin(chip_data_name_list)]
        df["scRNA-seq_data"] = data.loc[i, "scRNA-seq_data"]
        li.append(df)
    dataset = pd.concat(li, axis=0)
    
    return dataset


def plot_ROC(df, plot_value="max_score", plotRandom=True, by="scRNA-seq_data", plot=True):
    
    '''
    Coef; Ceficient [data frame]
    gene_list; gene list [list or numpy array]
    guide_name; gene name [str (when single guide) or list (multiple guides)]
    No_of_random_genes; [int]
    seed; random seed [int]
    shade; plotting with/without shade [bool]
    '''
    
    group = df[by].unique()
    aucs = []
    for i in group:
        df_ = df[df[by] == i]

        # data prep for ROC
        y_label = df_["chip_bind"].values
        vals = df_[plot_value].values
        fp, tp, _ = roc_curve(y_label, vals)
        auc_score = auc(fp, tp)
        aucs.append(auc_score)

        # plot curve
        if plot:
            plt.plot(fp, tp, lw=2, 
                     label=f"{i}, auc: {auc_score:.2g}")

    # plot landom line
    if plot:
        if plotRandom:
            plt.plot([0, 1], [0, 1], lw=2, linestyle='--', label="random")
    
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.0])
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.title('ROC')
        plt.legend(loc="lower right")
        
    aucs = np.array(aucs)
    aucs = pd.DataFrame(aucs, index=group).transpose()
    return aucs

def binarizeChip(df, threshold_chip_score_order):
    df_ = df.copy()
    df_["chip_bind"] = (df_.chip_score_order <= threshold_chip_score_order).astype("int") 
    return df_

def filterGRNscore(df, threshold_p_max_score):
    df_ = df.copy()
    
    if "-logp" in df.columns:
        df_["p-threshold"] = (df_["-logp"] >= threshold_p_max_score).astype("int")
    elif "max_score" in df.columns:
        df_["p-threshold"] = (df_["max_score"] >= threshold_p_max_score).astype("int")
    df_["coef_abs_filtered"] = df_["coef_abs"]*df_["p-threshold"]
    
    return df_
