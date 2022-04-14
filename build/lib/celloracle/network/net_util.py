# -*- coding: utf-8 -*-
'''
This is a series of custom functions for the inferring of GRN from single cell RNA-seq data.

Codes were written by Kenji Kamimoto.


'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing
import glob
import os

import numpy as np
import pandas as pd

# 0.3. custom libraries
from ..utility import load_pickled_object
from .net_core import Net

from velocyto.serialization import dump_hdf5, load_hdf5

##############################
### 1.  define main class  ###
##############################


def load_net(file_path):
    return load_hdf5(filename=file_path, obj_class=Net)


def load_net_from_patquets(folder_path):
    """
    Load a Net object that was saved with "save_as_compressed" function.

    Args:
        folder_path (str): Folder path
    """


    obj = load_pickled_object(os.path.join(folder_path, "transNet.pickle"))

    compressed_files_path = glob.glob(os.path.join(folder_path, "*.parquet"))
    compressed_file_name = [i.split("/")[-1] for i in  compressed_files_path]

    for name, path in zip(compressed_file_name, compressed_files_path):
        if name == "gem.parquet":
            obj.gem = pd.read_parquet(path)

        if name == "gem_standerdized.parquet":
            obj.gem_standerdized = pd.read_parquet(path)

        if name == "TFinfo.parquet":
            obj.TFinfo = pd.read_parquet(path)

        if name == "cellstate.parquet":
            obj.cellstate = pd.read_parquet(path)

        if name == "linkList.parquet":
            obj.linkList = pd.read_parquet(path)
    return obj



def getDF_TGxTF(net_object, value_of_interest):
    """
    Extract inferred GRN information and return as a pandas.DataFrame.
    The results was converted as Target gene x TF.


    Args:
        net_object (Net): Net object. GRN calculation have to be done in this object.
        value_of_interest (str): Kind of information to extract.

    """
    if not net_object.fitted_genes:
        print("No model found. Do fit first.")
        return None

    li = []
    genes = np.unique(net_object.fitted_genes)
    for i in genes:
        li.append(net_object.stats_dict[i][[value_of_interest]])

    df = pd.concat(li, axis=1, sort=False)
    df.columns = genes
    df = df.transpose().fillna(0)

    return df

def getDF_peakxTF(net_object, stat_value_of_interest):
    """
    Extract inferred GRN information and return as a pandas.DataFrame.
    The results was converted as peak (genome locus) x TF.

    Args:
        net_object (Net): Net object. GRN calculation have to be done in this object.
        value_of_interest (str): Kind of information to extract.

    """
    # load DF
    DF_grn = getDF_TGxTF(net_object, stat_value_of_interest)
    DF_tfi = trans_net_object.TFinfo.copy()

    # get index and column info
    index_info = DF_tfi[["peak_id", "gene_short_name"]]
    column_info = DF_grn.columns.values

    # convert DF_grn
    DF_grn = DF_grn.reindex(DF_tfi.gene_short_name.values).fillna(0)
    DF_grn = DF_grn.reset_index(drop=True)

    # convert DF_tfi
    DF_tfi = DF_tfi[column_info]

    merge = DF_grn * DF_tfi
    merge = pd.concat([index_info, merge], axis=1)

    merge = merge.set_index("peak_id")
    return merge
