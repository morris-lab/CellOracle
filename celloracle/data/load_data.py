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

import pandas as pd
import numpy as np
import scanpy as sc
import sys, os
from ..data import __path__ as parent_path
from ..data.config import CELLORACLE_DATA_DIR, WEB_PAR_DIR
from ..utility import load_hdf5
from ..utility.data_download_from_web import download_data_if_data_not_exist


def load_mouse_scATAC_atlas_base_GRN(version="0.10.0", force_download=False):
    """
    Load Transcription factor binding information made from mouse scATAC-seq atlas dataset.
    mm9 genome was used for the reference genome.

    Args:

    Returns:
        pandas.dataframe: TF binding info.
    """
    if version == "0.9.0":
        filename = "TFinfo_data/mm9_mouse_atac_atlas_data_TSS_and_cicero_0.9_accum_threshold_10.5_DF_peaks_by_TFs.parquet"
    elif version == "0.10.0":
        filename = "TFinfo_data/mm9_mouse_atac_atlas_data_TSS_and_cicero_0.9_accum_threshold_10.5_DF_peaks_by_TFs_v202204.parquet"
    else:
        raise ValueError("This version is not found.")

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)

load_TFinfo_df_mm9_mouse_atac_atlas = load_mouse_scATAC_atlas_base_GRN # Old function name

def load_Paul2015_data(version="0.10.0", force_download=False):
    """
    Load mouse hematopoiesis scRNA-seq data. The data was processed according to the standard CellOracle preprocessing method described in the tutorial.
    """
    if version == "0.9.0":
        filename = "anndata/Paul_etal_v201212.h5ad"
    elif version == "0.10.0":
        filename = "anndata/Paul_etal_v202204.h5ad"
    else:
        raise ValueError("This version is not found.")

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return sc.read_h5ad(path)

def load_tutorial_links_object(version="0.10.0", force_download=False):
    """
    """
    if version == "0.9.0":
        filename = "tutorial_data/links_louvain_190829.celloracle.links"
    elif version == "0.10.0":
        filename = "tutorial_data/links_louvain_v20220406.celloracle.links"
    else:
        raise ValueError("This version is not found.")
    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return load_hdf5.load_hdf5(path)

def load_tutorial_oracle_object(version="0.10.0", force_download=False):
    if version == "0.9.0":
        filename = "tutorial_data/Paul_etal_v20210704.celloracle.oracle"
    elif version == "0.10.0":
        filename = "tutorial_data/Paul_etal_v20220406.celloracle.oracle"
    else:
        raise ValueError("This version is not found.")

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)
    return load_hdf5.load_hdf5(path)
