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


def load_drosophila_promoter_base_GRN(version="dm6_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"dm3_CisBPv2_fpr1": "promoter_base_GRN/dm3_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "dm3_CisBPv2_fpr2": "promoter_base_GRN/dm3_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "dm6_CisBPv2_fpr1": "promoter_base_GRN/dm6_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "dm6_CisBPv2_fpr2": "promoter_base_GRN/dm6_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)

def load_rat_promoter_base_GRN(version="rn6_gimmemotifsv5_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"rn4_gimmemotifsv5_fpr1": "promoter_base_GRN/rn4_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "rn4_gimmemotifsv5_fpr2": "promoter_base_GRN/rn4_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               "rn5_gimmemotifsv5_fpr1": "promoter_base_GRN/rn5_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "rn5_gimmemotifsv5_fpr2": "promoter_base_GRN/rn5_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               "rn6_gimmemotifsv5_fpr1": "promoter_base_GRN/rn6_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "rn6_gimmemotifsv5_fpr2": "promoter_base_GRN/rn6_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)


def load_mouse_promoter_base_GRN(version="mm10_gimmemotifsv5_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"mm9_gimmemotifsv5_fpr1": "promoter_base_GRN/mm9_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "mm9_gimmemotifsv5_fpr2": "promoter_base_GRN/mm9_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               "mm10_gimmemotifsv5_fpr1": "promoter_base_GRN/mm10_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "mm10_gimmemotifsv5_fpr2": "promoter_base_GRN/mm10_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)

def load_human_promoter_base_GRN(version="hg19_gimmemotifsv5_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"hg19_gimmemotifsv5_fpr1": "promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "hg19_gimmemotifsv5_fpr2": "promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               "hg38_gimmemotifsv5_fpr1": "promoter_base_GRN/hg38_TFinfo_dataframe_gimmemotifsv5_fpr1_threshold_10_20210630.parquet",
               "hg38_gimmemotifsv5_fpr2": "promoter_base_GRN/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)


def load_chicken_promoter_base_GRN(version="galGal6_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"galGal4_CisBPv2_fpr1": "promoter_base_GRN/galGal4_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "galGal4_CisBPv2_fpr2": "promoter_base_GRN/galGal4_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "galGal5_CisBPv2_fpr1": "promoter_base_GRN/galGal5_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "galGal5_CisBPv2_fpr2": "promoter_base_GRN/galGal5_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "galGal6_CisBPv2_fpr1": "promoter_base_GRN/galGal6_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "galGal6_CisBPv2_fpr2": "promoter_base_GRN/galGal6_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)


def load_zebrafish_promoter_base_GRN(version="danRer11_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"danRer7_CisBPv2_fpr1": "promoter_base_GRN/danRer7_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "danRer7_CisBPv2_fpr2": "promoter_base_GRN/danRer7_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "danRer10_CisBPv2_fpr1": "promoter_base_GRN/danRer10_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "danRer10_CisBPv2_fpr2": "promoter_base_GRN/danRer10_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "danRer11_CisBPv2_fpr1": "promoter_base_GRN/danRer11_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "danRer11_CisBPv2_fpr2": "promoter_base_GRN/danRer11_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)

def load_xenopus_tropicalis_promoter_base_GRN(version="xenTro3_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"xenTro2_CisBPv2_fpr1": "promoter_base_GRN/xenTro2_TFinfo_dataframe_CisBPv2_Xenopus_tropicalis_and_Xenopus_laevis_fpr1_threshold_10_20210630.parquet",
               "xenTro2_CisBPv2_fpr2": "promoter_base_GRN/xenTro2_TFinfo_dataframe_CisBPv2_Xenopus_tropicalis_and_Xenopus_laevis_fpr2_threshold_10_20210630.parquet",
               "xenTro3_CisBPv2_fpr1": "promoter_base_GRN/xenTro3_TFinfo_dataframe_CisBPv2_Xenopus_tropicalis_and_Xenopus_laevis_fpr1_threshold_10_20210630.parquet",
               "xenTro3_CisBPv2_fpr2": "promoter_base_GRN/xenTro3_TFinfo_dataframe_CisBPv2_Xenopus_tropicalis_and_Xenopus_laevis_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)


def load_arabidopsis_promoter_base_GRN(version="TAIR10_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"TAIR10_CisBPv2_fpr1": "promoter_base_GRN/TAIR10_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "TAIR10_CisBPv2_fpr2": "promoter_base_GRN/TAIR10_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)



def load_Scerevisiae_promoter_base_GRN(version="sacCer3_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"sacCer2_CisBPv2_fpr1": "promoter_base_GRN/sacCer2_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "sacCer2_CisBPv2_fpr2": "promoter_base_GRN/sacCer2_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "sacCer3_CisBPv2_fpr1": "promoter_base_GRN/sacCer3_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "sacCer3_CisBPv2_fpr2": "promoter_base_GRN/sacCer3_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)

def load_Celegans_promoter_base_GRN(version="ce10_CisBPv2_fpr2", force_download=False):
    """
    Load Base GRN made from promoter DNA sequence and motif scan.

    Args:

    Returns:
        pandas.dataframe: Base GRN as a matrix.
    """

    options = {"ce6_CisBPv2_fpr1": "promoter_base_GRN/ce6_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "ce6_CisBPv2_fpr2": "promoter_base_GRN/ce6_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               "ce10_CisBPv2_fpr1": "promoter_base_GRN/ce10_TFinfo_dataframe_CisBPv2_fpr1_threshold_10_20210630.parquet",
               "ce10_CisBPv2_fpr2": "promoter_base_GRN/ce10_TFinfo_dataframe_CisBPv2_fpr2_threshold_10_20210630.parquet",
               }
    if version in options.keys():
        print(f"Loading prebuilt promoter base-GRN. Version: {version}")
        filename = options[version]
    else:
        print(f"Version error. {version} is not in the list.")
        print("Available option: ", list(options.keys()))

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], filename)
    if (force_download == False) & os.path.isfile(path):
        pass
    else:
        path = os.path.join(CELLORACLE_DATA_DIR, filename)
        backup_url = os.path.join(WEB_PAR_DIR, filename)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)

    return pd.read_parquet(path)
