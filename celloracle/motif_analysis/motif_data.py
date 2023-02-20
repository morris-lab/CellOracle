# -*- coding: utf-8 -*-
'''

'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing

import pandas as pd
import numpy as np

import sys, os

from tqdm.auto import tqdm
from glob import glob

# 0.2. libraries for DNA and genome data wrangling and Motif analysis

from gimmemotifs.motif import read_motifs

from ..data import __path__ as parent_path
from ..data.config import CELLORACLE_DATA_DIR, WEB_PAR_DIR
from ..utility.data_download_from_web import download_data_if_data_not_exist



MOTIFS_LIST = ['CisBP_ver2_Arabidopsis_thaliana.pfm',
               'CisBP_ver2_Arabidopsis_thaliana_GENE_ID.pfm',
               'CisBP_ver2_Caenorhabditis_elegans.pfm',
               'CisBP_ver2_Cavia_porcellus.pfm',
               'CisBP_ver2_Danio_rerio.pfm',
               'CisBP_ver2_Drosophila_ananassae.pfm',
               'CisBP_ver2_Drosophila_erecta.pfm',
               'CisBP_ver2_Drosophila_grimshawi.pfm',
               'CisBP_ver2_Drosophila_melanogaster.pfm',
               'CisBP_ver2_Drosophila_mix.pfm',
               'CisBP_ver2_Drosophila_mojavensis.pfm',
               'CisBP_ver2_Drosophila_persimilis.pfm',
               'CisBP_ver2_Drosophila_pseudoobscura.pfm',
               'CisBP_ver2_Drosophila_sechellia.pfm',
               'CisBP_ver2_Drosophila_simulans.pfm',
               'CisBP_ver2_Drosophila_virilis.pfm',
               'CisBP_ver2_Drosophila_willistoni.pfm',
               'CisBP_ver2_Drosophila_yakuba.pfm',
               'CisBP_ver2_Gallus_gallus.pfm',
               'CisBP_ver2_Homo_sapiens.pfm',
               'CisBP_ver2_Mus_musculus.pfm',
               'CisBP_ver2_Rattus_norvegicus.pfm',
               'CisBP_ver2_Saccharomyces_cerevisiae.pfm',
               'CisBP_ver2_Sus_scrofa.pfm',
               'CisBP_ver2_Xenopus_laevis.pfm',
               'CisBP_ver2_Xenopus_tropicalis.pfm',
               'CisBP_ver2_Xenopus_tropicalis_and_Xenopus_laevis.pfm']

def load_motifs(motifs_name, force_download=False):
    """

    Load motifs from celloracle motif database

    Args:
        motifs_name (str) : Name of motifs.

    Returns:
        list : List of gimmemotifs.motif object.


    """

    if motifs_name not in MOTIFS_LIST:
        raise ValueError("The motifs name was not in the list. Available motifs: ", MOTIFS_LIST)

    # Load data from local directory if file exits.
    path = os.path.join(parent_path[0], "motif_data", motifs_name)
    path_factor_data = os.path.join(parent_path[0], "motif_data", motifs_name.replace(".pfm", ".motif2factors.txt"))
    if (force_download == False) & os.path.isfile(path) & os.path.isfile(path_factor_data):
        pass
    else:
        # Download pfm file
        path = os.path.join(CELLORACLE_DATA_DIR, "motif_data", motifs_name)
        backup_url = os.path.join(WEB_PAR_DIR, "motif_data", motifs_name)
        download_data_if_data_not_exist(path=path, backup_url=backup_url)
        # Download motif2factors.txt file
        path_factor_data = os.path.join(CELLORACLE_DATA_DIR, "motif_data", motifs_name.replace(".pfm", ".motif2factors.txt"))
        backup_url = os.path.join(WEB_PAR_DIR, "motif_data", motifs_name.replace(".pfm", ".motif2factors.txt"))
        download_data_if_data_not_exist(path=path_factor_data, backup_url=backup_url)

    motifs = read_motifs(path)

    return motifs


def _make_pwm_file_list(PATH_TO_MOTIF_DATA=None):
    if PATH_TO_MOTIF_DATA is None:
        from ..data import __path__ as data_folder_path
        PATH_TO_MOTIF_DATA = os.path.join(list(data_folder_path)[0], "motif_data")

    PFM_PATHS = glob(os.path.join(PATH_TO_MOTIF_DATA, "*.pfm"))
    PFM_PATHS.sort()
    MOTIFS_PATH_DICT = {os.path.basename(path):path for path in PFM_PATHS}
    MOTIFS_LIST = list(MOTIFS_PATH_DICT.keys())

    return MOTIFS_LIST

####
###
