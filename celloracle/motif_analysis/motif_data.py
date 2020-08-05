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

from tqdm.notebook import tqdm
from glob import glob

# 0.2. libraries for DNA and genome data wrangling and Motif analysis


from gimmemotifs.motif import read_motifs

from ..data.motif_data import __path__ as PATH_TO_MOTIF_DATA

PATH_TO_MOTIF_DATA = list(PATH_TO_MOTIF_DATA)[0]
PFM_PATHS = glob(os.path.join(PATH_TO_MOTIF_DATA, "*.pfm"))
PFM_PATHS.sort()
MOTIFS_PATH_DICT = {os.path.basename(path):path for path in PFM_PATHS}
MOTIFS_LIST = list(MOTIFS_PATH_DICT.keys())

def load_motifs(motifs_name):
    """

    Load motifs from celloracle motif database

    Args:
        motifs_name (str) : Name of motifs.

    Returns:
        list : List of gimmemotifs.motif object.


    """

    if motifs_name not in MOTIFS_LIST:
        raise ValueError("The motifs name was not in the list. Available motifs: ", MOTIFS_LIST)

    path = MOTIFS_PATH_DICT[motifs_name]
    motifs = read_motifs(path)

    return motifs


####
###
