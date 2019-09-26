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

import sys, os
from ..data import __path__ as parent_path

def load_TFinfo_df_mm9_mouse_atac_atlas():
    """
    Load Transcription factor binding information made from mouse scATAC-seq atlas dataset.
    mm9 genome was used for the reference genome.

    Args:

    Returns:
        pandas.dataframe: TF binding info.
    """

    path = os.path.join(parent_path[0], "TFinfo_data", "mm9_mouse_atac_atlas_data_TSS_and_cicero_0.9_accum_threshold_10.5_DF_peaks_by_TFs.parquet")
    return pd.read_parquet(path)
