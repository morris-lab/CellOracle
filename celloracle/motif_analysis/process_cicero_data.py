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
#import matplotlib.pyplot as plt
#import seaborn as sns

import sys
import os
import pickle
import glob
import logging
from copy import deepcopy
from datetime import datetime

from tqdm.auto import tqdm


from .process_bed_file import df_to_list_peakstr
from .tss_annotation import get_tss_info

def integrate_tss_peak_with_cicero(tss_peak, cicero_connections):
    """
    Process output of cicero data and returns DNA peak information for motif analysis in celloracle.
    Please see the celloracle tutorial for more information.

    Args:
        tss_peak (pandas.dataframe): dataframe about TSS information. Please use the function, "get_tss_info" to get this dataframe.

        cicero_connections (dataframe): dataframe that stores the results of cicero analysis.

    Returns:
        pandas.dataframe: DNA peak about promoter/enhancer and its annotation about target gene.

    """

    # 1. check tss data format and convert if needed
    tss = tss_peak.copy()
    if np.all([i in tss.columns for i in ["chr", "start", "end"]]):
        tss = pd.DataFrame({"peak_id": df_to_list_peakstr(tss),
                            "gene_short_name": tss.gene_short_name.values})
    elif "peak_id" in tss.columns:
        pass
    else:
        raise ValueError("tss_peak format error")

    # 2. process cicero coaccess data
    cicero_tss = pd.merge(cicero_connections, tss.rename(columns={"peak_id": "Peak1"}), on="Peak1", how="inner")
    cicero_tss = cicero_tss.rename(columns={"Peak2": "peak_id"})[["peak_id", "gene_short_name", "coaccess"]]
    cicero_tss = cicero_tss[cicero_tss.coaccess>0]

    # 3. merge tss data and cicero data
    tss["coaccess"] = 1 # Set coaccecss stcore of peaks in tss as 1 before merging data
    merged = pd.concat([cicero_tss, tss], axis=0)
    merged = merged.groupby(by=["peak_id", "gene_short_name"]).max()
    merged = merged.reset_index()

    return merged
