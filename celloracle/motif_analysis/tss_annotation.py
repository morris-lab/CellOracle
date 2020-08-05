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

from tqdm.notebook import tqdm

from pybedtools import BedTool

from ..motif_analysis import __path__ as parent_path
from .process_bed_file import list_peakstr_to_df
from .tfinfo_core import SUPPORTED_REF_GENOME

def _load_tss_ref_data(ref_genome):

    """
    Args:
        ref_genome (str): Reference genome name.
            Please contact us through github issue page if you have a request for another referene genome.
    """
    path = os.path.join(parent_path[0], "tss_ref_data", f"{ref_genome}_tss_info.bed")
    return BedTool(fn=path)




def get_tss_info(peak_str_list, ref_genome, verbose=True):
    """
    Get annotation about Transcription Starting Site (TSS).

    Args:
        peak_str_list (list of str): list of peak_id. e.g., [“chr5_0930303_9499409”, “chr11_123445555_123445577”]
        ref_genome (str): reference genome name.
        verbose (bool): verbosity.
    """
    SUPPORTED_REF_GENOME_LIST = []
    for refs in SUPPORTED_REF_GENOME.values():
        SUPPORTED_REF_GENOME_LIST += refs

    if ref_genome not in SUPPORTED_REF_GENOME_LIST:
        raise ValueError(ref_genome, " is not supported currently. Supported refgenomes are ", SUPPORTED_REF_GENOME)

    ref = _load_tss_ref_data(ref_genome=ref_genome)

    queue = list_peakstr_to_df(peak_str_list)
    queue = BedTool.from_dataframe(queue)

    annotated = annotate_tss(tss_ref_bed=ref, queue_bed=queue, verbose=verbose)

    return annotated

def annotate_tss(tss_ref_bed, queue_bed, verbose=True):
    """
    Search for peaks that exist in TSS.
    If a peak overlap TSS peak, that peak will be annotated with gene name of the TSS.
    If a peak does not overlap any TSS peak, the peak will be removed and will not be contained in returned bed file.

    Args:
        tss_fer_bed (Bedtool): bedtool object of tss data.
            This bed file should be a bed file with 5 columns: ["chrom", "start", "end", "name", "score", "strand"].
            Gene name of TSS should be stored in "name" column.

        queue_bed (Bedtool): bedtool object.
            This bed file should be a bed file with 3 columns: ["chrom", "start", "end"].

    Returns:
        dataframe: dataframe of peaks.
    """

    # check data structure
    if (tss_ref_bed.to_dataframe().columns != ["chrom", "start", "end", "name", "score", "strand"]).all():
        raise ValueError("tss annotation error")

    if (queue_bed.to_dataframe().columns != ["chrom", "start", "end"]).all():
        raise ValueError("queue bed file error")

    # intersect
    intersected = tss_ref_bed.intersect(queue_bed, wb=True)
    intersected = intersected.to_dataframe()

    # change name
    intersected = intersected.rename(columns={"start": "start_intersected",
                                          "end": "end_intersected",
                                          "name": "gene_short_name"})
    intersected = intersected.rename(columns={"thickStart": "chr", "thickEnd": "start",
                                          "itemRgb": "end"})

    # select data
    intersected = intersected[["chr", "start", "end", "gene_short_name", "strand"]]

    if verbose:
        print(f"que bed peaks: {len(queue_bed.to_dataframe())}")
        print(f"tss peaks in que: {len(intersected)}")
    return intersected
