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

from tqdm.notebook import tqdm

# 0.2. libraries for DNA and genome data wrangling and Motif analysis
from genomepy import Genome

#from gimmemotifs.motif import Motif
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta

from pybedtools import BedTool


####
### bed f
def decompose_chrstr(peak_str):
    """
    Take peak name as input and return splitted strs.

    Args:
        peak_str (str): peak name.

    Returns:
        tuple: splitted peak name.

    Examples:
       >>> decompose_chrstr("chr1_111111_222222")
       "chr1", "111111", "222222"
    """
    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)

    return chr_, start, end

def list_peakstr_to_df(x):
    """
    Convert list of peaks(str) into data frame.

    Args:
       x (list of str): list of peak names

    Returns:
       pandas.dataframe: peak info as DataFrame

    Examples:
       >>> x = ['chr1_3094484_3095479', 'chr1_3113499_3113979', 'chr1_3119478_3121690']
       >>> list_peakstr_to_df(x)
                   chr	start	end
            0	chr1	3094484	3095479
            1	chr1	3113499	3113979
            2	chr1	3119478	3121690
    """
    df = np.array([decompose_chrstr(i) for i in x])
    df = pd.DataFrame(df, columns=["chr", "start", "end"])
    df["start"] = df["start"].astype(np.int)
    df["end"] = df["end"].astype(np.int)
    return df

def df_to_list_peakstr(x):

    x = x.rename(columns={"chrom":"chr"})

    peak_str = x.chr + "_" + x.start.astype(str) + "_" + x.end.astype(str)
    peak_str = peak_str.values

    return peak_str

####
###



def peak_M1(peak_id):
    """
    Take a peak_id (index of bed file) as input,
    then subtract 1 from Start position.

    Args:
        peak_id (str): Index of bed file. It should be made of "chromosome name", "start position", "end position"
            e.g. "chr11_123445555_123445577"
    Returns:
        str: Processed peak_id.

    Examples:
        >>> a = "chr11_123445555_123445577"
        >>> peak_M1(a)
        "chr11_123445554_123445577"
    """
    chr_, start, end = decompose_chrstr(peak_id)
    return chr_ + "_" + str(int(start)-1) + "_" + end


def peak2fasta(peak_ids, ref_genome):

    '''
    Convert peak_id into fasta object.

    Args:
        peak_id (str or list of str): Peak_id.  e.g. "chr5_0930303_9499409"
            or it can be a list of peak_id.  e.g. ["chr5_0930303_9499409", "chr11_123445555_123445577"]

        ref_genome (str): Reference genome name.   e.g. "mm9", "mm10", "hg19" etc

    Returns:
        gimmemotifs fasta object: DNA sequence in fasta format

    '''
    genome_data = Genome(ref_genome)

    def peak2seq(peak_id):
        chromosome_name, start, end = decompose_chrstr(peak_id)
        locus = (int(start),int(end))

        tmp = genome_data[chromosome_name][locus[0]:locus[1]]
        name = f"{tmp.name}_{tmp.start}_{tmp.end}"
        seq = tmp.seq
        return (name, seq)


    if type(peak_ids) is str:
        peak_ids = [peak_ids]

    fasta = Fasta()
    for peak_id in peak_ids:
        name, seq = peak2seq(peak_id)
        fasta.add(name, seq)

    return fasta

def remove_zero_seq(fasta_object):
    """
    Remove DNA sequence with zero length
    """
    fasta = Fasta()
    for i, seq in enumerate(fasta_object.seqs):
        if seq:
            name = fasta_object.ids[i]
            fasta.add(name, seq)
    return fasta


def read_bed(bed_path):
    """
    Load bed file and return as dataframe.

    Args:
        bed_path (str): File path.

    Returns:
        pandas.dataframe: bed file in dataframe.

    """
    tt = BedTool(bed_path).to_dataframe().dropna(axis=0)
    tt["seqname"] = tt.chrom + "_" + tt.start.astype("int").astype("str") + "_" + tt.end.astype("int").astype("str")
    return tt
