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
import functools

#from tqdm.notebook import tqdm
from tqdm import tqdm
# 0.2. libraries for DNA and genome data wrangling and Motif analysis
from genomepy import Genome

#from gimmemotifs.motif import Motif
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta
from gimmemotifs.config import DIRECT_NAME, INDIRECT_NAME


from .process_bed_file import peak_M1

from .reference_genomes import SUPPORTED_REF_GENOME


####
###

def get_available_ref_genome_info(species, provider="UCSC"):

    df = pd.DataFrame(np.array(list(genomepy.list_available_genomes(provider))))
    species_ = df[2].apply(lambda x: x.split(" ")[0])

    return df[species_.isin([species])]



def is_genome_installed(ref_genome):
    """
    Celloracle motif_analysis module uses gimmemotifs and genomepy internally.
    Reference genome files should be installed in the PC to use gimmemotifs and genomepy.
    This function checks the installation status of the reference genome.

    Args:
        ref_genome (str): names of reference genome. i.e., "mm10", "hg19"

    """
    try:
        genome_data = Genome(ref_genome)

        return True

    except:
        if  ref_genome in SUPPORTED_REF_GENOME.ref_genome.values:
            provider = SUPPORTED_REF_GENOME["provider"][SUPPORTED_REF_GENOME.ref_genome==ref_genome].values[0]
        else:
            provider = "PROVIDER"

        print(f"genome {ref_genome} is not installed in this environment.")
        print("Please install genome using genomepy.")
        print(f'e.g.\n    >>> import genomepy\n    >>> genomepy.install_genome("{ref_genome}", "{provider}")')

    return False
        #raise ValueError(f"Ref_Genome: {ref_genome} is not available.")
def list2str(li):
    """
    Convert list of str into merged one str with " ,".
    See example below for detail.

    Args:
        li (list of str): list

    Returns:
        str: merged str.

    Examples:
        >>> a = ["a", "b", "c"]
        >>> list2str(a)
        "a, b, c"
    """
    re = ""
    for i, j in enumerate(li):
        if i>0:
            re += ", "
        re += j
    return re

def _list2str(li):
    """
    Convert list of str into merged one str with " ,".
    See example below for detail.

    Args:
        li (list of str): list

    Returns:
        str: merged str.

    Examples:
        >>> a = ["a", "b", "c"]
        >>> list2str(a)
        "a, b, c"
    """
    return ", ".join(li)



def _gimmemotifs_wrapper(sequence_object, scanner_object, motifs_object):
    li = []
    for i, result in enumerate(scanner_object.scan(sequence_object)):
        seqname = sequence_object.ids[i]
        for m,matches in enumerate(result):
            motif = motifs_object[m]
            for score, pos, strand in matches:
                li.append([seqname,
                        motif.id,
                        _list2str(motif.factors[DIRECT_NAME]),
                        _list2str(motif.factors[INDIRECT_NAME]),
                                    score, pos, strand])
    return li

def run_function_by_small_batch(function=print, input_data_iterable=range(100), batch_size=30, return_method="list"):
    """
    Run a function little by little.
    This can be applied to any function.
    This may be useful if a function can cause problem with large data.
    """


    i = 0
    results = []
    while len(input_data_iterable) > batch_size*i:
        small_sequences = input_data_iterable[i*batch_size:batch_size*(i + 1)]
        out = function(small_sequences)
        if return_method == "add":
            results += out
        elif return_method == "list":
            results.append(out)
        i += 1
    return results


def scan_dna_for_motifs_by_batch(scanner_object, motifs_object, sequence_object, batch_size=100, divide=100000, verbose=True):
    '''
    This is a wrapper function to scan DNA sequences searchig for Gene motifs.

    Args:

        scanner_object (gimmemotifs.scanner): Object that do motif scan.

        motifs_object (gimmemotifs.motifs): Object that stores motif data.

        sequence_object (gimmemotifs.fasta): Object that stores sequence data.

    Returns:
        pandas.dataframe: scan results is stored in data frame.

    '''
    func_ = functools.partial(_gimmemotifs_wrapper,
                              scanner_object=scanner_object,
                              motifs_object=motifs_object)
    li = run_function_by_small_batch(function=func_,
                                     input_data_iterable=sequence_object,
                                     batch_size=batch_size,
                                     return_method="add")



    if verbose:
        print("Motif Scan finished. Start post processing.")

    if len(li)==0:
        df = pd.DataFrame(columns=["seqname",
                               "motif_id",
                               "factors_direct",
                               "factors_indirect",
                               "score", "pos", "strand"])
    else:# Convert results as a dataframe little by little, and concat together later.

        remaining = 1
        LI = []
        k = 0
        while remaining == 1:
            #print(k)
            tmp_li = li[divide*k:min(len(li), divide*(k+1))]

            df = pd.DataFrame(tmp_li,
                              columns=["seqname",
                                       "motif_id",
                                       "factors_direct",
                                       "factors_indirect",
                                       "score", "pos", "strand"])
            df.score = df.score.astype(np.float)
            df.pos = df.pos.astype(np.int)
            df.strand = df.strand.astype(np.int)

            # keep peak id name
            df.seqname = list(map(peak_M1, df.seqname.values))

            LI.append(df)


            if divide*(k+1) >= len(li):
                remaining = 0
            k += 1

        df = pd.concat(LI, axis=0)
        df = df.reset_index(drop=True)

    return df


def scan_dna_for_motifs(scanner_object, motifs_object, sequence_object, divide=100000, verbose=True):
    '''
    This is a wrapper function to scan DNA sequences searchig for Gene motifs.

    Args:

        scanner_object (gimmemotifs.scanner): Object that do motif scan.

        motifs_object (gimmemotifs.motifs): Object that stores motif data.

        sequence_object (gimmemotifs.fasta): Object that stores sequence data.

    Returns:
        pandas.dataframe: scan results is stored in data frame.

    '''

    li  = []
    if verbose:
        for i, result in enumerate(tqdm(scanner_object.scan(sequence_object))):
            seqname = sequence_object.ids[i]
            for m,matches in enumerate(result):
                motif = motifs_object[m]
                for score, pos, strand in matches:
                    li.append(np.array([seqname,
                                        motif.id,
                                        list2str(motif.factors[DIRECT_NAME]),
                                        list2str(motif.factors[INDIRECT_NAME]),
                                        score, pos, strand]))
    else:
        for i, result in enumerate(scanner_object.scan(sequence_object)):
            seqname = sequence_object.ids[i]
            for m,matches in enumerate(result):
                motif = motifs_object[m]
                for score, pos, strand in matches:
                    li.append(np.array([seqname,
                                        motif.id,
                                        list2str(motif.factors[DIRECT_NAME]),
                                        list2str(motif.factors[INDIRECT_NAME]),
                                        score, pos, strand]))

    #save_as_pickled_object(li, "./tmp_li.pickle")
    #print("saved tmp list")

    if verbose:
        print("Motif Scan finished. Start post processing.")

    if len(li)==0:
        df = pd.DataFrame(columns=["seqname",
                               "motif_id",
                               "factors_direct",
                               "factors_indirect",
                               "score", "pos", "strand"])
    else:

        remaining = 1
        LI = []
        k = 0
        while remaining == 1:
            #print(k)
            tmp_li = li[divide*k:min(len(li), divide*(k+1))]

            tmp_li = np.stack(tmp_li)
            df = pd.DataFrame(tmp_li,
                              columns=["seqname",
                                       "motif_id",
                                       "factors_direct",
                                       "factors_indirect",
                                       "score", "pos", "strand"])
            df.score = df.score.astype(np.float)
            df.pos = df.pos.astype(np.int)
            df.strand = df.strand.astype(np.int)

            # keep peak id name
            df.seqname = list(map(peak_M1, df.seqname.values))

            LI.append(df)


            if divide*(k+1) >= len(li):
                remaining = 0
            k += 1

        df = pd.concat(LI, axis=0)
    return df
