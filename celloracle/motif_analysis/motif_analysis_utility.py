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

#from tqdm.auto import tqdm
from tqdm.auto import tqdm
# 0.2. libraries for DNA and genome data wrangling and Motif analysis
from genomepy import Genome

#from gimmemotifs.motif import Motif
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta
from gimmemotifs.config import DIRECT_NAME, INDIRECT_NAME
from gimmemotifs import __version__ as gmotif_version

from .modified_gimmemotif_function import custom_scan


from .process_bed_file import peak_M1

from .reference_genomes import SUPPORTED_REF_GENOME
from ..utility.package_version_checker import _is_version_OK

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
        if ref_genome == "AmexG_v6.0-DD":
            print(f"genome {ref_genome} is not installed in this environment.")
            print(f"Please install {ref_genome} data with genomepy using the following command in command line (terminal).")
            print(f"Please make sure to use recent version of genomepy and have enough space on your PC. Otherwise you may get an error.")
            print(f'genomepy install -p url https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz')

        else:
            if  ref_genome in SUPPORTED_REF_GENOME.ref_genome.values:    
                provider = SUPPORTED_REF_GENOME["provider"][SUPPORTED_REF_GENOME.ref_genome==ref_genome].values[0]
            else:
                provider = "PROVIDER"

            print(f"genome {ref_genome} is not installed in this environment.")
            print("Please install genome using genomepy.")
            print(f'e.g.\n    >>> import genomepy\n    >>> genomepy.install_genome(name="{ref_genome}", provider="{provider}")')

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
    return ", ".join(li)






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
    # If the gimmemotifs version is larger than 0.17.0, it will automatically return the progress bar,
    # So we don't need to make tqdm qbar here.
    if _is_version_OK(que=gmotif_version, ref='0.17.0'):
        verbose=False

    pbar = tqdm(
        desc="scanning",
        unit=" sequences",
        total=len(sequence_object),
        disable=(verbose==False),  # can be silenced
    )
    iter = scanner_object.scan(sequence_object)

    for i, result in enumerate(iter):
        seqname = sequence_object.ids[i]
        for m,matches in enumerate(result):
            motif = motifs_object[m]
            for score, pos, strand in matches:
                li.append(np.array([seqname,
                                    motif.id,
                                    list2str(motif.factors[DIRECT_NAME]),
                                    list2str(motif.factors[INDIRECT_NAME]),
                                    score, pos, strand]))
        pbar.update(1)
    pbar.close()

    #save_as_pickled_object(li, "./tmp_li.pickle")
    #print("saved tmp list")


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
            df.score = df.score.astype(float)
            df.pos = df.pos.astype(int)
            df.strand = df.strand.astype(int)

            # keep peak id name
            df.seqname = list(map(peak_M1, df.seqname.values))

            LI.append(df)


            if divide*(k+1) >= len(li):
                remaining = 0
            k += 1

        df = pd.concat(LI, axis=0)
    return df


# This is the preliminary version of mini-batch scan function. Keep this function just in case.
'''
def scan_dna_for_motifs_by_batch(scanner_object, motifs_object, sequence_object, batch_size=100, divide=100000, verbose=True):
    """
    This is a wrapper function to scan DNA sequences searchig for Gene motifs.

    Args:

        scanner_object (gimmemotifs.scanner): Object that do motif scan.

        motifs_object (gimmemotifs.motifs): Object that stores motif data.

        sequence_object (gimmemotifs.fasta): Object that stores sequence data.

    Returns:
        pandas.dataframe: scan results is stored in data frame.

    """
    li = scan_by_small_batch(sequence_object=sequence_object,
                             scanner_object=scanner_object,
                             motifs_object=motifs_object,
                             batch_size=batch_size,
                             verbose=verbose)



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
            df.score = df.score.astype(float)
            df.pos = df.pos.astype(int)
            df.strand = df.strand.astype(int)

            # keep peak id name
            df.seqname = list(map(peak_M1, df.seqname.values))

            LI.append(df)


            if divide*(k+1) >= len(li):
                remaining = 0
            k += 1

        df = pd.concat(LI, axis=0)
        df = df.reset_index(drop=True)

    return df

def _gimmemotifs_wrapper(sequence_object, scanner_object, motifs_object):
    li = []
    for i, result in enumerate(custom_scan(self=scanner_object, seqs=sequence_object, verbose=False)):
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

def scan_by_small_batch(sequence_object, scanner_object, motifs_object, batch_size=100, verbose=True):
    func = functools.partial(_gimmemotifs_wrapper, scanner_object=scanner_object, motifs_object=motifs_object)

    li = run_function_by_small_batch(function=func,
                                     input_data_iterable=sequence_object,
                                     batch_size=batch_size,
                                     return_method="add",
                                     verbose=verbose
                                     )
    return li

def run_function_by_small_batch(function=print, input_data_iterable=range(100), batch_size=30, return_method="list", verbose=True):
    """
    Run a function little by little.
    This can be applied to any function.
    This may be useful if a function can cause problem with large data.
    """
    batch_total = len(input_data_iterable) // batch_size
    if len(input_data_iterable) % batch_size != 0:
        batch_total += 1

    pbar = tqdm(
        desc="scanning",
        unit=" sequence_batch",
        total=batch_total,
        disable=(verbose==False),  # can be silenced
    )


    i = 0
    results = []
    while len(input_data_iterable) > batch_size*i:
        small_sequences = input_data_iterable[i*batch_size:batch_size*(i + 1)]
        out = function(small_sequences)
        if return_method == "add":
            results += out
        elif return_method == "list":
            results.append(out)
        pbar.update(1)
        i += 1
    pbar.close()

    return results
'''
