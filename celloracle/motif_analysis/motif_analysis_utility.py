# -*- coding: utf-8 -*-
'''

'''

###########################
### 0. Import libralies ###
###########################


# 0.1. libraries for fundamental data science and data processing
import inspect

import pandas as pd
import numpy as np

import logging

#from tqdm.auto import tqdm
from tqdm.auto import tqdm
# 0.2. libraries for DNA and genome data wrangling and Motif analysis
from genomepy import Genome
import genomepy

#from gimmemotifs.motif import Motif
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta
from gimmemotifs.config import DIRECT_NAME, INDIRECT_NAME
from gimmemotifs import __version__ as gmotif_version

from .process_bed_file import peak_M1

from .reference_genomes import SUPPORTED_REF_GENOME
from ..utility.package_version_checker import _is_version_OK

####
###
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
    

def get_available_ref_genome_info(species, provider="UCSC"):

    df = pd.DataFrame(np.array(list(genomepy.list_available_genomes(provider))))
    species_ = df[2].apply(lambda x: x.split(" ")[0])

    return df[species_.isin([species])]



def is_genome_installed(ref_genome, genomes_dir=None):
    """
    Celloracle motif_analysis module uses gimmemotifs and genomepy internally.
    Reference genome files should be installed in the PC to use gimmemotifs and genomepy.
    This function checks the installation status of the reference genome.

    Args:
        ref_genome (str): names of reference genome. i.e., "mm10", "hg19"
        genomes_dir (str): Installation directory of Genomepy reference genome data.

    """
    try:
        genome_data = Genome(name=ref_genome, genomes_dir=genomes_dir)

        return True

    except:
        if ref_genome == "AmexG_v6.0-DD":
            print(f"genome {ref_genome} is not installed in this environment.")
            print(f"Please install {ref_genome} data with genomepy using the following command in command line (terminal).")
            print(f"Please make sure to use recent version of genomepy and have enough space on your PC. Otherwise you may get an error.")
            print(f'genomepy install https://www.axolotl-omics.org/dl/AmexG_v6.0-DD.fa.gz url')

        else:
            if  ref_genome in SUPPORTED_REF_GENOME.ref_genome.values:
                provider = SUPPORTED_REF_GENOME["provider"][SUPPORTED_REF_GENOME.ref_genome==ref_genome].values[0]
            else:
                provider = "PROVIDER"

            print(f"genome {ref_genome} is not installed in this environment.")
            print("Please install genome using genomepy.")
            if genomes_dir is None:
                print(f'e.g.\n    >>> import genomepy\n    >>> genomepy.install_genome(name="{ref_genome}", provider="{provider}")')
            else:
                print(f'e.g.\n    >>> import genomepy\n    >>> genomepy.install_genome(name="{ref_genome}", provider="{provider}", genomes_dir={genomes_dir})')

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

def scan_dna_for_motifs(scanner_object, motifs_object, sequence_object, divide=100000, verbose=True, batch_size=None):
    '''
    This is a wrapper function to scan DNA sequences searchig for Gene motifs.

    Args:

        scanner_object (gimmemotifs.scanner): Object that do motif scan.

        motifs_object (gimmemotifs.motifs): Object that stores motif data.

        sequence_object (gimmemotifs.fasta): Object that stores sequence data.

    Returns:
        pandas.dataframe: scan results is stored in data frame.

    '''
    logging.debug(f"batch_size: {batch_size}")

    seqnames = []
    motif_ids = []
    factors_directs = []
    factors_indirects = []
    scores = []
    poss = []
    strands = []
   
    
    if batch_size is None:
        batch_size = 50000

    # This is alternative version of scanner_object.scan() function.
    iter = scan_modified(scanner_object=scanner_object, seqs=sequence_object, verbose=verbose, batch_size=batch_size)

    logging.debug(f"Scan started.")
    for i, result in enumerate(iter):
        seqname = sequence_object.ids[i]
        for m,matches in enumerate(result):
            motif = motifs_object[m]
            for score, pos, strand in matches:
                seqnames.append(seqname)
                motif_ids.append(motif.id)
                factors_directs.append(list2str(motif.factors[DIRECT_NAME]))
                factors_indirects.append(list2str(motif.factors[INDIRECT_NAME]))
                scores.append(score)
                poss.append(pos)
                strands.append(strand)
                
    logging.debug(f"Scan finished. Making dataframe..")

    if len(seqnames)==0:
        df = pd.DataFrame(columns=["seqname",
                               "motif_id",
                               "factors_direct",
                               "factors_indirect",
                               "score", "pos", "strand"])
    else:

        df = pd.DataFrame({"seqname": seqnames,
                           "motif_id": motif_ids,
                           "factors_direct": factors_directs,
                           "factors_indirect": factors_indirects,
                           "score": scores,
                           "pos": poss,
                           "strand": strands})
        
        logging.debug(f"Adjusting data types..")
        df.score = df.score.astype(float)
        df.pos = df.pos.astype(int)
        df.strand = df.strand.astype(int)

        # keep peak id name
        df.seqname = list(map(peak_M1, df.seqname.values))
        df.reset_index(drop=True, inplace=True)

    return df



from gimmemotifs.utils import as_fasta

def scan_modified(scanner_object, seqs, nreport=100, scan_rc=True, zscore=False, gc=False, verbose=True, batch_size=50000):
    """
    This is modified version of gimmemotifs.scanner.Scanner.scan() function.
    1. Added batch_size control.
    2. Added progressbar control.
    """
    seqs = as_fasta(seqs, genome=scanner_object.genome)
    if zscore:
        if gc:
            if len(scanner_object.meanstd) <= 1:
                scanner_object.set_meanstd(gc=gc)
        else:
            if len(scanner_object.meanstd) != 1:
                scanner_object.set_meanstd(gc=gc)

    # progress bar
    pbar = tqdm(
        desc="Scanning",
        unit=" sequences",
        total=len(seqs),
        disable=(not verbose),  # can be silenced
    )

    # Check version of gimmemotifs.
    sig = inspect.signature(scanner_object._scan_sequences)

    # For the recent gimmemotifs version.
    if ['nreport', 'scan_rc', 'seqs', 'zscore'] ==  sorted(list(sig.parameters.keys())):
        for batch_idx in range(0, len(seqs), batch_size):
            it = scanner_object._scan_sequences(
                seqs=seqs.seqs[batch_idx : batch_idx + batch_size],
                nreport=nreport,
                scan_rc=scan_rc,
                zscore=zscore,
            )
            for result in it:
                yield result
                pbar.update(1)

    # For the old gimmemotifs version, 0.14.x
    elif ['nreport', 'scan_rc', 'seqs'] ==  sorted(list(sig.parameters.keys())):
        for batch_idx in range(0, len(seqs), batch_size):
            it = scanner_object._scan_sequences(
                seqs=seqs.seqs[batch_idx : batch_idx + batch_size],
                nreport=nreport,
                scan_rc=scan_rc,
            )
            for result in it:
                yield result
                pbar.update(1)
    else:
        raise ValueError(f"Please check gimmemotifs version. Your version is {gmotif_version}.")
    pbar.close()



