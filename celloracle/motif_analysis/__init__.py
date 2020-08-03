# -*- coding: utf-8 -*-
"""
The :mod:`.motif_analysis` module implements transcription factor motif scan.

Genomic activity information (peak of ATAC-seq or Chip-seq) is extracted first.
Then the peak DNA sequence will be subjected to TF motif scan.
Finally we will get list of TFs that potentially binds to a specific gene.

"""

from .motif_analysis_utility import is_genome_installed
from .process_bed_file import peak2fasta, read_bed, remove_zero_seq
from .tfinfo_core import (load_TFinfo, load_TFinfo_from_parquets,
                          make_TFinfo_from_scanned_file, TFinfo, scan_dna_for_motifs,
                          SUPPORTED_REF_GENOME)
from .tss_annotation import get_tss_info
from .process_cicero_data import integrate_tss_peak_with_cicero
from . import process_cicero_data
from .motif_data import load_motifs, PATH_TO_MOTIF_DATA, PFM_PATHS, MOTIFS_PATH_DICT, MOTIFS_LIST


__all__ = ["is_genome_installed", "peak2fasta", "read_bed", "remove_zero_seq",
           "load_TFinfo", "load_TFinfo_from_parquets",
           "make_TFinfo_from_scanned_file",
           "TFinfo", "scan_dna_for_motifs", "SUPPORTED_REF_GENOME",
           "get_tss_info", "process_cicero_data",
           "integrate_tss_peak_with_cicero",
           "load_motifs", "PATH_TO_MOTIF_DATA", "PFM_PATHS", "MOTIFS_PATH_DICT", "MOTIFS_LIST"]
