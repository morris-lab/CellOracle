# -*- coding: utf-8 -*-
"""
The :mod:`.chip_analysis` module implements transcription factor motif scan.

"""
from .integrate_chip_data_with_grn_results import intersect_chip_bed_with_atac_bed, load_and_merge_chip_bed_files, integrate_Chip_and_TN, loadAndProcessData, integrate_Chip_and_TN_for_HM_analysus, loadAndProcessData_for_HM_analysis
from .make_plots import get_meta_data, pickUpData, plot_ROC, binarizeChip, filterGRNscore


__all__ = ["intersect_chip_bed_with_atac_bed",
           "load_and_merge_chip_bed_files",
           "integrate_Chip_and_TN",
           "loadAndProcessData",
           "integrate_Chip_and_TN_for_HM_analysus",
           "loadAndProcessData_for_HM_analysis",
           "get_meta_data",
           "pickUpData",
           "plot_ROC",
           "binarizeChip",
           "filterGRNscore"]
