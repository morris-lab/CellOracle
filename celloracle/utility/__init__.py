# -*- coding: utf-8 -*-
"""
The :mod:`.utility` module has several functions that support celloracle.

"""

from .make_log import makelog

from .utility import (save_as_pickled_object, load_pickled_object,
                      intersect,
                      exec_process,
                      standard, inverse_dictionary,
                      adata_to_color_dict,
                      transfer_all_colors_between_anndata,
                      transfer_color_between_anndata,
                      knn_data_transferer,
                      update_adata)
from .load_hdf5 import load_hdf5

from .pandas_utility_for_jupyternotebook import init_datatable_mode

__all__ = [
           "makelog",
           "save_as_pickled_object",
           "load_pickled_object",
           "intersect",
           "exec_process",
           "standard",
           "load_hdf5",
           "inverse_dictionary"
           "adata_to_color_dict",
           "transfer_all_colors_between_anndata",
           "transfer_color_between_anndata",
           "knn_data_transferer",
           "update_adata"
           ]
