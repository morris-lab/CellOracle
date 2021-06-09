# -*- coding: utf-8 -*-
"""
The :mod:`.data` module implements data download and loading.

"""

from .load_data import (load_TFinfo_df_mm9_mouse_atac_atlas,
                        load_Paul2015_data,
                        load_tutorial_links_object,
                        load_tutorial_oracle_object)

__all__ = ["load_TFinfo_df_mm9_mouse_atac_atlas", "load_Paul2015_data",
           "load_tutorial_links_object",
           "load_tutorial_oracle_object"]
