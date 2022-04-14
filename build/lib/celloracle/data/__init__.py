# -*- coding: utf-8 -*-
"""
The :mod:`.data` module implements data download and loading.

"""

from .load_data import (load_TFinfo_df_mm9_mouse_atac_atlas,
                        load_mouse_scATAC_atlas_base_GRN,
                        load_Paul2015_data,
                        load_tutorial_links_object,
                        load_tutorial_oracle_object)

from .load_promoter_base_GRN import (load_drosophila_promoter_base_GRN,
                                     load_rat_promoter_base_GRN,
                                     load_mouse_promoter_base_GRN,
                                     load_human_promoter_base_GRN,
                                     load_chicken_promoter_base_GRN,
                                     load_zebrafish_promoter_base_GRN,
                                     load_xenopus_tropicalis_promoter_base_GRN,
                                     load_arabidopsis_promoter_base_GRN,
                                     load_Scerevisiae_promoter_base_GRN,
                                     load_Celegans_promoter_base_GRN)

__all__ = ["load_TFinfo_df_mm9_mouse_atac_atlas",
           "load_mouse_scATAC_atlas_base_GRN"
           "load_Paul2015_data",
           "load_tutorial_links_object",
           "load_tutorial_oracle_object",
           "load_drosophila_promoter_base_GRN",
           "load_rat_promoter_base_GRN",
           "load_mouse_promoter_base_GRN",
           "load_human_promoter_base_GRN",
           "load_chicken_promoter_base_GRN",
           "load_zebrafish_promoter_base_GRN",
           "load_xenopus_tropicalis_promoter_base_GRN",
           "load_arabidopsis_promoter_base_GRN",
           "load_Scerevisiae_promoter_base_GRN",
           "load_Celegans_promoter_base_GRN"]
