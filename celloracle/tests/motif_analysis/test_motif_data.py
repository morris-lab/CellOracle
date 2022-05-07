# -*- coding: utf-8 -*-
import pandas as pd

import celloracle as co

def test_load_motifs():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=False)
    assert isinstance(motifs, list)

def test_load_motifs_dl():
    motifs_name = 'CisBP_ver2_Arabidopsis_thaliana.pfm'
    motifs = co.motif_analysis.motif_data.load_motifs(motifs_name=motifs_name,
                                                      force_download=True)
    assert isinstance(motifs, list)
