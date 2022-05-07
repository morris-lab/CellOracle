# -*- coding: utf-8 -*-
import pandas as pd
from scanpy import AnnData

import celloracle as co

# Mouse scATAC base GRN loading function
def test_load_mouse_scATAC_atlas_base_GRN():
    df = co.data.load_mouse_scATAC_atlas_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)
