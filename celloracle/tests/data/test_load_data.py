# -*- coding: utf-8 -*-
import pandas as pd

import celloracle as co

def test_base_GRN_loading_functions():
    df = co.data.load_mouse_scATAC_atlas_base_GRN()
    assert isinstance(df, pd.core.frame.DataFrame)
