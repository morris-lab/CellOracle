# -*- coding: utf-8 -*-
import pandas as pd
from scanpy import AnnData

import celloracle as co

# Mouse scATAC base GRN loading function
def test_load_mouse_scATAC_atlas_base_GRN():
    df = co.data.load_mouse_scATAC_atlas_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_mouse_scATAC_atlas_base_GRN_dl():
    df = co.data.load_mouse_scATAC_atlas_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

# Paul dataset anndata loading function
def test_load_Paul2015_data():
    adata = co.data.load_Paul2015_data(force_download=False)
    assert isinstance(adata, AnnData)

def test_load_Paul2015_data_dl():
    adata = co.data.load_Paul2015_data(force_download=True)
    assert isinstance(adata, AnnData)

# Tutorial Links object loading function
def test_load_tutorial_links_object():
    links = co.data.load_tutorial_links_object(force_download=False)
    assert isinstance(links, co.Links)

def test_load_tutorial_links_object_dl():
    links = co.data.load_tutorial_links_object(force_download=True)
    assert isinstance(links, co.Links)

# Tutorial oracle object loading function
def test_load_tutorial_oracle_object():
    oracle = co.data.load_tutorial_oracle_object(force_download=False)
    assert isinstance(oracle, co.Oracle)

def test_load_tutorial_oracle_object_dl():
    oracle = co.data.load_tutorial_oracle_object(force_download=True)
    assert isinstance(oracle, co.Oracle)
