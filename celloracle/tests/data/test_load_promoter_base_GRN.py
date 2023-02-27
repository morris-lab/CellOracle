# -*- coding: utf-8 -*-
import pandas as pd
from scanpy import AnnData

import celloracle as co

#
def test_load_drosophila_promoter_base_GRN_dl():
    df = co.data.load_drosophila_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_drosophila_promoter_base_GRN():
    df = co.data.load_drosophila_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_rat_promoter_base_GRN_dl():
    df = co.data.load_rat_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_rat_promoter_base_GRN():
    df = co.data.load_rat_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_mouse_promoter_base_GRN_dl():
    df = co.data.load_mouse_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_mouse_promoter_base_GRN():
    df = co.data.load_mouse_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_human_promoter_base_GRN_dl():
    df = co.data.load_human_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_human_promoter_base_GRN():
    df = co.data.load_human_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_chicken_promoter_base_GRN_dl():
    df = co.data.load_chicken_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_chicken_promoter_base_GRN():
    df = co.data.load_chicken_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_xenopus_tropicalis_promoter_base_GRN_dl():
    df = co.data.load_xenopus_tropicalis_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_xenopus_tropicalis_promoter_base_GRN():
    df = co.data.load_xenopus_tropicalis_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_xenopus_laevis_promoter_base_GRN_dl():
    df = co.data.load_xenopus_laevis_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_xenopus_laevis_promoter_base_GRN():
    df = co.data.load_xenopus_laevis_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_arabidopsis_promoter_base_GRN_dl():
    df = co.data.load_arabidopsis_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_arabidopsis_promoter_base_GRN():
    df = co.data.load_arabidopsis_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_Scerevisiae_promoter_base_GRN_dl():
    df = co.data.load_Scerevisiae_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_Scerevisiae_promoter_base_GRN():
    df = co.data.load_Scerevisiae_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_Celegans_promoter_base_GRN_dl():
    df = co.data.load_Celegans_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_Celegans_promoter_base_GRN():
    df = co.data.load_Celegans_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)

#
def test_load_Pig_promoter_base_GRN_dl():
    df = co.data.load_Pig_promoter_base_GRN(force_download=True)
    assert isinstance(df, pd.core.frame.DataFrame)

def test_load_Pig_promoter_base_GRN():
    df = co.data.load_Pig_promoter_base_GRN(force_download=False)
    assert isinstance(df, pd.core.frame.DataFrame)
