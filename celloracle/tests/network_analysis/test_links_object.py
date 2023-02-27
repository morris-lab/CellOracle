# -*- coding: utf-8 -*-


import celloracle as co
import pandas as pd
import igraph

def test_get_network_score():
    """
    Several version of the igraph package has incompatibility with Celloracle.
    The function will check whether the igraph package is working.
    """
    print("igraph version: ", igraph.__version__)
    links = co.data.load_tutorial_links_object()
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
    links.get_network_score()
    assert isinstance(links.merged_score, pd.core.frame.DataFrame)
