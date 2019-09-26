# -*- coding: utf-8 -*-
"""
The :mod:`.network` module implements GRN inference.

"""
from .net_core import Net
from .net_util import load_net_from_patquets, getDF_TGxTF, getDF_peakxTF, load_net


__all__ = ["Net",
           "load_net_from_patquets",
           "getDF_TGxTF",
           "getDF_peakxTF"
           "load_net"]
