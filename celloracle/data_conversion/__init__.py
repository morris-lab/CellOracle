# -*- coding: utf-8 -*-
"""
The :mod:`.data_conversion` module implements data conversion between different platform.

"""

from .process_seurat_object import seurat_object_to_anndata


__all__ = ["seurat_object_to_anndata"]
