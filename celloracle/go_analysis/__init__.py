# -*- coding: utf-8 -*-
"""
The :mod:`.go_analysis` module implements Gene Ontology analysis.
This module use goatools internally.

"""

from .goatools_wrapper import geneSymbol2ID, geneID2Symbol, get_GO


__all__ = ["geneSymbol2ID", "geneID2Symbol", "get_GO"]
