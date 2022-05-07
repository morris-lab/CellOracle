# -*- coding: utf-8 -*-

import sys
import re
import warnings
import logging

from . import utility, network, network_analysis, go_analysis, data, data_conversion, oracle_utility
from .trajectory.oracle_core import Oracle
from .network import Net
from .network_analysis import Links
from .utility.load_hdf5 import load_hdf5
from .utility import check_python_requirements
from .network_analysis import test_R_libraries_installation
from .version import __version__
#from . import motif_analysis

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


# Make sure that DeprecationWarning within this package always gets printed
warnings.filterwarnings('always', category=DeprecationWarning,
                        module=r'^{0}\.'.format(re.escape(__name__)))

__copyright__    = 'Copyright (C) 2020 Kenji Kamimoto'
__license__      = 'Apache License Version 2.0'
__author__       = 'Kenji Kamimoto'
__author_email__ = 'kamimoto@wustl.edu'
__url__          = 'https://github.com/morris-lab/CellOracle'
#__version__ = "0.10.2"
__all__ = ["utility", "motif_analysis", "network", "network_analysis",
           "go_analysis", "data", "data_conversion",
           "Oracle",
           "Links",
           "Net",
           "load_hdf5",
           "test_R_libraries_installation",
           "check_python_requirements"
           ]

check_python_requirements(return_detail=False, print_warning=True)
