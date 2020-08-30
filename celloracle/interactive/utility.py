# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys

import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt

#import scanpy as sc
#import seaborn as sns
import h5py

ATTRS = ["delta_embedding", "delta_embedding_random"]

class Oracle_data_strage():
    def __init__(self):
        self.names = []
        self.path = None

    def create_hdf_path(self, path):
        self.path = path
        with h5py.File(self.path, mode='w') as f:
            pass

    def set_existing_hdf_path(self, path):
        self.path = path
        with h5py.File(self.path, mode='r') as f:
            self.names = [i[0] for i in list(f.items())]

    def save_data(self, oracle, name):
        with h5py.File(self.path, mode='r+') as f:
            if name in self.names:
                del f[name]
            for i, j in enumerate(ATTRS):
                f[f"{name}/{j}"] = getattr(oracle, j)
        self.names.append(name)
        self.names = list(set(self.names))

    def load_data(self, oracle, name):
        with h5py.File(self.path, mode='r') as f:
            for i, j in enumerate(ATTRS):
                val = f[f"{name}/{j}"][...]
                setattr(oracle, j, val)
        oracle.corrcoef_random = "dummy"
