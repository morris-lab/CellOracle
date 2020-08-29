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

class Oracle_data_strage():
    def __init__(self):
        self.strage = {}
        self.names = []
        self.id = 0
    def extract_data(self, oracle, name=None):
        if name is None:
            name = f"data_{self.id}"
            self.id += 1

        self.strage[name] = [#oracle.embedding.copy(),
                             oracle.delta_embedding.copy(),
                             oracle.delta_embedding_random.copy()
                             #oracle.perturb_condition.copy(),
                            ]
        self.names.append(name)

    def get_back_data(self, oracle, name):
        oracle.delta_embedding, oracle.delta_embedding_random = \
            self.strage[name]

        oracle.corrcoef_random = "dummy"

    def create_hdf_path(self, path):
        self.path = path
        with h5py.File(self.path, mode='w') as f:
            pass

    def set_existing_hdf_path(self, path):
        self.path = path
        with h5py.File(self.path, mode='r') as f:
            self.names = [i[0] for i in list(f.items())]

    def save_one_data(self, name):
        with h5py.File(self.path, mode='r+') as f:
            for i, j in enumerate(["delta_embedding", "delta_embedding_random"]):
                f[f"{name}/{j}"] = self.strage[name][i]


    def save_whole(self):
        with h5py.File(self.path, mode='r+') as f:
            for name in self.strage.keys():
                for i, j in enumerate(["delta_embedding", "delta_embedding_random"]):
                    f[f"{name}/{j}"] = self.strage[name][i]

    def load_one_data(self, name):
        with h5py.File(self.path, mode='r') as f:
            li = []
            for i, j in enumerate(["delta_embedding", "delta_embedding_random"]):
                data = f[f"{name}/{j}"][...]
                li.append(data)
            self.strage[name] = li

    def load_whole_data(self):
        with h5py.File(self.path, mode='r') as f:
            for name in self.names:
                li = []
                for i, j in enumerate(["delta_embedding", "delta_embedding_random"]):
                    data = f[f"{name}/{j}"]
                    li.append(data)
                self.strage[name] = li
