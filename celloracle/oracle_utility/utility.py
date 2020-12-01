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

#ATTRS = ["delta_embedding", "delta_embedding_random"]

class Oracle_data_strage():
    def __init__(self):
        pass

    def create_hdf_path(self, path):
        if not path.endswith(".hdf5"):
            raise ValueError("path should ends with '.hdf5'")
        self.path = path
        self.path_df = self.path.replace(".hdf5", "_df.hdf5")

        self.names = []
        with h5py.File(self.path, mode='w') as f:
            pass

    def set_hdf_path(self, path, create_if_not_exist=True):
        if not path.endswith(".hdf5"):
            raise ValueError("path should ends with '.hdf5'")
        self.path = path
        self.path_df = self.path.replace(".hdf5", "_df.hdf5")

        try:
            with h5py.File(self.path, mode='r') as f:
                self.names = []
                f.visit(self.names.append)

            unique_genes = np.unique([i.split("/")[0] for i in self.names])
            #print(f"hdf file with {len(unique_genes)} data was found.")
        except:
            print("No hdf file found in the path. New hdf5 file was created.")
            self.create_hdf_path(path=path)

    def save_data(self, oracle, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to save.
        """

        with h5py.File(self.path, mode='r+') as f:
            for i, j in enumerate(attributes):
                name_ = f"{place}/{j}"
                if name_ in self.names:
                    del f[name_]
                att = getattr(oracle, j)
                try:
                    f[name_] = att
                except:
                    f[name_] = att.astype(h5py.string_dtype(encoding='utf-8'))
                self.names.append(name_)
        self.names = list(set(self.names))

    def load_data(self, oracle, place, attributes):
        with h5py.File(self.path, mode='r') as f:
            for i, j in enumerate(attributes):
                val = f[f"{place}/{j}"][...]
                setattr(oracle, j, val)

    def save_dfs(self, oracle, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to save.
        """
        for i, j in enumerate(attributes):
            name_ = f"{place}/{j}"
            getattr(oracle, j).to_hdf(self.path_df, key=name_)

    def load_dfs(self, oracle, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to load.
        """
        for i, j in enumerate(attributes):
            name_ = f"{place}/{j}"
            setattr(oracle, j, pd.read_hdf(self.path_df, key=name_))
