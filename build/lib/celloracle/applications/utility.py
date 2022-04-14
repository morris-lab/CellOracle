# -*- coding: utf-8 -*-



import io
import logging
import os
import sys

import pandas as pd
import numpy as np

import h5py

class Data_strage:
    """
    Custom class for quick save and load.
    Attributes can be easily and quickly saved and loaded with hdf5 file.
    Supported attribute types: [int, float, str, list(of int, float, str), np.ndarray (or int, float, str), pd.Dataframe]

    1. If an attribute is stored as an unsupperted type, it will be ignored.
    2. If an attribute name starts with "_", it will be ignored.
    """

    def __init__(self):
        self._exemptions_when_del_attrs = []

    def _set_hdf_path(self, path, key, create_if_not_exist=True):

        if not path.endswith(".hdf5"):
            raise ValueError("path should ends with '.hdf5'")

        self._path = path
        self._key = key
        self._names = []


        try:
            with h5py.File(self._path, mode='r') as f:
                f.visit(self._names.append)

        except:
            if create_if_not_exist:
                print("No hdf file found in the path. New hdf5 file was created.")
                with h5py.File(self._path, mode='w') as f:
                    pass

    def _save_attrs(self, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to save.
        """

        with h5py.File(self._path, mode='r+') as f:
            for i, j in enumerate(attributes):
                name_ = f"{place}/{j}"
                if name_ in self._names:
                    del f[name_]
                att = getattr(self, j)
                try:
                    f[name_] = att
                except:
                    print(j)
                    f[name_] = att.astype(h5py.string_dtype(encoding='utf-8'))
                self._names.append(name_)
        self._names = list(set(self._names))

    def _load_attrs(self, place, attributes):
        with h5py.File(self._path, mode='r') as f:
            for i, j in enumerate(attributes):
                val = f[f"{place}/{j}"][...]
                setattr(self, j, val)

    def _save_attrs_list(self, place, attributes):
        with h5py.File(self._path, mode='r+') as f:
            for i, j in enumerate(attributes):
                name_ = f"{place}/{j}"
                if name_ in self._names:
                    del f[name_]
                att = getattr(self, j)
                if len(att) == 0:
                    f[name_] = np.array(att)
                elif isinstance(att[0], str):
                    f[name_] = np.array(att).astype(h5py.string_dtype(encoding='utf-8'))
                    self._names.append(name_)
                elif type(att[0]) in [int, float]:
                    f[name_] = np.array(att)
                    self._names.append(name_)
        self._names = list(set(self._names))

    def _load_attrs_list(self, place, attributes):
        with h5py.File(self._path, mode='r') as f:
            for i, j in enumerate(attributes):
                val = f[f"{place}/{j}"][...]
                if len(val) > 0:
                    if isinstance(val[0], bytes):
                        setattr(self, j, [i.decode() for i in val])
                    else:
                        setattr(self, j, list(val))
                else:
                    setattr(self, j, [])

    def _save_attrs_misc(self, place, attributes):
        with h5py.File(self._path, mode='r+') as f:
            for i, j in enumerate(attributes):
                name_ = f"{place}/{j}"
                if name_ in self._names:
                    del f[name_]
                att = getattr(self, j)
                if isinstance(att, str):
                    f[name_] = np.array([att]).astype(h5py.string_dtype(encoding='utf-8'))
                    self._names.append(name_)
                elif type(att) in [int, float]:
                    f[name_] = np.array([att])
                    self._names.append(name_)
        self._names = list(set(self._names))

    def _load_attrs_misc(self, place, attributes):
        with h5py.File(self._path, mode='r') as f:
            for i, j in enumerate(attributes):
                val = f[f"{place}/{j}"][...]
                setattr(self, j, val[0])

    def _save_attrs_df(self, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to save.
        """
        for i, j in enumerate(attributes):
            key = f"{place}/{j}"
            att = getattr(self, j)
            att.to_hdf(self._path, key=key)

    def _load_attrs_df(self, place, attributes):
        """
        Args:
            place (str): Directry in hdf5 file.
            attributes (list of str): attributes to load.
        """
        for i, j in enumerate(attributes):
            key = f"{place}/{j}"
            att = pd.read_hdf(self._path, key=key)
            setattr(self, j, att)

    def _load_attrs_None(self, place, attributes):

        for i, j in enumerate(attributes):
            setattr(self, j, None)

    def _dump_hdf5(self, place=None):
        """
        Save all attribues in hdf5 file.

        Args:
            place (str): The key for saving data into hdf5 file.

        """
        if place is None:
            place = self._key


        attrs_df = []
        attrs_np = []
        attrs_list = []
        attrs_misc = []
        attrs_None = []

        for key, val in self.__dict__.items():
            if not key.startswith("_"):
                if type(val) in [pd.core.frame.DataFrame]:
                    attrs_df.append(key)
                elif type(val) in [np.ndarray]:
                    attrs_np.append(key)
                elif type(val) in [list]:
                    attrs_list.append(key)
                elif type(val) in [int, str, float]:
                    attrs_misc.append(key)
                elif val is None:
                    attrs_None.append(key)

        self.attrs_df = attrs_df
        self.attrs_np = attrs_np
        self.attrs_list = attrs_list
        self.attrs_misc = attrs_misc
        self.attrs_None = attrs_None

        self._save_attrs_df(place, attrs_df)
        self._save_attrs_list(place, attrs_list + ["attrs_df", "attrs_np", "attrs_list", "attrs_misc", "attrs_None"])
        self._save_attrs_misc(place, attrs_misc)
        self._save_attrs(place, attrs_np)


    def _load_hdf5(self, place=None, specify_attributes=None):
        """
        Load all attributes that was saved in hdf5 file.

        If you enter the list of attributes for specify_attributes, you can select the attribute to be loaded.
        This feature might be useful if you don't waste time and memory to load whold attributes.

        Args:
            place (str): The key in the hdf5 file that was used when save attributes.
            specify_attributes (list): Attribute names to be loaded.

        """

        if place is None:
            place = self._key


        self._load_attrs_list(place, ["attrs_df", "attrs_np", "attrs_list", "attrs_misc", "attrs_None"])

        if specify_attributes is not None:
            attrs_np = [i for i in self.attrs_np if i in specify_attributes]
            attrs_misc = [i for i in self.attrs_misc if i in specify_attributes]
            attrs_list = [i for i in self.attrs_list if i in specify_attributes]
            attrs_df = [i for i in self.attrs_df if i in specify_attributes]
            attrs_None = [i for i in self.attrs_None if i in specify_attributes]
        else:
            attrs_np = self.attrs_np
            attrs_misc = self.attrs_misc
            attrs_list = self.attrs_list
            attrs_df = self.attrs_df
            attrs_None = self.attrs_None

        self._load_attrs(place, attrs_np)
        self._load_attrs_misc(place, attrs_misc)
        self._load_attrs_list(place, attrs_list)
        self._load_attrs_df(place, attrs_df)
        self._load_attrs_None(place, attrs_None)

    def _del_attrs(self, exemptions=[]):
        """
        Delete attributes in the object.
        """

        del_exemptions = self._exemptions_when_del_attrs + exemptions
        for key, val in self.__dict__.items():
            if (not key.startswith("_")) & (not key in del_exemptions):
                setattr(self, key, None)
