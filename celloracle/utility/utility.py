# -*- coding: utf-8 -*-



import io
import logging
import os
import pickle
import subprocess
import sys

import pandas as pd
import numpy as np
from tqdm import tqdm_notebook as tqdm

from sklearn.preprocessing import StandardScaler


def save_as_pickled_object(obj, filepath):
    """
    Save any object using pickle.

    Args:
        obj (any python object): python object.

        filepath (str): file path.
    """
    max_bytes = 2**31 - 1
    bytes_out = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    n_bytes = sys.getsizeof(bytes_out)
    with open(filepath, 'wb') as f_out:
        for idx in range(0, n_bytes, max_bytes):
            f_out.write(bytes_out[idx:idx+max_bytes])


def load_pickled_object(filepath):
    """
    Load pickled object.

    Args:
        filepath (str): file path.

    Returns:
        python object: loaded object.

    """
    max_bytes = 2**31 - 1
    input_size = os.path.getsize(filepath)
    bytes_in = bytearray(0)
    with open(filepath, 'rb') as f_in:
        for _ in range(0, input_size, max_bytes):
            bytes_in += f_in.read(max_bytes)
    obj = pickle.loads(bytes_in)

    return obj



def intersect(list1, list2):
    """
    Intersect two list and get components that exists in both list.

    Args:
        list1 (list): input list.
        list2 (list): input list.

    Returns:
        list: intersected list.

    """
    inter_list = list(set(list1).intersection(list2))
    return(inter_list)


def inverse_dictionary(dictionary, verbose=True, return_value_as_numpy=False):
    """
    Make inversed dictionary.
    See examples below for detail.

    Args:
        dictionary (dict): python dictionary

        verbose (bool): Whether to show progress bar.

        return_value_as_numpy (bool): Whether to convert values into numpy array.

    Returns:
        dict: Python dictionary.

    Examples:
        >>> dic = {"a": [1, 2, 3], "b": [2, 3, 4]}
        >>> inverse_dictionary(dic)
        {1: ['a'], 2: ['a', 'b'], 3: ['a', 'b'], 4: ['b']}

        >>> dic = {"a": [1, 2, 3], "b": [2, 3, 4]}
        >>> inverse_dictionary(dic, return_value_as_numpy=True)
        {1: array(['a'], dtype='<U1'),
         2: array(['a', 'b'], dtype='<U1'),
         3: array(['a', 'b'], dtype='<U1'),
         4: array(['b'], dtype='<U1')}

    """
    def _get_key_list_that_contain_voi_in_its_values(value_of_interest, dictionary):
        keys_list = []
        for i in dictionary:
            values = dictionary[i]
            if value_of_interest in values:
                keys_list.append(i)
        return keys_list

    # get list of all TFs
    all_values = list(dictionary.values())
    all_values = np.unique(np.concatenate(all_values))

    keys_inversed_dict = all_values

    inversed_dictionary = {}

    if verbose:
        loop = tqdm(keys_inversed_dict)
    else:
        loop = keys_inversed_dict

    for key in loop:
        values_in_inversed_dict = _get_key_list_that_contain_voi_in_its_values(key, dictionary)
        values_in_inversed_dict = np.unique(values_in_inversed_dict)
        if not return_value_as_numpy:
            values_in_inversed_dict = list(values_in_inversed_dict)
        inversed_dictionary[key] = values_in_inversed_dict

    return  inversed_dictionary


def exec_process(commands, message=True, wait_finished=True, return_process=True):
    """
    Excute a command. This is a wrapper of "subprocess.Popen"

    Args:
        commands (str): command.
        message (bool): Whether to return a message or not.
        wait_finished (bool): Whether to wait for the process finished. If False, the function finish immediately.
        return_process (bool): Whether to return "process".


    """
    my_env = os.environ.copy()
    #print(my_env["PATH"])


    process = subprocess.Popen(commands, stdout=subprocess.PIPE, bufsize=-1,shell=True,
                               env=my_env)

    if message:
        with io.open(process.stdout.fileno(), closefd=False) as stream:
            [print(line.rstrip('\n')) for line in stream]
    if wait_finished:
        # プロセス終了まで待ち、結果を判定する
        process.wait()
        if process.returncode != 0:
            print('Build process aborts.')
            sys.exit(1)
    if return_process:
        return process

def standard(df):
    """
    Standerdize value.

    Args:
        df (padas.dataframe): dataframe.

    Returns:
        pandas.dataframe: data after standerdization.
        
    """
    dt = df.copy()
    #dt = dt.astype("float64")
    # function for data standerdization
    gene_names = dt.columns
    cell_names = dt.index

    scaler = StandardScaler(with_mean=False)
    scaler.fit(dt)
    dt = scaler.transform(dt)

    dt = pd.DataFrame(dt)
    dt.columns = gene_names
    dt.index = cell_names

    return(dt)
