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

from sklearn.preprocessing import StandardScaler

from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score

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
    Make inverse dictionary.
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
        wait_finished (bool): Whether or not to wait for the process to finish. If false, the process will be perfomed in background and the function will finish immediately
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
    Standardize value.

    Args:
        df (padas.dataframe): dataframe.

    Returns:
        pandas.dataframe: Data after standardization.

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

## Anndata processing functions

def update_adata(adata):
    # Update Anndata
    # Anndata generated with Scanpy 1.4 or less should be updated with this function
    # This function will be depricated in the future.

    try:
        lo = adata.uns['draw_graph']['params']['layout']
        if isinstance(lo, np.ndarray):
            lo = lo[0]
        adata.uns['draw_graph']['params']['layout'] = lo
    except:
        pass

def adata_to_color_dict(adata, cluster_use):
    """
    Extract color information from adata and returns as dictionary.

    Args:
        adata (anndata): anndata

        cluster_use (str): column name in anndata.obs

    Returns:
        dictionary: python dictionary, key is cluster name, value is clor name
    """
    color_dict = {}
    for i,j in enumerate(adata.obs[cluster_use].cat.categories):
        color_dict[j] = adata.uns[f"{cluster_use}_colors"][i]
    return color_dict



def transfer_all_colors_between_anndata(adata_ref, adata_que):
    """
    Extract all color information from reference anndata and transfer the color into query anndata.

    Args:
        adata_ref (anndata): reference anndata

        adata_que (anndata): query anndata

    """

    # Get all clusters that have colors in adata_ref
    keys = list(adata_ref.uns.keys())
    keys = [i.replace("_colors", "") for i in keys if "_colors" in i]

    keys = [i for i in keys if i in adata_que.obs.columns]

    for cluster_name in keys:
        transfer_color_between_anndata(adata_ref=adata_ref,
                                      adata_que=adata_que,
                                      cluster_name=cluster_name)

    print("Color meta data were transfered for \n")
    print(" ", keys)


def transfer_color_between_anndata(adata_ref, adata_que, cluster_name):
    """
    Extract color information from reference anndata and transfer the color into query anndata.

    Args:
        adata_ref (anndata): reference anndata

        adata_que (anndata): query anndata

        cluster_name (str): cluster name. This information should exist in the anndata.obs.

    """

    # Get color as a dictionary
    dic_ref = adata_to_color_dict(adata_ref, cluster_name)


    # Get color keys from que data
    adata_que.obs[cluster_name] = adata_que.obs[cluster_name].astype("category")
    color_keys = list(adata_que.obs[cluster_name].cat.categories)

    # Convert color list
    colors_after = []
    for i in color_keys:
        color = dic_ref[i]
        colors_after.append(color)

    colors_after = np.array(colors_after)

    # Update color
    adata_que.uns[cluster_name + "_colors"] = colors_after


def knn_data_transferer(adata_ref, adata_que,
                        n_neighbors=20, cluster_name=None, embedding_name=None, adata_true=None,
                        transfer_color=True, n_PCA=30, use_PCA_in_adata=False):
    """
    Extract categorical information from adata.obs or embedding information from adata.obsm and transfer these information into query anndata.
    In the calculation, KNN is used after PCA.

    Args:
        adata_ref (anndata): reference anndata

        adata_que (anndata): query anndata

        cluster_name (str or list of str): cluster name(s) to be transfered. If you want to transfer multiple data, you can set the cluster names as a list.

        embedding_name (str or list of str): embedding name(s) to be transfered. If you want to transfer multiple data, you can set the embedding names as a list.

        adata_true (str): This argument can be used for the validataion of this algorithm. If you have true answer, you can set it. If you set true answer, the function will return some metrics for benchmarking.

        transfer_color (bool): Whether or not to transfer color data in addition to cluster information.

        n_PCA (int): Number of PCs that will be used for the input of KNN algorithm.

    """
    # 1. Prepare data
    X_train = adata_ref.to_df()
    X_test = adata_que.to_df()

    genes = np.intersect1d(X_train.columns, X_test.columns)
    X_train = X_train[genes]
    X_test = X_test[genes]

    # 2. PCA
    if use_PCA_in_adata:
        X_train_PCA = adata_ref.obsm["X_pca"]
        X_test_PCA = adata_que.obsm["X_pca"]
    else:
        model_PCA = PCA(n_components=n_PCA)
        X_train_PCA = model_PCA.fit_transform(X_train)
        X_test_PCA = model_PCA.transform(X_test)

    # 3. Learning and prediction
    if cluster_name is not None:
        if isinstance(cluster_name, str):
            cluster_name = [cluster_name]
        if isinstance(cluster_name, list):
            for i in  cluster_name:
                model_kncl = KNeighborsClassifier(n_neighbors=n_neighbors)
                model_kncl.fit(X_train_PCA, adata_ref.obs[i])
                pred = model_kncl.predict(X_test_PCA)
                adata_que.obs[i] = pred
                adata_que.obs[i] = adata_que.obs[i].astype("category")

                if transfer_color:
                    transfer_color_between_anndata(adata_ref=adata_ref, adata_que=adata_que, cluster_name=i)

                if adata_true is not None:
                    true = adata_true.obs[i]
                    print(f"accuracy for {i}: ", np.mean(pred == true))
        else:
            raise ValueError("cluster name format error")

    if embedding_name is not None:
        if isinstance(embedding_name, str):
            embedding_name = [embedding_name]
        if isinstance(embedding_name, list):
            for i in  embedding_name:
                model_knreg = KNeighborsRegressor(n_neighbors=n_neighbors)
                model_knreg.fit(X_train_PCA, adata_ref.obsm[i])
                pred = model_knreg.predict(X_test_PCA)
                adata_que.obsm[i] = pred

                if adata_true is not None:
                    true = adata_true.obsm[i]

                    fig, ax = plt.subplots(1, 2)
                    for k in [0, 1]:
                        x = true[:, k]
                        y = pred[:, k]
                        ax[k].scatter(x, y)
                        ax[k].set_xlabel("true")
                        ax[k].set_ylabel("pred")
                        ax[k].set_title(f"r: {r2_score(x, y): .3g}")
                    plt.show()
        else:
            raise ValueError("embedding name format error")
