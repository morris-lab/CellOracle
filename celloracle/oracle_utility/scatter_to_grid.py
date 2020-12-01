# -*- coding: utf-8 -*-



import logging
import os
import sys

import pandas as pd
import numpy as np

from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.preprocessing import PolynomialFeatures


def scatter_value_to_grid_value(embedding, grid, value, method="knn", n_knn=30, n_poly=3):

    """
    Transfer value on 2D scatter into 2d Grid space.

    Args:

        embedding (np.array): shape = (n_cell, 2)
        grid (np.array): shape = (n_grid**2, 2)
        value (np.array): shape = (n_cell, )

    Returns:
        np.array : shape = (n_grid**2, 1)

    """

    x, y = embedding[:, 0], embedding[:, 1]
    x_new, y_new = grid[:, 0], grid[:, 1]

    if method == "poly":
        value_on_grid =  _polynomial_regression_old_ver(x, y, x_new, y_new, value, n_degree=n_poly)
    if method == "polynomial":
        value_on_grid =  _polynomial_regression_sklearn(x, y, x_new, y_new, value, n_degree=n_poly)
    elif method == "knn":
        value_on_grid = _knn_regression(x, y, x_new, y_new, value, n_knn=n_knn)
    elif method == "knn_class":
        value_on_grid = _knn_classification(x, y, x_new, y_new, value, n_knn=n_knn)


    return value_on_grid


def _polynomial_regression_sklearn(x, y, x_new, y_new, value, n_degree=3):

    # Make polynomial features
    data = np.stack([x, y], axis=1)
    data_new = np.stack([x_new, y_new], axis=1)

    pol = PolynomialFeatures(degree=n_degree, include_bias=False)
    data = pol.fit_transform(data)
    data_new = pol.transform(data_new)


    model = Ridge(random_state=123)
    model.fit(data, value)

    return model.predict(data_new)

def _polynomial_regression_old_ver(x, y, x_new, y_new, value, n_degree=3):

    def __conv(x, y, n_degree=3): # Make polynomial data for polynomial ridge regression
        dic = {}
        for d in range(1, n_degree + 1):
            for i, name in zip([x, y], ["x", "y"]):
                dic[name + str(d)] = i**d
        return pd.DataFrame(dic)

    data = __conv(x=x, y=y, n_degree=n_degree)

    model = Ridge()
    model.fit(data, value)

    data_new = __conv(x=x_new, y=y_new, n_degree=n_degree)

    return model.predict(data_new)

def _knn_regression(x, y, x_new, y_new, value, n_knn=30):

    data = np.stack([x, y], axis=1)

    model = KNeighborsRegressor(n_neighbors=n_knn)
    model.fit(data, value)

    data_new = np.stack([x_new, y_new], axis=1)

    return model.predict(data_new)

def _knn_classification(x, y, x_new, y_new, value, n_knn=30):

    data = np.stack([x, y], axis=1)

    model = KNeighborsClassifier(n_neighbors=n_knn)
    model.fit(data, value)

    data_new = np.stack([x_new, y_new], axis=1)

    return model.predict(data_new)
