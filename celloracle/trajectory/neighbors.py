import numpy as np
from numba import jit
from sklearn.neighbors import kneighbors_graph, NearestNeighbors
from scipy import sparse
import logging
from typing import *

# Mutual KNN functions


@jit(signature_or_function="Tuple((float64[:,:], int64[:,:], int64[:]))(int64[:,:], float64[:, :], int64[:], int64, int64, boolean)",
     nopython=True)
def balance_knn_loop(dsi: np.ndarray, dist: np.ndarray, lsi: np.ndarray, maxl: int, k: int, return_distance: bool) -> Tuple:
    """Fast and greedy algorythm to balance a K-NN graph so that no node is the NN to more than maxl other nodes

        Arguments
        ---------
        dsi : np.ndarray  (samples, K)
            distance sorted indexes (as returned by sklearn NN)
        dist : np.ndarray  (samples, K)
            the actual distance corresponding to the sorted indexes
        lsi : np.ndarray (samples,)
            degree of connectivity (l) sorted indexes
        maxl : int
            max degree of connectivity (from others to self) allowed
        k : int
            number of neighbours in the final graph
        return_distance : bool
            whether to return distance

        Returns
        -------
        dsi_new : np.ndarray (samples, k+1)
            indexes of the NN, first column is the sample itself
        dist_new : np.ndarray (samples, k+1)
            distances to the NN
        l: np.ndarray (samples)
            l[i] is the number of connections from other samples to the sample i

    """
    assert dsi.shape[1] >= k, "sight needs to be bigger than k"
    # numba signature "Tuple((int64[:,:], float32[:, :], int64[:]))(int64[:,:], int64[:], int64, int64, bool)"
    dsi_new = -1 * np.ones((dsi.shape[0], k + 1), np.int64)  # maybe d.shape[0]
    l = np.zeros(dsi.shape[0], np.int64)
    if return_distance:
        dist_new = np.zeros(dsi_new.shape, np.float64)
    for i in range(dsi.shape[0]):  # For every node
        el = lsi[i]
        p = 0
        j = 0
        for j in range(dsi.shape[1]):  # For every other node it is connected (sight)
            if p >= k:
                break
            m = dsi[el, j]
            if el == m:
                dsi_new[el, 0] = el
                continue
            if l[m] >= maxl:
                continue
            dsi_new[el, p + 1] = m
            l[m] = l[m] + 1
            if return_distance:
                dist_new[el, p + 1] = dist[el, j]
            p += 1
        if (j == dsi.shape[1] - 1) and (p < k):
            while p < k:
                dsi_new[el, p + 1] = el
                dist_new[el, p + 1] = dist[el, 0]
                p += 1
    if not return_distance:
        dist_new = np.ones_like(dsi_new, np.float64)
    return dist_new, dsi_new, l


@jit(signature_or_function="Tuple((float64[:,:], int64[:,:], int64[:]))(int64[:,:], float64[:, :], int64[:], int64[:], int64, int64, boolean)",
     nopython=True)
def balance_knn_loop_constrained(dsi: np.ndarray, dist: np.ndarray, lsi: np.ndarray, groups: np.ndarray, maxl: int, k: int, return_distance: bool) -> Tuple:
    """Fast and greedy algorythm to balance a K-NN graph so that no node is the NN to more than maxl other nodes

        Arguments
        ---------
        dsi : np.ndarray  (samples, K)
            distance sorted indexes (as returned by sklearn NN)
        dist : np.ndarray  (samples, K)
            the actual distance corresponding to the sorted indexes
        lsi : np.ndarray (samples,)
            degree of connectivity (l) sorted indexes
        groups: np.ndarray (samples,)
            labels of the samples that constrain the connectivity
        maxl : int
            max degree of connectivity (from others to self) allowed
        k : int
            number of neighbours in the final graph
        return_distance : bool
            whether to return distance

        Returns
        -------
        dsi_new : np.ndarray (samples, k+1)
            indexes of the NN, first column is the sample itself
        dist_new : np.ndarray (samples, k+1)
            distances to the NN
        l: np.ndarray (samples)
            l[i] is the number of connections from other samples to the sample i

    """
    assert dsi.shape[1] >= k, "sight needs to be bigger than k"
    # numba signature "Tuple((int64[:,:], float32[:, :], int64[:]))(int64[:,:], int64[:], int64, int64, bool)"
    dsi_new = -1 * np.ones((dsi.shape[0], k + 1), np.int64)  # maybe d.shape[0]
    l = np.zeros(dsi.shape[0], np.int64)
    if return_distance:
        dist_new = np.zeros(dsi_new.shape, np.float64)
    for i in range(dsi.shape[0]):  # For every node
        el = lsi[i]
        p = 0
        j = 0
        for j in range(dsi.shape[1]):  # For every other node it is connected (sight)
            if p >= k:  # if k-nn were found
                break
            m = dsi[el, j]
            if el == m:
                dsi_new[el, 0] = el
                continue
            if groups[el] != groups[m]:  # This is the constraint!
                continue
            if l[m] >= maxl:
                continue
            dsi_new[el, p + 1] = m
            l[m] = l[m] + 1
            if return_distance:
                dist_new[el, p + 1] = dist[el, j]
            p += 1
        if (j == dsi.shape[1] - 1) and (p < k):
            while p < k:
                dsi_new[el, p + 1] = el
                dist_new[el, p + 1] = dist[el, 0]
                p += 1
    if not return_distance:
        dist_new = np.ones_like(dsi_new, np.float64)
    return dist_new, dsi_new, l


def knn_balance(dsi: np.ndarray, dist: np.ndarray=None, maxl: int=200, k: int=60, constraint: np.ndarray=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Balance a K-NN graph so that no node is the NN to more than maxl other nodes

        Arguments
        ---------
        dsi : np.ndarray  (samples, K)
            distance sorted indexes (as returned by sklearn NN)
        dist : np.ndarray  (samples, K)
            the actual distance corresponding to the sorted indexes
        maxl : int
            max degree of connectivity allowed
        k : int
            number of neighbours in the final graph
        constraint: np.ndarray (samples,)
            labels of the samples that constrain the connectivity

        Returns
        -------
        dist_new : np.ndarray (samples, k+1)
            distances to the NN
        dsi_new : np.ndarray (samples, k+1)
            indexes of the NN, first column is the sample itself
        l: np.ndarray (samples)
            l[i] is the number of connections from other samples to the sample i


    """

    dsi = dsi.astype("int64").copy(order="C")
    l = np.bincount(dsi.flat[:], minlength=dsi.shape[0])
    lsi = np.argsort(l, kind="mergesort")[::-1]

    if dist is None:
        dist = np.ones(dsi.shape, dtype="float64").copy(order="C")
        dist[:, 0] = 0

        # Change data type and contingency for numba calculation. Added by Kenji 2021/1/21
        dsi = dsi.astype("int64").copy(order="C")
        dist = dist.astype("float64").copy(order="C")
        lsi = lsi.astype("int64").copy(order="C")

        if constraint is not None:
            return balance_knn_loop_constrained(dsi, dist, lsi, constraint.astype("int64"), maxl, k, return_distance=False)
        else:
            return balance_knn_loop(dsi, dist, lsi, maxl, k, return_distance=False)
    else:
        # Change data type and contingency for numba calculation. Added by Kenji 2021/1/21
        dsi = dsi.astype("int64").copy(order="C")
        dist = dist.astype("float64").copy(order="C")
        lsi = lsi.astype("int64").copy(order="C")

        if constraint is not None:
            return balance_knn_loop_constrained(dsi, dist, lsi, constraint.astype("int64"), maxl, k, return_distance=True)
        else:
            return balance_knn_loop(dsi, dist, lsi, maxl, k, return_distance=True)


class BalancedKNN:
    """Greedy algorythm to balance a K-nearest neighbour graph

    It has an API similar to scikit-learn

    Parameters
    ----------
    k : int  (default=50)
        the number of neighbours in the final graph
    sight_k : int  (default=100)
        the number of neighbours in the initialization graph
        It correspondent to the farthest neighbour that a sample is allowed to connect to
        when no closest neighbours are allowed. If sight_k is reached then the matrix is filled
        with the sample itself
    maxl : int  (default=200)
         max degree of connectivity allowed. Avoids the presence of hubs in the graph, it is the
         maximum number of neighbours that are allowed to contact a node before the node is blocked
    constraint: np.ndarray (default=None)
        a numpy array defining the dirrent groups within wich connectivity is allowed
        if "clusters" it uses the clusters as in self.clusters_ix
    mode : str (default="connectivity")
        decide wich kind of utput "distance" or "connectivity"
    n_jobs : int  (default=4)
        parallelization of the standard KNN search preformed at initialization
    """
    def __init__(self, k: int=50, sight_k: int=100, maxl: int=200, constraint: np.ndarray=None, mode: str="distance", metric: str="euclidean", n_jobs: int=4) -> None:
        self.k = k
        self.sight_k = sight_k
        self.maxl = maxl
        self.mode = mode
        self.metric = metric
        self.n_jobs = n_jobs
        self.dist_new = self.dsi_new = self.l = None  # type: np.ndarray
        self.bknn = None  # type: sparse.csr_matrix
        self.constraint = constraint

    @property
    def n_samples(self) -> int:
        return self.data.shape[0]

    def fit(self, data: np.ndarray, sight_k: int=None) -> Any:
        """Fits the model

        data: np.ndarray (samples, features)
            np
        sight_k: int
            the farthest point that a node is allowed to connect to when its closest neighbours are not allowed
        """
        self.data = data
        self.fitdata = data
        if sight_k is not None:
            self.sight_k = sight_k
        logging.debug(f"First search the {self.sight_k} nearest neighbours for {self.n_samples}")
        if self.metric == "correlation":
            self.nn = NearestNeighbors(n_neighbors=self.sight_k + 1, metric=self.metric, n_jobs=self.n_jobs, algorithm="brute")
        else:
            self.nn = NearestNeighbors(n_neighbors=self.sight_k + 1, metric=self.metric, n_jobs=self.n_jobs, leaf_size=30)
        self.nn.fit(self.fitdata)
        return self

    def kneighbors(self, X: np.ndarray=None, maxl: int=None, mode: str="distance") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Finds the K-neighbors of a point.

            Returns indices of and distances to the neighbors of each point.

            Parameters
            ----------
            X : array-like, shape (n_query, n_features),
                The query point or points.
                If not provided, neighbors of each indexed point are returned.
                In this case, the query point is not considered its own neighbor.

            maxl: int
                max degree of connectivity allowed

            mode : "distance" or "connectivity"
                Decides the kind of output

            Returns
            -------
            dist_new : np.ndarray (samples, k+1)
                distances to the NN
            dsi_new : np.ndarray (samples, k+1)
                indexes of the NN, first column is the sample itself
            l: np.ndarray (samples)
                l[i] is the number of connections from other samples to the sample i

            NOTE:
            First column (0) correspond to the sample itself, the nearest nenigbour is at the second column (1)

        """
        if X is not None:
            self.data = X
        if maxl is not None:
            self.maxl = maxl

        self.dist, self.dsi = self.nn.kneighbors(self.data, return_distance=True)
        logging.debug(f"Using the initialization network to find a {self.k}-NN graph with maximum connectivity of {self.maxl}")
        self.dist_new, self.dsi_new, self.l = knn_balance(self.dsi, self.dist, maxl=self.maxl, k=self.k, constraint=self.constraint)

        if mode == "connectivity":
            self.dist = np.ones_like(self.dsi)
            self.dist[:, 0] = 0
        return self.dist_new, self.dsi_new, self.l

    def kneighbors_graph(self, X: np.ndarray=None, maxl: int=None, mode: str="distance") -> sparse.csr_matrix:
        """Retrun the K-neighbors graph as a sparse csr matrix

            Parameters
            ----------
            X : array-like, shape (n_query, n_features),
                The query point or points.
                If not provided, neighbors of each indexed point are returned.
                In this case, the query point is not considered its own neighbor.

            maxl: int
                max degree of connectivity allowed

            mode : "distance" or "connectivity"
                Decides the kind of output

            Returns
            -------
            neighbor_graph : scipy.sparse.csr_matrix
                The values are either distances or connectivity dependig of the mode parameter

            NOTE: The diagonal will be zero even though the value 0 is actually stored

        """
        dist_new, dsi_new, l = self.kneighbors(X=X, maxl=maxl, mode=mode)
        logging.debug("Returning sparse matrix")
        self.bknn = sparse.csr_matrix((np.ravel(dist_new),
                                       np.ravel(dsi_new),
                                       np.arange(0, dist_new.shape[0] * dist_new.shape[1] + 1, dist_new.shape[1])),
                                      (self.n_samples,
                                       self.n_samples))
        return self.bknn

    def smooth_data(self, data_to_smooth: np.ndarray, X: np.ndarray=None, maxl: int=None,
                    mutual: bool=False, only_increase: bool=True) -> np.ndarray:
        """Use the wights learned from knn to smooth any data matrix

        Arguments
        ---------
        data_to_smooth: (features, samples) !! NOTE !! this is different from the input (for speed issues)
            if the data is provided (samples, features), this will be detected and
            the correct operation performed at cost of some effciency
            In the case where samples == samples then the shape (features, samples) will be assumed

        """
        if self.bknn is None:
            assert (X is None) and (maxl is None), "graph was already fit with different parameters"
            self.kneighbors_graph(X=X, maxl=maxl, mode=self.mode)
        if mutual:
            connectivity = make_mutual(self.bknn > 0)
        else:
            connectivity = self.bknn.T > 0
        connectivity = connectivity.tolil()
        connectivity.setdiag(1)
        w = connectivity_to_weights(connectivity).T
        assert np.allclose(w.sum(0), 1), "weight matrix need to sum to one over the columns"
        if data_to_smooth.shape[1] == w.shape[0]:
            result = sparse.csr_matrix.dot(data_to_smooth, w)
        elif data_to_smooth.shape[0] == w.shape[0]:
            result = sparse.csr_matrix.dot(data_to_smooth.T, w).T
        else:
            raise ValueError(f"Incorrect size of matrix, none of the axis correspond to the one of graph. {w.shape}")

        if only_increase:
            return np.maximum(result, data_to_smooth)
        else:
            return result


# Mutual KNN version


def knn_distance_matrix(data: np.ndarray, metric: str=None, k: int=40, mode: str='connectivity', n_jobs: int=4) -> sparse.csr_matrix:
    """Calculate a nearest neighbour distance matrix

    Notice that k is meant as the actual number of neighbors NOT INCLUDING itself
    To achieve that we call kneighbors_graph with X = None
    """
    if metric == "correlation":
        nn = NearestNeighbors(n_neighbors=k, metric="correlation", algorithm="brute", n_jobs=n_jobs)
        nn.fit(data)
        return nn.kneighbors_graph(X=None, mode=mode)
    else:
        nn = NearestNeighbors(n_neighbors=k, n_jobs=n_jobs, )
        nn.fit(data)
        return nn.kneighbors_graph(X=None, mode=mode)


def make_mutual(knn: sparse.csr.csr_matrix) -> sparse.coo_matrix:
    """Removes edges between neighbours that are not mutual
    """
    return knn.minimum(knn.T)


def connectivity_to_weights(mknn: sparse.csr.csr_matrix, axis: int=1) -> sparse.lil_matrix:
    """Convert a binary connectivity matrix to weights ready to be multiplied to smooth a data matrix
    """
    if type(mknn) is not sparse.csr.csr_matrix:
        mknn = mknn.tocsr()
    return mknn.multiply(1. / sparse.csr_matrix.sum(mknn, axis=axis))


def min_n(row_data: np.ndarray, row_indices: np.ndarray, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Find the smallest entry and smallest indices of a row
    """
    i = row_data.argsort()[:n]
    # i = row_data.argpartition(-n)[-n:]
    top_values = row_data[i]
    top_indices = row_indices[i]  # do the sparse indices matter?
    return top_values, top_indices


def take_top(matrix: sparse.spmatrix, n: int) -> sparse.lil_matrix:
    """Filter the top nearest neighbours from a sprse distance matrix
    """
    arr_ll = matrix.tolil(copy=True)
    for i in range(arr_ll.shape[0]):
        d, r = min_n(np.array(arr_ll.data[i]), np.array(arr_ll.rows[i]), n)
        arr_ll.data[i] = d.tolist()
        arr_ll.rows[i] = r.tolist()
    return arr_ll


# Common functions

def convolve_by_sparse_weights(data: np.ndarray, w: sparse.csr_matrix) -> np.ndarray:
    """Use the wights learned from knn to convolve any data matrix

    NOTE: A improved implementation could detect wich one is sparse and wich kind of sparse and perform faster computation
    """
    w_ = w.T
    assert np.allclose(w_.sum(0), 1), "weight matrix need to sum to one over the columns"
    return sparse.csr_matrix.dot(data, w_)


def knn_smooth_weights(matrix: np.ndarray, metric: str="euclidean", k_search: int=20,
                       k_mutual: int=10, n_jobs: int=10) -> Tuple[sparse.spmatrix, sparse.csr_matrix]:
    """Find the weights to smooth the dataset using efficient sparse matrix operations

    Arguments:
        matrix: (genes, cells)
            expression matrix
        metric
        k_search : int
            the first k nearest neighbour search number of neighbours
        k_mutual : int
            the number of mutual neighbours to select
        n_jobs
        return_knn

    Retruns
        weights (, knn)
    """
    assert k_search >= k_mutual, "k_search needs to be bigger than k_mutual"
    knn = knn_distance_matrix(matrix.T, metric=metric, k=k_search, mode="distance", n_jobs=n_jobs)
    mknn = make_mutual(knn)
    top_mknn = take_top(mknn, k_mutual)
    top_mknn.setdiag(1)
    connectivity = top_mknn > 0
    w = connectivity_to_weights(connectivity)
    return w, knn
