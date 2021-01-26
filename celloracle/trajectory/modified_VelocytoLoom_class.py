import logging
import warnings
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Any, List, Union, Tuple
import pandas as pd
import scipy.stats
from numba import jit
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm as normal
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors

from velocyto.diffusion import Diffusion
from velocyto.estimation import (colDeltaCor, colDeltaCorLog10,
                                 colDeltaCorLog10partial, colDeltaCorpartial,
                                 colDeltaCorSqrt, colDeltaCorSqrtpartial)
#from velocyto.neighbors import (BalancedKNN, connectivity_to_weights,
#                                convolve_by_sparse_weights,
#                                knn_distance_matrix)
from .neighbors import (BalancedKNN, connectivity_to_weights,
                                convolve_by_sparse_weights,
                                knn_distance_matrix)
from velocyto.serialization import dump_hdf5, load_hdf5


#from tqdm.notebook import tqdm
from .oracle_utility import _adata_to_matrix, _adata_to_df, _get_clustercolor_from_anndata


class modified_VelocytoLoom():

    def __init__(self):

        pass



    def score_cv_vs_mean(self, N: int=3000, min_expr_cells: int=2, max_expr_avg: float=20, min_expr_avg: int=0, svr_gamma: float=None,
                         winsorize: bool=False, winsor_perc: Tuple[float, float]=(1, 99.5), sort_inverse: bool=False, plot: bool=False) -> np.ndarray:
        from sklearn.svm import SVR
        """Rank genes on the basis of a CV vs mean fit, it uses a nonparametric fit (Support Vector Regression)

        Arguments
        ---------
        N: int
            the number to select
        min_expr_cells: int, (default=2)
            minimum number of cells that express that gene for it to be considered in the fit
        min_expr_avg: int, (default=0)
            The minimum average accepted before discarding from the the gene as not expressed
        max_expr_avg: float, (default=20)
            The maximum average accepted before discarding from the the gene as house-keeping/outlier
        svr_gamma: float
            the gamma hyper-parameter of the SVR
        winsorize: bool
            Wether to winsorize the data for the cv vs mean model
        winsor_perc: tuple, default=(1, 99.5)
            the up and lower bound of the winsorization
        sort_inverse: bool, (default=False)
            if True it sorts genes from less noisy to more noisy (to use for size estimation not for feature selection)
        which: bool, (default="S")
            it performs the same cv_vs mean procedure on spliced "S" or unspliced "U" count
            "both" is NOT supported here because most often S the two procedure would have different parameters
            (notice that default parameters are good heuristics only for S)
        plot: bool, default=False
            whether to show a plot

        Returns
        -------
        Nothing but it creates the attributes
        cv_mean_score: np.ndarray
            How much the observed CV is higher than the one predicted by a noise model fit to the data
        cv_mean_selected: np.ndarray bool
            on the basis of the N parameter

        Note: genes excluded from the fit will have in the output the same score as the lowest scoring gene in the dataset.

        To perform the filtering use the method `filter_genes`
        """

        X = _adata_to_matrix(self.adata, "raw_count")
        if winsorize:
            if min_expr_cells <= ((100 - winsor_perc[1]) * X.shape[1] * 0.01):
                min_expr_cells = int(np.ceil((100 - winsor_perc[1]) * X.shape[0] * 0.01)) + 2
                logging.debug(f"min_expr_cells is too low for winsorization with upper_perc ={winsor_perc[1]}, upgrading to min_expr_cells ={min_expr_cells}")

        detected_bool = ((X > 0).sum(1) > min_expr_cells) & (X.mean(1) < max_expr_avg) & (X.mean(1) > min_expr_avg)
        Sf = X[detected_bool, :]
        if winsorize:
            down, up = np.percentile(Sf, winsor_perc, 1)
            Sfw = np.clip(Sf, down[:, None], up[:, None])
            mu = Sfw.mean(1)
            sigma = Sfw.std(1, ddof=1)
        else:
            mu = Sf.mean(1)
            sigma = Sf.std(1, ddof=1)

        cv = sigma / mu
        log_m = np.log2(mu)
        log_cv = np.log2(cv)

        if svr_gamma is None:
            svr_gamma = 150. / len(mu)
            logging.debug(f"svr_gamma set to {svr_gamma}")
        # Fit the Support Vector Regression
        clf = SVR(gamma=svr_gamma)
        clf.fit(log_m[:, None], log_cv)
        fitted_fun = clf.predict
        ff = fitted_fun(log_m[:, None])
        score = log_cv - ff
        if sort_inverse:
            score = - score
        nth_score = np.sort(score)[::-1][N]
        if plot:
            scatter_viz(log_m[score > nth_score], log_cv[score > nth_score], s=3, alpha=0.4, c="tab:red")
            scatter_viz(log_m[score <= nth_score], log_cv[score <= nth_score], s=3, alpha=0.4, c="tab:blue")
            mu_linspace = np.linspace(np.min(log_m), np.max(log_m))
            plt.plot(mu_linspace, fitted_fun(mu_linspace[:, None]), c="k")
            plt.xlabel("log2 mean S")
            plt.ylabel("log2 CV S")
        self.cv_mean_score = np.zeros(detected_bool.shape)
        self.cv_mean_score[~detected_bool] = np.min(score) - 1e-16
        self.cv_mean_score[detected_bool] = score
        self.cv_mean_selected = self.cv_mean_score >= nth_score
        self.cv_mean_selected_genes = self.adata.var.index[self.cv_mean_selected].values



    def perform_PCA(self, n_components: int=None, div_by_std: bool=False) -> None:
        """Perform PCA (cells as samples)

        Arguments
        ---------
        which: str, default="S_norm"
            The name of the attribute to use for the calculation (e.g. S_norm or Sx_norm)
        n_components: int, default=None
            Number of components to keep. If None all the components will be kept.
        div_by_std: bool, default=False
            Wether to divide by standard deviation

        Returns
        -------
        Returns nothing but it creates the attributes:
        pca: np.ndarray
            a numpy array of shape (cells, npcs)

        """
        X = _adata_to_matrix(self.adata, "normalized_count")

        self.pca = PCA(n_components=n_components)
        if div_by_std:
            self.pcs = self.pca.fit_transform(X.T / X.std(0))
        else:
            self.pcs = self.pca.fit_transform(X.T)



    def plot_pca(self, dim: List[int]=[0, 1, 2], elev: float=60, azim: float=-140) -> None:
        """Plot 3d PCA
        """

        # update color information
        col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=self.cluster_column_name,
                                                  return_as="dict")
        self.colorandum = np.array([col_dict[i] for i in self.adata.obs[self.cluster_column_name]])

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.pcs[:, dim[0]],
                   self.pcs[:, dim[1]],
                   self.pcs[:, dim[2]],
                   c=self.colorandum)
        ax.view_init(elev=elev, azim=azim)



    def knn_imputation(self, k: int=None, metric: str="euclidean", diag: float=1,
                       n_pca_dims: int=None, maximum: bool=False,
                       balanced: bool=False, b_sight: int=None, b_maxl: int=None,
                       group_constraint: Union[str, np.ndarray]=None, n_jobs: int=8) -> None:
        """Performs k-nn smoothing of the data matrix

        Arguments
        ---------
        k: int
            number of neighbors. If None the default it is chosen to be `0.025 * Ncells`
        metric: str
            "euclidean" or "correlation"
        diag: int, default=1
            before smoothing this value is substituted in the diagonal of the knn contiguity matrix
            Resulting in a reduction of the smoothing effect.
            E.g. if diag=8 and k=10 value of Si = (8 * S_i + sum(S_n, with n in 5nn of i)) / (8+5)
        maximum: bool, default=False
            If True the maximum value of the smoothing and the original matrix entry is taken.
        n_pca_dims: int, default=None
            number of pca to use for the knn distance metric. If None all pcs will be used. (used only if pca_space == True)
        balanced: bool
            whether to use BalancedKNN version
        b_sight: int
            the sight parameter of BalancedKNN (used only if balanced == True)
        b_maxl: int
            the maxl parameter of BalancedKNN (used only if balanced == True)

        n_jobs: int, default 8
            number of parallel jobs in knn calculation

        Returns
        -------
        Nothing but it creates the attributes:
        knn: scipy.sparse.csr_matrix
            knn contiguity matrix
        knn_smoothing_w: scipy.sparse.lil_matrix
            the weights used for the smoothing
        Sx: np.ndarray
            smoothed spliced
        Ux: np.ndarray
            smoothed unspliced

        """
        X = _adata_to_matrix(self.adata, "normalized_count")


        N = self.adata.shape[0] # cell number

        if k is None:
            k = int(N * 0.025)
        if b_sight is None and balanced:
            b_sight = int(k * 8)
        if b_maxl is None and balanced:
            b_maxl = int(k * 4)

        space = self.pcs[:, :n_pca_dims]

        if balanced:
            bknn = BalancedKNN(k=k, sight_k=b_sight, maxl=b_maxl,
                               metric=metric, mode="distance", n_jobs=n_jobs)
            bknn.fit(space)
            self.knn = bknn.kneighbors_graph(mode="distance")
        else:

            self.knn = knn_distance_matrix(space, metric=metric, k=k,
                                           mode="distance", n_jobs=n_jobs)
        connectivity = (self.knn > 0).astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # SparseEfficiencyWarning: Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.
            connectivity.setdiag(diag)
        self.knn_smoothing_w = connectivity_to_weights(connectivity)

        ###
        Xx = convolve_by_sparse_weights(X, self.knn_smoothing_w)
        self.adata.layers["imputed_count"] = Xx.transpose().copy()


    def estimate_transition_prob(self,
                                 n_neighbors: int=None,
                                 knn_random: bool=True, sampled_fraction: float=0.3,
                                 sampling_probs: Tuple[float, float]=(0.5, 0.1),
                                 n_jobs: int=4, threads: int=None, calculate_randomized: bool=True,
                                 random_seed: int=15071990) -> None:
        """Use correlation to estimate transition probabilities for every cells to its embedding neighborhood

        Arguments
        ---------
        embed: str, default="ts"
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
        transform: str, default="sqrt"
            The transformation that is applies on the high dimensional space.
            If None the raw data will be used

        n_sight: int, default=None (also n_neighbors)
            The number of neighbors to take into account when performing the projection
        knn_random: bool, default=True
            whether to random sample the neighborhoods to speedup calculation
        sampling_probs: Tuple, default=(0.5, 1)
        max_dist_embed: float, default=None
            CURRENTLY NOT USED
            The maximum distance allowed
            If None it will be set to 0.25 * average_distance_two_points_taken_at_random
        n_jobs: int, default=4
            number of jobs to calculate knn
            this only applies to the knn search, for the more time consuming correlation computation see threads
        threads: int, default=None
            The threads will be used for the actual correlation computation by default half of the total.
        calculate_randomized: bool, default=True
            Calculate the transition probabilities with randomized residuals.
            This can be plotted downstream as a negative control and can be used to adjust the visualization scale of the velocity field.
        random_seed: int, default=15071990
            Random seed to make knn_random mode reproducible

        Returns
        -------
        """

        numba_random_seed(random_seed)

        X = _adata_to_matrix(self.adata, "imputed_count")  # [:, :ndims]
        delta_X = _adata_to_matrix(self.adata, "delta_X")
        embedding = self.adata.obsm[self.embedding_name]
        self.embedding = embedding

        if n_neighbors is None:
            n_neighbors = int(self.adata.shape[0] / 5)


        if knn_random:
            np.random.seed(random_seed)
            self.corr_calc = "knn_random"


            if calculate_randomized:
                delta_X_rndm = np.copy(delta_X)
                permute_rows_nsign(delta_X_rndm)

            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
            nn.fit(embedding)  # NOTE should support knn in high dimensions
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity")

            # Pick random neighbours and prune the rest
            neigh_ixs = self.embedding_knn.indices.reshape((-1, n_neighbors + 1))
            p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
            p = p / p.sum()

            # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
            # resulting self.corrcoeff with different number of nonzero entry per row.
            # Not updated yet not to break previous analyses
            # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
            # I change it here since I am doing some breaking changes
            sampling_ixs = np.stack([np.random.choice(neigh_ixs.shape[1],
                                                      size=(int(sampled_fraction * (n_neighbors + 1)),),
                                                      replace=False,
                                                      p=p) for i in range(neigh_ixs.shape[0])], 0)
            self.sampling_ixs = sampling_ixs
            neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
            nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
            self.embedding_knn = sparse.csr_matrix((np.ones(nonzero),
                                                    neigh_ixs.ravel(),
                                                    np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                                   shape=(neigh_ixs.shape[0],
                                                          neigh_ixs.shape[0]))

            logging.debug(f"Correlation Calculation '{self.corr_calc}'")

            ###
            ###
            self.corrcoef = colDeltaCorpartial(X, delta_X, neigh_ixs, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                self.corrcoef_random = colDeltaCorpartial(X, delta_X_rndm, neigh_ixs, threads=threads)
            ######

            if np.any(np.isnan(self.corrcoef)):
                self.corrcoef[np.isnan(self.corrcoef)] = 1
                logging.debug("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
                #logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)
                if np.any(np.isnan(self.corrcoef_random)):
                    self.corrcoef_random[np.isnan(self.corrcoef_random)] = 1
                    #logging.warning("Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
                    logging.debug("Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
            logging.debug(f"Done Correlation Calculation")
        else:
            self.corr_calc = "full"

            if calculate_randomized:
                delta_X_rndm = np.copy(delta_X)
                permute_rows_nsign(delta_X_rndm)

            logging.debug("Calculate KNN in the embedding space")
            nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
            nn.fit(embedding)
            self.embedding_knn = nn.kneighbors_graph(mode="connectivity")

            logging.debug("Correlation Calculation 'full'")
            #####
            self.corrcoef = colDeltaCor(X, delta_X, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                self.corrcoef_random = colDeltaCor(X, delta_X_rndm, threads=threads)

            #####
            np.fill_diagonal(self.corrcoef, 0)
            if calculate_randomized:
                np.fill_diagonal(self.corrcoef_random, 0)

    def calculate_embedding_shift(self, sigma_corr: float=0.05) -> None:
        """Use the transition probability to project the velocity direction on the embedding

        Arguments
        ---------
        sigma_corr: float, default=0.05
            the kernel scaling

        Returns
        -------
        Nothing but it creates the following attributes:
        transition_prob: np.ndarray
            the transition probability calculated using the exponential kernel on the correlation coefficient
        delta_embedding: np.ndarray
            The resulting vector
        """
        # Kernel evaluation
        logging.debug("Calculate transition probability")

        # NOTE maybe sparse matrix here are slower than dense
        # NOTE if knn_random this could be made much faster either using sparse matrix or neigh_ixs
        self.transition_prob = np.exp(self.corrcoef / sigma_corr) * self.embedding_knn.A  # naive
        self.transition_prob /= self.transition_prob.sum(1)[:, None]
        if hasattr(self, "corrcoef_random"):
            logging.debug("Calculate transition probability for negative control")
            self.transition_prob_random = np.exp(self.corrcoef_random / sigma_corr) * self.embedding_knn.A  # naive
            self.transition_prob_random /= self.transition_prob_random.sum(1)[:, None]

        unitary_vectors = self.embedding.T[:, None, :] - self.embedding.T[:, :, None]  # shape (2,ncells,ncells)
        with np.errstate(divide='ignore', invalid='ignore'):
            unitary_vectors /= np.linalg.norm(unitary_vectors, ord=2, axis=0)  # divide by L2
            np.fill_diagonal(unitary_vectors[0, ...], 0)  # fix nans
            np.fill_diagonal(unitary_vectors[1, ...], 0)

        self.delta_embedding = (self.transition_prob * unitary_vectors).sum(2)
        self.delta_embedding -= (self.embedding_knn.A * unitary_vectors).sum(2) / self.embedding_knn.sum(1).A.T
        self.delta_embedding = self.delta_embedding.T


        if hasattr(self, "corrcoef_random"):
            self.delta_embedding_random = (self.transition_prob_random * unitary_vectors).sum(2)
            self.delta_embedding_random -= (self.embedding_knn.A * unitary_vectors).sum(2) / self.embedding_knn.sum(1).A.T
            self.delta_embedding_random = self.delta_embedding_random.T


    def calculate_grid_arrows(self, smooth: float=0.5, steps: Tuple=(40, 40),
                              n_neighbors: int=100, n_jobs: int=4, xylim: Tuple=((None, None), (None, None))) -> None:
        """Calculate the velocity using a points on a regular grid and a gaussian kernel

        Note: the function should work also for n-dimensional grid

        Arguments
        ---------
        embed: str, default=embedding
            The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
            The difference vector is getattr(self, 'delta' + '_' + embed)
        smooth: float, smooth=0.5
            Higher value correspond to taking in consideration further points
            the standard deviation of the gaussian kernel is smooth * stepsize
        steps: tuple, default
            the number of steps in the grid for each axis
        n_neighbors:
            number of neighbors to use in the calculation, bigger number should not change too much the results..
            ...as soon as smooth is small
            Higher value correspond to slower execution time
        n_jobs:
            number of processes for parallel computing
        xymin:
            ((xmin, xmax), (ymin, ymax))

        Returns
        -------
        Nothing but it sets the attributes:
        flow_embedding: np.ndarray
            the coordinates of the embedding
        flow_grid: np.ndarray
            the gridpoints
        flow: np.ndarray
            vector field coordinates
        flow_magnitude: np.ndarray
            magnitude of each vector on the grid
        total_p_mass: np.ndarray
            density at each point of the grid

        """
        embedding = self.embedding
        delta_embedding = getattr(self, f"delta_embedding")

        if hasattr(self, "corrcoef_random"):
            delta_embedding_random = getattr(self, f"delta_embedding_random")

        # Prepare the grid
        grs = []
        for dim_i in range(embedding.shape[1]):
            m, M = np.min(embedding[:, dim_i]), np.max(embedding[:, dim_i])

            if xylim[dim_i][0] is not None:
                m = xylim[dim_i][0]
            if xylim[dim_i][1] is not None:
                M = xylim[dim_i][1]

            m = m - 0.025 * np.abs(M - m)
            M = M + 0.025 * np.abs(M - m)
            gr = np.linspace(m, M, steps[dim_i])
            grs.append(gr)

        meshes_tuple = np.meshgrid(*grs)
        gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

        nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
        nn.fit(embedding)
        dists, neighs = nn.kneighbors(gridpoints_coordinates)

        std = np.mean([(g[1] - g[0]) for g in grs])
        # isotropic gaussian kernel
        gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
        self.total_p_mass = gaussian_w.sum(1)

        UZ = (delta_embedding[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average
        magnitude = np.linalg.norm(UZ, axis=1)
        # Assign attributes
        self.flow_embedding = embedding
        self.flow_grid = gridpoints_coordinates
        self.flow = UZ
        self.flow_norm = UZ / np.percentile(magnitude, 99.5)
        self.flow_norm_magnitude = np.linalg.norm(self.flow_norm, axis=1)

        if hasattr(self, "corrcoef_random"):
            UZ_rndm = (delta_embedding_random[neighs] * gaussian_w[:, :, None]).sum(1) / np.maximum(1, self.total_p_mass)[:, None]  # weighed average
            magnitude_rndm = np.linalg.norm(UZ, axis=1)
            # Assign attributes
            self.flow_rndm = UZ_rndm
            self.flow_norm_rndm = UZ_rndm / np.percentile(magnitude_rndm, 99.5)
            self.flow_norm_magnitude_rndm = np.linalg.norm(self.flow_norm_rndm, axis=1)

    def prepare_markov(self, sigma_D: np.ndarray, sigma_W: np.ndarray, direction: str="forward", cells_ixs: np.ndarray=None) -> None:
        """Prepare a transition probability for Markov process

        Arguments
        ---------
        sigma_D: float
            the standard deviation used on the locality-limiting component
        sigma_W: float
            the standard deviation used on the noise component
        direction: str, default="backwards"
            whether to diffuse forward of backwards
        cells_ixs: np.ndarray, default=None
            Cells to use, if None all the cells will be considered.

        Returns
        -------
        Nothing but it creates the following attributes:
        tr: np.ndarray
            the transition probability matrix

        """
        if cells_ixs is None:
            cells_ixs = np.arange(self.transition_prob.shape[0])

        # NOTE: This implementation is not speed optimized to improve the speed of the implementation:
        # - the C/Fortran contiguity of the transition matrix should be taken into account
        # - a knn implementation would reduce computation
        # - should avoid transformation to and from dense-sparse formats
        if direction == "forward":
            self.tr = np.array(self.transition_prob[cells_ixs, :][:, cells_ixs])
        elif direction == "backwards":
            self.tr = np.array((self.transition_prob[cells_ixs, :][:, cells_ixs]).T, order="C")
        else:
            raise NotImplementedError(f"{direction} is not an implemented direction")
        dist_matrix = squareform(pdist(self.embedding[cells_ixs, :]))
        K_D = gaussian_kernel(dist_matrix, sigma=sigma_D)
        self.tr = self.tr * K_D
        # Fill diagonal with max or the row and sum=1 normalize
        np.fill_diagonal(self.tr, self.tr.max(1))
        self.tr = self.tr / self.tr.sum(1)[:, None]

        K_W = gaussian_kernel(dist_matrix, sigma=sigma_W)
        K_W = K_W / K_W.sum(1)[:, None]
        self.tr = 0.8 * self.tr + 0.2 * K_W
        self.tr = self.tr / self.tr.sum(1)[:, None]
        self.tr = scipy.sparse.csr_matrix(self.tr)

        if hasattr(self, "corrcoef_random"):
            if direction == "forward":
                self.tr_random = np.array(self.transition_prob_random[cells_ixs, :][:, cells_ixs])
            elif direction == "backwards":
                self.tr_random = np.array((self.transition_prob_random[cells_ixs, :][:, cells_ixs]).T, order="C")
            else:
                raise NotImplementedError(f"{direction} is not an implemented direction")
            #dist_matrix = squareform(pdist(self.embedding[cells_ixs, :]))
            #K_D = gaussian_kernel(dist_matrix, sigma=sigma_D)
            self.tr_random = self.tr_random * K_D
            # Fill diagonal with max or the row and sum=1 normalize
            np.fill_diagonal(self.tr_random, self.tr_random.max(1))
            self.tr_random = self.tr_random / self.tr_random.sum(1)[:, None]

            #K_W = gaussian_kernel(dist_matrix, sigma=sigma_W)
            #K_W = K_W / K_W.sum(1)[:, None]
            self.tr_random = 0.8 * self.tr_random + 0.2 * K_W
            self.tr_random = self.tr_random / self.tr_random.sum(1)[:, None]
            self.tr_random = scipy.sparse.csr_matrix(self.tr_random)

    def run_markov(self, starting_p: np.ndarray=None, n_steps: int=2500, mode: str="time_evolution") -> None:
        """Run a Markov process

        Arguments
        ---------
        starting_p: np.ndarray, default=None
            specifies the starting density
            if None is passed an array of 1/self.tr.shape[0] will be created
        n_steps: np.ndarray, default=2500
            Numbers of steps to be performed
        mode: str, default="time_evolution"
            this argument is passed to the Diffusion.diffuse call

        Returns
        -------
        Nothing but it creates the attribute:
        diffused: np.ndarray
            The probability to be found at any of the states
        """
        self.prepare_markov_simulation()

        if starting_p is None:
            starting_p = np.ones(self.tr.shape[0]) / self.tr.shape[0]
        diffusor = Diffusion()
        self.diffused = diffusor.diffuse(starting_p, self.tr, n_steps=n_steps, mode=mode)[0]

    def plot_mc_resutls_as_density(self, args={}):
        """
        this function plot the results of mc chain simulation in velocyto method.
        """
        diffused_n = self.diffused - np.percentile(self.diffused, 3)
        diffused_n /= np.percentile(diffused_n, 97)
        diffused_n = np.clip(diffused_n, 0, 1)
        #diffused_n = tv.diffused
        #plt.figure(None,(7,7))
        args_ = {"alpha": 0.5, "s": 50, "lw": 0., "edgecolor":"", "cmap": "viridis_r", "rasterized": True}
        args_.update(args)

        plt.scatter(self.embedding[self.ixs_mcmc, 0], self.embedding[self.ixs_mcmc, 1],
                        c=diffused_n, **args_)
        plt.axis("off")
        cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        plt.colorbar()
        plt.axis("off")



    def plot_grid_arrows(self, quiver_scale: Union[str, float]="auto", scale_type: str= "relative", min_mass: float=1, min_magnitude: float=None,
                         scatter_kwargs_dict: Dict= None, plot_dots: bool=False, plot_random: bool=False, **quiver_kwargs: Any) -> None:
        """Plots vector field averaging velocity vectors on a grid

        Arguments
        ---------
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
            NOTE: In the situation where "auto" is resulting in very small or big velocities, pass a float to this parameter
            The float will be interpreted as a scaling, importantly both the data and the control will be scaled
            in this way you can rescale the velocity arbitrarily without the risk of observing just an overfit of the noise
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function (not recommended if you are not sure what this implies)
        min_mass: float, default=1
            the minimum density around a grid point for it to be considered and plotted
        min_magnitude: float, default=None
            the minimum magnitude of the velocity for it to be considered and plotted
        scatter_kwargs_dict: dict, default=None
            a dictionary of keyword arguments to pass to scatter
            by default the following are passed: s=20, zorder=-1, alpha=0.2, lw=0, c=self.colorandum. But they can be overridden.
        plot_dots: bool, default=True
            whether to plot dots in correspondence of all low velocity grid points
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.
        """
        # update color information
        col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=self.cluster_column_name,
                                                  return_as="dict")
        self.colorandum = np.array([col_dict[i] for i in self.adata.obs[self.cluster_column_name]])

        # plt.figure(figsize=(10, 10))
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _quiver_kwargs.update(quiver_kwargs)
        scatter_dict = {"s": 20, "zorder": -1, "alpha": 0.2, "lw": 0, "c": self.colorandum}
        if scatter_kwargs_dict is not None:
            scatter_dict.update(scatter_kwargs_dict)

        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "flow_rndm"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.flow_rndm[self.total_p_mass >= min_mass, :], 2, 1), 90)  # Tipical lenght of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.0025)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.0025)
            else:
                raise ValueError(""""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

        mass_filter = self.total_p_mass < min_mass
        if min_magnitude is None:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow)
            if not plot_dots:
                UV = UV[~mass_filter, :]
                XY = XY[~mass_filter, :]
            else:
                UV[mass_filter, :] = 0
        else:
            XY, UV = np.copy(self.flow_grid), np.copy(self.flow_norm)
            if not plot_dots:
                UV = UV[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
                XY = XY[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
            else:
                UV[mass_filter | (self.flow_norm_magnitude < min_magnitude), :] = 0


        if min_magnitude is None:
            XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_rndm)
            if not plot_dots:
                UV_rndm = UV_rndm[~mass_filter, :]
                XY = XY[~mass_filter, :]
            else:
                UV_rndm[mass_filter, :] = 0
        else:
            XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_norm_rndm)
            if not plot_dots:
                UV_rndm = UV_rndm[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                XY = XY[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
            else:
                UV_rndm[mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude), :] = 0

        if plot_random:
            plt.subplot(122)
            plt.title("Randomized")
            plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
            plt.quiver(XY[:, 0], XY[:, 1], UV_rndm[:, 0], UV_rndm[:, 1],
                       scale=quiver_scale, zorder=20000, **_quiver_kwargs)
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")

        plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
        plt.quiver(XY[:, 0], XY[:, 1], UV[:, 0], UV[:, 1],
                   scale=quiver_scale, zorder=20000, **_quiver_kwargs)
        plt.axis("off")

    def plot_arrows_embedding(self, choice: Union[str, int]="auto", quiver_scale: Union[str, float]="auto", scale_type: str="relative",
                              plot_scatter: bool=False, scatter_kwargs: Dict={}, color_arrow: str="cluster",
                              new_fig: bool=False, plot_random: bool=True, **quiver_kwargs: Any) -> None:
        """Plots velocity on the embedding cell-wise

        Arguments
        ---------
        choice: int, default = "auto"
            the number of cells to randomly pick to plot the arrows (To avoid overcrowding)
        quiver_scale: float, default="auto"
            Rescaling factor applied to the arrow field to enhance visibility
            If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
            If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
            NOTE: Despite a similar option than plot_grid_arrows, here there is no strong motivation to calculate the scale relative to the randomized control
            This is because the randomized doesn't have to have smaller velocity cell-wise, there might be for example
            scattered cells that will have strong velocity but they will, correctly just average out when calculating the average velocity field.
        scale_type: str, default="relative"
            How to interpret `quiver_scale`:
            If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
            (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
            If "absolute" the value will be passed to the matplotlib quiver function
        plot_scatter: bool, default = False
            whether to plot the points
        scatter_kwargs: Dict
            A dictionary containing all the keywords arguments to pass to matplotlib scatter
            by default the following are passed: c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3. But they can be overridden.
        color_arrow: str, default = "cluster"
            the color of the arrows, if "cluster" the arrows are colored the same as the cluster
        new_fig: bool, default=False
            whether to create a new figure
        plot_random: bool, default=True
            whether to plot the randomized control next to the plot
        **quiver_kwargs: dict
            keyword arguments to pass to quiver
            By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.

        Returns
        -------
        Nothing, just plots the tsne with arrows
        """
        # update color information
        col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=self.cluster_column_name,
                                                  return_as="dict")
        self.colorandum = np.array([col_dict[i] for i in self.adata.obs[self.cluster_column_name]])

        if choice == "auto":
            choice = int(self.S.shape[1] / 3)
            logging.warning(f"Only {choice} arrows will be shown to avoid overcrowding, you can choose the exact number setting the `choice` argument")
        _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
        _scatter_kwargs = dict(c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3)
        _scatter_kwargs.update(scatter_kwargs)
        if new_fig:
            if plot_random and hasattr(self, "delta_embedding_random"):
                plt.figure(figsize=(22, 12))
            else:
                plt.figure(figsize=(14, 14))

        ix_choice = np.random.choice(self.embedding.shape[0], size=choice, replace=False)

        # Determine quiver scale
        if scale_type == "relative":
            if hasattr(self, "delta_embedding_random"):
                plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0), 2)  # Diagonal of the plot
                arrows_scale = np.percentile(np.linalg.norm(self.delta_embedding_random, 2, 1), 80)  # Tipical length of an arrow
                if quiver_scale == "auto":
                    quiver_scale = arrows_scale / (plot_scale * 0.005)
                else:
                    quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.005)
            else:
                raise ValueError("""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
                Please run estimate_transition_prob or set `scale_type` to `absolute`""")
        else:
            logging.warning("The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

        if color_arrow == "cluster":
            colorandum = self.colorandum[ix_choice, :]
        else:
            colorandum = color_arrow

        _quiver_kwargs.update({"color": colorandum})
        _quiver_kwargs.update(quiver_kwargs)

        if plot_random and hasattr(self, "delta_embedding_random"):
            plt.subplot(122)
            plt.title("Randomized")
            if plot_scatter:
                plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)
            plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                       self.delta_embedding_random[ix_choice, 0], self.delta_embedding_random[ix_choice, 1],
                       scale=quiver_scale, **_quiver_kwargs)
            plt.axis("off")
            plt.subplot(121)
            plt.title("Data")

        if plot_scatter:
            plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)

        plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                   self.delta_embedding[ix_choice, 0], self.delta_embedding[ix_choice, 1],
                   scale=quiver_scale, **_quiver_kwargs)
        plt.axis("off")

    def plot_cell_transitions(self, cell_ix: int=0, alpha: float=0.1, alpha_neigh: float=0.2,
                              cmap_name: str="RdBu_r", plot_arrow: bool=True,
                              mark_cell: bool=True, head_width: int=3) -> None:
        """Plot the probability of a cell to transition to any other cell

        This function is untested
        """
        cmap = plt.cm.get_cmap(name=cmap_name)
        colorandum = np.ones((self.embedding.shape[0], 4))
        colorandum *= 0.3
        colorandum[:, -1] = alpha

        plt.scatter(self.embedding[:, 0], self.embedding[:, 1],
                    c=colorandum, s=50, edgecolor="")
        if mark_cell:
            plt.scatter(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                        facecolor="none", s=100, edgecolor="k")
        if plot_arrow:
            plt.arrow(self.embedding[cell_ix, 0], self.embedding[cell_ix, 1],
                      self.delta_embedding[cell_ix, 0], self.delta_embedding[cell_ix, 1],
                      head_width=head_width, length_includes_head=True)





def scatter_viz(x: np.ndarray, y: np.ndarray, *args: Any, **kwargs: Any) -> Any:
    """A wrapper of scatter plot that guarantees that every point is visible in a very crowded scatterplot

    Args
    ----
    x: np.ndarray
        x axis coordinates
    y: np.ndarray
        y axis coordinates
    args and kwargs:
        positional and keyword arguments as in matplotplib.pyplot.scatter

    Returns
    -------
    Plots the graph and returns the axes object
    """
    ix_x_sort = np.argsort(x, kind="mergesort")
    ix_yx_sort = np.argsort(y[ix_x_sort], kind="mergesort")
    args_new = []
    kwargs_new = {}
    for arg in args:
        if type(arg) is np.ndarray:
            args_new.append(arg[ix_x_sort][ix_yx_sort])
        else:
            args_new.append(arg)
    for karg, varg in kwargs.items():
        if type(varg) is np.ndarray:
            kwargs_new[karg] = varg[ix_x_sort][ix_yx_sort]
        else:
            kwargs_new[karg] = varg
    ax = plt.scatter(x[ix_x_sort][ix_yx_sort], y[ix_x_sort][ix_yx_sort], *args_new, **kwargs_new)
    return ax





@jit(nopython=True)
def numba_random_seed(value: int) -> None:
    """Same as np.random.seed but for numba"""
    np.random.seed(value)


@jit(nopython=True)
def permute_rows_nsign(A: np.ndarray) -> None:
    """Permute in place the entries and randomly switch the sign for each row of a matrix independently.
    """
    plmi = np.array([+1, -1])
    for i in range(A.shape[0]):
        np.random.shuffle(A[i, :])
        A[i, :] = A[i, :] * np.random.choice(plmi, size=A.shape[1])



def gaussian_kernel(X: np.ndarray, mu: float=0, sigma: float=1) -> np.ndarray:
    """Compute gaussian kernel"""
    return np.exp(-(X - mu)**2 / (2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)
