
import numpy as np

from tqdm.auto import tqdm

#from network_fitting import TransNet as tn
from numba import jit, f8, i8, void


#############################################
### custom function for markov simulation ###
#############################################


#@jit(i8[:,:](i8[:], f8[:,:], i8))
def _walk(start_cell_id_array, transition_prob, n_steps):
    """
    Do markov simulation based on given transition probability.
    Transition probability suppose to be NxN 2d numpy array that
    stands for transition prob between N states.

    Args:
        start_cell_id_array (1d numpy array (int)): position for initiation

        transition_prob (2d numpy array (int)): transition probabiliry matrix

        n_steps (int) : numbers of steps

    Returns:
        2d numpy array: trajectory of idices saved in a matrix.
    """
    n_cells = transition_prob.shape[0]

    # define starting cell index.
    ids_now = start_cell_id_array
    li = []

    li.append(list(ids_now))
    ids_unique = list(range(n_cells))

    # walk for n_steps
    for i in range(n_steps):
        # calculate next position for all cells
        li_ = []
        for j in ids_now:
            choiced = ids_unique[np.searchsorted(np.cumsum(transition_prob[j,:]), np.random.random(), side="right")]
            #choiced = np.random.choice(a=ids_unique, replace=True, p=transition_prob[j,:])
            li_.append(choiced)
        # update current position
        ids_now = np.array(li_)
        li.append(list(ids_now))

    #li = tuple(li)
    trajectory_matrix = np.array(li).transpose()
    return trajectory_matrix



@jit(nopython=True)
def numba_random_seed(value: int) -> None:
    """Same as np.random.seed but for numba"""
    np.random.seed(value)
