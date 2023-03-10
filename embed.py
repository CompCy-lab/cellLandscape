import numpy as np
import anndata
from typing import Union
import scipy
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

def random_feats(X: np.ndarray,
                gamma: Union[int, float] = 1,
                D: int = 2000, 
                frequency_seed: int = None):
    """Computes random Fourier frequency features: https://papers.nips.cc/paper/2007/hash/013a006f03dbc5392effeb8f18fda755-Abstract.html

    Parameters
    X: np.ndarray
        array of input data (dimensions = cells x features)
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution  
    D: int (default = 2000):
        dimensionality of the random Fourier frequency features, D/2 sin and D/2 cos basis  
    frequency_seed: int (default = None):
        random state parameter  
    ----------

    Returns
    phi: np.ndarray
        random Fourier frequency features (dimensions = cells x D)
    ----------
    """
    scale = 1 / gamma
    d = int(np.floor(D / 2))

    if (frequency_seed is not None):
        np.random.seed(frequency_seed)
        W = np.random.normal(scale = scale, size = (X.shape[1], d))
    else:
        W = np.random.normal(scale = scale, size = (X.shape[1], d))

    XW = np.dot(X, W)
    sin_XW = np.sin(XW)
    cos_XW = np.cos(XW)
    phi = np.concatenate((cos_XW, sin_XW), axis=1)

    return phi

def mean_embed(adata,
            sample_set_key: str = None,
            gamma: Union[int, float] = 1,
            frequency_seed: int = None,
            D: int = 2000,
            n_jobs: int = -1):
    """Computes the mean embedding of the random Fourier frequency features

    Parameters
    adata: anndata.Anndata
        annotated data object (dimensions = cells x features)
    sample_set_key: str (default = None)
        string referring to the key within adata.obs that contains the samples to compute the embedding
            ~ if sample_set_key is None, will use all cells as a single sample 
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution within random Fourier frequency feature computation
    frequency_seed: int (default = None):
        random state parameter   
    D: int (default = 2000)
        dimensionality of the random Fourier frequency features, D/2 sin and D/2 cos basis 
    n_jobs (default = -1)
        number of tasks
    ----------
    
    Returns
    X_embed: np.array
        mean embedding (dimensions = sample_sets x D)
    ----------
    """
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    elif n_jobs < -1:
        n_jobs = mp.cpu_count() + 1 + n_jobs

    if sample_set_key is not None:
        sample_set_id, idx = np.unique(adata.obs[sample_set_key], return_index = True)
        sample_set_id = sample_set_id[np.argsort(idx)]
        sample_set_inds = [np.where(adata.obs[sample_set_key] == i)[0] for i in sample_set_id]
    else:
        sample_set_inds = [np.arange(0, adata.X.shape[0])]

    X = _parse_input(adata)

    n_sample_sets = len(sample_set_inds)
    p = mp.Pool(n_jobs)
        
    X_embed = []
    for result in tqdm(p.imap(partial(_run_mean_embed, X = X, gamma = gamma, frequency_seed = frequency_seed, D = D), sample_set_inds), total = n_sample_sets, desc = 'computing mean embedding'):
        X_embed.append(result)
    
    X_embed = np.concatenate(X_embed, 0)

    return X_embed

def _run_mean_embed(sample_set_ind,
                X: np.ndarray = None,
                gamma: Union[int, float] = 1,
                frequency_seed: int = None,
                D: int = 2000):
    """Computes the mean embedding of the Fourier frequency features for a single sample_set

    Parameters
    X: np.array
        data matrix (dimensions = cells x features)
    sample_set_ind: np.array
        indicies of cells from sample_sets of interest
    gamma: Union([int, float]) (default = 1)
        scale for standard deviation of the normal distribution within random Fourier frequency feature computation
    frequency_seed: int (default = None):
        random state parameter   
    D: int (default = 2000)
        dimensionality of the random Fourier frequency features, D/2 sin and D/2 cos bases 
    ----------

    Returns
    X_embed: np.array
        mean embedding for sample_set of interest (dimensions = 1 x D)
    ----------
    """
    X_ = X[sample_set_ind, :].copy()
    phi = random_feats(X_, gamma = gamma, frequency_seed = frequency_seed, D = D)
    f = np.mean(phi, 0)
    X_embed = f[None, :]

    return X_embed

def _parse_input(adata: anndata.AnnData):
    """accesses and parses data from adata object

    Parameters
    adata: anndata.AnnData
        annotated data object where adata.X is the attribute for preprocessed data
    ----------

    Returns
    X: np.ndarray
        array of data (dimensions = cells x features)
    ----------
    """
    try:
        if isinstance(adata, anndata.AnnData):
            X = adata.X.copy()
        if isinstance(X, scipy.sparse.csr_matrix):
            X = np.asarray(X.todense())
    except NameError:
        pass
    
    return X