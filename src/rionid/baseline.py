import numpy as np
import math
import scipy.special as sp
from scipy import sparse
from scipy.sparse.linalg import spsolve

class NONPARAMS_EST(object):
    """
    Wrapper class for non-parametric baseline estimation algorithms.

    This class provides a unified interface to access various baseline subtraction
    methods, primarily Penalized Least Squares (PLS) variants.

    Parameters
    ----------
    data : array-like
        The input 1D spectrum array (e.g., Power Spectral Density) from which
        the baseline should be estimated.
    """
    def __init__(self, data):
        self.data = data

    def pls(self, method, l, **kwargs):
        """
        Executes a Penalized Least Squares (PLS) baseline estimation.

        Parameters
        ----------
        method : str
            The specific PLS algorithm to use. Currently supports 'BrPLS'.
        l : float
            The smoothness parameter (lambda). Higher values result in a stiffer,
            smoother baseline. Typical values range from 10^5 to 10^8 for Schottky spectra.
        **kwargs : dict
            Additional arguments passed to the specific PLS method (e.g., 'ratio', 'nitermax').

        Returns
        -------
        numpy.ndarray
            The estimated baseline array, with the same shape as the input data.

        Raises
        ------
        ValueError
            If the specified `method` is not implemented.
        """
        pls_method = PLS(self.data)
        if method == 'BrPLS':
            return pls_method.BrPLS(l=l, **kwargs)
        else:
            raise ValueError(f"Method {method} not implemented.")

class PLS(object):
    """
    Implementation of Penalized Least Squares (PLS) algorithms.

    Parameters
    ----------
    data : array-like
        The input 1D spectrum array.
    """
    def __init__(self, data):
        self.data = data

    def BrPLS(self, l, ratio=1e-6, nitermax=50):
        """
        Bayesian reweighted Penalized Least Squares (BrPLS).

        This algorithm estimates the baseline by iteratively reweighting the data points.
        It assumes that the baseline is smooth and that peaks are positive deviations
        from this baseline. It uses a Bayesian approach to update weights based on
        the probability that a point belongs to the background noise versus a peak.

        Reference: Q. Wang et al., NUCL SCI TECH, 33: 148 (2022).

        Parameters
        ----------
        l : float
            Smoothness parameter (lambda). Controls the trade-off between fidelity to
            the data and smoothness of the baseline.
        ratio : float, optional
            Convergence threshold. The iteration stops when the relative change in
            the baseline vector is less than this value. Default is 1e-6.
        nitermax : int, optional
            Maximum number of iterations. Default is 50.

        Returns
        -------
        numpy.ndarray
            The calculated baseline array.
        """
        L, beta = len(self.data), 0.5
        
        # Construct the difference matrix (2nd order derivative approximation)
        D = sparse.diags([1,-2,1], [0,-1,-2], shape=(L, L-2))
        D = l * D.dot(D.transpose())
        
        w, z = np.ones(L), self.data.copy()
        
        for i in range(nitermax):
            # Solve the linear system (W + D)z = Wy
            W = sparse.spdiags(w, 0, L, L)
            Z = W + D
            zt = spsolve(Z, w*self.data)
            
            # Calculate residuals
            d = self.data - zt
            
            # Separate positive and negative residuals to estimate noise statistics
            d_pos = d[d > 0]
            d_neg = d[d < 0]
            
            # Robust mean and sigma estimation (avoid division by zero)
            d_m = np.mean(d_pos) if len(d_pos) > 0 else 1e-6
            d_sigma = np.sqrt(np.mean(d_neg**2)) if len(d_neg) > 0 else 1e-6
            
            # Bayesian weight update formula using Error Function (erf)
            # This calculates the probability that a point is part of the baseline
            w = 1 / (1 + beta / (1 - beta) * np.sqrt(np.pi / 2) * d_sigma / d_m * 
                     (1 + sp.erf((d / d_sigma - d_sigma / d_m) / np.sqrt(2))) * 
                     np.exp((d / d_sigma - d_sigma / d_m)**2 / 2))
            
            # Check convergence of the baseline vector z
            if np.sqrt(np.sum((z - zt)**2) / np.sum(z**2)) < ratio: break
            z = zt
            
            # Update prior probability beta
            if np.abs(beta + np.mean(w) - 1.) < ratio: break
            beta = 1 - np.mean(w)
            
        return z