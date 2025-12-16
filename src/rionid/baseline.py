import numpy as np
import math
import scipy.special as sp
from scipy import sparse
from scipy.sparse.linalg import spsolve

class NONPARAMS_EST(object):
    def __init__(self, data):
        self.data = data

    def pls(self, method, l, **kwargs):
        pls_method = PLS(self.data)
        if method == 'BrPLS':
            return pls_method.BrPLS(l=l, **kwargs)
        else:
            raise ValueError(f"Method {method} not implemented.")

class PLS(object):
    def __init__(self, data):
        self.data = data

    def BrPLS(self, l, ratio=1e-6, nitermax=50):
        L, beta = len(self.data), 0.5
        D = sparse.diags([1,-2,1], [0,-1,-2], shape=(L, L-2))
        D = l * D.dot(D.transpose())
        w, z = np.ones(L), self.data.copy()
        
        for i in range(nitermax):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + D
            zt = spsolve(Z, w*self.data)
            d = self.data - zt
            d_pos = d[d > 0]
            d_neg = d[d < 0]
            d_m = np.mean(d_pos) if len(d_pos) > 0 else 1e-6
            d_sigma = np.sqrt(np.mean(d_neg**2)) if len(d_neg) > 0 else 1e-6
            
            w = 1 / (1 + beta / (1 - beta) * np.sqrt(np.pi / 2) * d_sigma / d_m * 
                     (1 + sp.erf((d / d_sigma - d_sigma / d_m) / np.sqrt(2))) * 
                     np.exp((d / d_sigma - d_sigma / d_m)**2 / 2))
            
            if np.sqrt(np.sum((z - zt)**2) / np.sum(z**2)) < ratio: break
            z = zt
            if np.abs(beta + np.mean(w) - 1.) < ratio: break
            beta = 1 - np.mean(w)
        return z