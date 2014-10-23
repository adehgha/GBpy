# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np


def eclose(a,b,rtol=1.0000000000000001e-05, atol=1e-08):
    if a.ndim == 1:
        return np.abs(a - b) <= (atol + rtol * np.abs(b))
        # return abs(a - b) < atol
    else:
        return np.all((np.abs(a - b) <= (atol + rtol * np.abs(b))), axis = 1)
        # return np.all(abs(a - b) < atol, axis=1)
        # return abs(a - b) <= atol + rtol*abs(b)


def unique_rows_tol(xin, tol=1e-12, return_index=False, return_inverse=False):
    """

    :param xin:
    :param tol:
    :param return_index:
    :param return_inverse:
    :return:
    """
    y = xin.copy()
    sz1 = np.shape(xin)[0]
    if xin.ndim == 1:
        ind = np.argsort(y)
        y = y[ind]
    else:
        sz2 = np.shape(xin)[1]
        sort_tup = []
        for i in range(sz2):
            sort_tup.append(y[:, sz2-i-1])

        ind = np.lexsort(sort_tup)
        y = y[ind, :]

    ci = 0
    u = np.empty((0, sz2), dtype=y.dtype)

    if return_index:
        ret_ind = []
    if return_index:
        ret_inv = []

    while ci < sz1:
        if xin.ndim == 1:
            y1 = np.tile(y[ci], (sz1, ))
        else:
            y1 = np.tile(y[ci, :], (sz1, 1))

        ii = eclose(y1, y, tol)
        mi = np.max(ii.nonzero())
        if return_index:
            ret_ind.append(ind[mi])
        if return_inverse:
            ret_inv.append(ind[ii.nonzero()])
        u = np.concatenate((u, [y[mi]]))
        ci = mi + 1

    if not return_index and not return_inverse:
        return u
    else:
        if return_index and return_inverse:
            return u, ret_ind, ret_inv
        elif return_index:
            return u, ret_ind
        elif return_inverse:
            return u, ret_inv



# A1 = np.random.rand(10, 4)
# A2 = A1 + 1e-8
# A3 = np.concatenate((A1, A2))
# A4 = unique_rows_tol(A3, 1e-5, return_index=True)

# print A4

