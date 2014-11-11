# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>
#               Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np


def eclose(a, b, rtol=1.0000000000000001e-05, atol=1e-08):
    """
    """
    if a.ndim == 1:
        return np.abs(a - b) <= (atol + rtol * np.abs(b))
    else:
        return np.all((np.abs(a - b) <= (atol + rtol * np.abs(b))),
                      axis=1)


def unique_rows_tol(xin, tol=1e-12, return_index=False, return_inverse=False):
    """

    :param xin:
    :param tol:
    :param return_index:
    :param return_inverse:
    :return:
    """
    sz1 = np.shape(xin)[0]
    if xin.ndim == 1:
        ind = np.argsort(xin)
        y = xin[ind]
        y = np.reshape(y, (sz1, 1))
    else:
        y = xin.copy()
        sz2 = np.shape(y)[1]
        # sort_tup = []
        # for i in range(sz2):
        #     sort_tup.append(y[:, sz2 - i - 1])
        # ind = np.lexsort((y[:, 1], y[:, 0]))
        # y1 = xin[ind, 0:2]
        # ind = np.lexsort(sort_tup)
        # y1 = xin[ind, :]
        ind = np.lexsort(xin.T)
        y1 = xin[np.lexsort(xin.T)]

    ci = 0
    u = np.empty((0, sz2), dtype=y.dtype)

    if return_index:
        ret_ind = []
    if return_inverse:
        ret_inv = np.zeros((sz1), dtype=int)

    ct1 = 0
    while ci < sz1:
        # print ci
        y2 = np.tile(y1[ci, :], (sz1, 1))
        ii = eclose(y2, y1, tol)
        # ii = eclose(y[ci, :], y)
        mi = np.max(ii.nonzero())
        u = np.concatenate((u, [y[ind[mi], :]]))
        if return_index:
            ret_ind.append(ind[mi])
        if return_inverse:
            ret_inv[ind[ii]] = ct1
        ct1 += 1
        ci = mi + 1

    if not return_index and not return_inverse:
        return u
    else:
        if return_index and return_inverse:
            return u, ret_ind, ret_inv.tolist()
        elif return_index:
            return u, ret_ind
        elif return_inverse:
            return u, ret_inv.tolist()

# def unique_rows_tol(data, prec):
#     d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
#     b = np.ascontiguousarray(d_r).view(np.dtype((np.void,
#                                                  d_r.dtype.itemsize
#                                                  * d_r.shape[1])))
#     return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1])

# import pickle

# pkl_file = ('bpnormals.pkl')
# jar1 = open(pkl_file, 'rb')
# bpnormals_go1 = pickle.load(jar1)

# unq_bpn_go1, ia, ic = unique_rows_tol(bpnormals_go1, 1e-12, True, True)

# unq_bpn_go2 = unique_rows(bpnormals_go1, 12)


# print unq_bpn_go1

# print unq_bpn_go2
