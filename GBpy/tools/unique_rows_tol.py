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
    sz1 = np.shape(xin)[0]
    if xin.ndim == 1:
        ind = np.argsort(xin)
        y = xin[ind]
        y = np.reshape(y, (sz1, 1))
    else:
        y = xin.copy()
        sz2 = np.shape(y)[1]
        sort_tup = []
        for i in range(sz2):
            sort_tup.append(y[:, sz2-i-1])
        ind = np.lexsort(sort_tup)
        y1 = y[ind, :]

    ci = 0
    u = np.empty((0, sz2), dtype=y.dtype)

    if return_index:
        ret_ind = []
    if return_inverse:
        ret_inv = np.zeros((sz1), dtype=int)

    ct1 = 0
    while ci < sz1:
        ## print ci
        y2 = np.tile(y1[ci, :], (sz1, 1))
        ii = eclose(y2, y1, tol)
        ## ii = eclose(y[ci, :],y)
        mi = np.max(ii.nonzero())
        u = np.concatenate((u,[y[ind[mi], :]]))
        ret_ind.append(ind[mi])
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


# def unique_rows_tol(xin, tol=1e-12, return_index=False, return_inverse=False):
#     """

#     :param xin:
#     :param tol:
#     :param return_index:
#     :param return_inverse:
#     :return:
#     """

#     # if xin.ndim == 1:
#     #     ind = np.argsort(y)
#     #     y = y[ind]
#     # else:
#     #     sz2 = np.shape(xin)[1]
#     #     sort_tup = []
#     #     for i in range(sz2):
#     #         sort_tup.append(y[:, sz2-i-1])

#     #     ind = np.lexsort(sort_tup)
#     #     y = y[ind, :]


#     sz1 = np.shape(xin)[0]

#     if xin.ndim == 1:
#         y = np.reshape(xin, (sz1, 1))
#     else:
#         y = xin.copy()

#     ci = 0
#     # u = np.empty((0, sz2), dtype=y.dtype)

#     if return_index:
#         ret_ind = []
#     if return_inverse:
#         ret_inv = []

#     y1 = y.copy()
#     y3 = np.zeros(np.shape(y1))

#     while sz1 > 0:
#         y2 = np.tile(y1[0, :], (sz1, 1))
#         ii = eclose(y2, y1, tol)
#         mi = ii.nonzero()[0]
#         y3[ci, :] = y1[mi[0], :]
#         ci += 1
#         y1 = np.delete(y1, mi, 0)
#         sz1 = np.shape(y1)[0]
        
#         u = np.delete(y3, np.arange(ci, np.shape(y)[0]), 0)
    
#     if return_index:
#         sz1 = np.shape(u)[0]
#         sz2 = np.shape(y)[0]
#         ret_ind = np.zeros((sz1, ))
#         for ct1 in range(sz1):
#             y2 = np.tile(u[ct1, :], (sz2, 1))
#             ii = eclose(y, y2, tol)
#             mi = ii.nonzero()[0]
#             ret_ind[ct1] = mi[0]

#     if return_inverse:
#         sz1 = np.shape(u)[0]
#         sz2 = np.shape(y)[0]
#         ret_inv = np.zeros((sz2, ))
#         for ct1 in range(sz1):
#             y2 = np.tile(u[ct1, :], (sz2, 1))
#             ii = eclose(y, y2, tol)
#             mi = ii.nonzero()[0]
#             ret_inv[mi] = ct1


#     if not return_index and not return_inverse:
#         return u
#     else:
#         if return_index and return_inverse:
#             return u, ret_ind, ret_inv
#         elif return_index:
#             return u, ret_ind
#         elif return_inverse:
#             return u, ret_inv



# A1 = np.random.rand(10, 4)
# A2 = A1 + 1e-8
# A3 = np.concatenate((A1, A2))
# A4 = unique_rows_tol(A3, 1e-5, return_index=True)

# print A4

