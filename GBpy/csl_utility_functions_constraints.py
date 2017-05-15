# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
# import sys
# import os
import integer_manipulations as int_man
import misorient_fz as mis_fz
import tools as trans
# -----------------------------------------------------------------------------------------------------------


def proper_ptgrp(cryst_ptgrp):
    """
    Returns the proper point group corresponding to a crystallographic point
    group

    Parameters
    ----------------
    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    Returns
    ----------
    proper_ptgrp: string
        Proper point group in Schoenflies notation
    """
    if cryst_ptgrp in ['D3', 'D3d']:
        proper_ptgrp = 'D3'
    if cryst_ptgrp in ['D4', 'D4h']:
        proper_ptgrp = 'D4'
    if cryst_ptgrp in ['D6', 'D6h']:
        proper_ptgrp = 'D6'
    if cryst_ptgrp in ['O', 'Oh']:
        proper_ptgrp = 'O'

    # prop_grps = ['C1', 'C2', 'C3', 'C4', 'C6', 'D2', 'D3', 'D4', 'D6',
    #              'T', 'O']
    # laue_grps = ['Ci', 'C2h', 'C3i', 'C4h', 'C6h', 'D2h', 'D3d', 'D4h', 'D6h',
    #              'Th', 'Oh']

    # if cryst_ptgrp in laue_grps:
    #     proper_ptgrp = 
    # elif cryst_ptgrp in prop_grps:
    #     proper_ptgrp = cryst_ptgrp

    return proper_ptgrp
# -----------------------------------------------------------------------------------------------------------


def largest_odd_factor(var_arr):
    """
    Function that computes the larges odd factors of an array of integers

    Parameters
    -----------------
    var_arr: numpy array
        Array of integers whose largest odd factors needs to be computed

    Returns
    ------------
    odd_d: numpy array
        Array of largest odd factors of each integer in var_arr
    """
    if var_arr.ndim == 1:
        odd_d = np.empty(np.shape(var_arr))
        odd_d[:] = np.NaN

        ind1 = np.where((np.remainder(var_arr, 2) != 0) | (var_arr == 0))[0]
        if np.size(ind1) != 0:
            odd_d[ind1] = var_arr[ind1]

        ind2 = np.where((np.remainder(var_arr, 2) == 0) & (var_arr != 0))[0]
        if np.size(ind2) != 0:
            odd_d[ind2] = largest_odd_factor(var_arr[ind2] / 2.0)
        return odd_d
    else:
        raise Exception('Wrong Input Type')
# -----------------------------------------------------------------------------------------------------------


def compute_inp_params(lattice, sig_type):
    """
    tau and kmax necessary for possible integer quadruple combinations
    are computed

    Parameters
    ----------------
    lattice: Lattice class
        Attributes of the underlying lattice

    sig_type: {'common', 'specific'}

    Returns
    -----------
    tau: float
        tau is a rational number :math:`= \\frac{\\nu}{\\mu}`

    kmax: float
        kmax is an integer that depends on :math:`\\mu \\ , \\nu`
    """
    lat_params = lattice.lat_params
    cryst_ptgrp = proper_ptgrp(lattice.cryst_ptgrp)

    if cryst_ptgrp == 'D3':
        c_alpha = np.cos(lat_params['alpha'])
        tau = c_alpha / (1 + 2 * c_alpha)
        if sig_type == 'specific':
            [nu, mu] = int_man.rat(tau)
            rho = mu - 3 * nu
            kmax = 4 * mu * rho
        elif sig_type == 'common':
            kmax = []

    if cryst_ptgrp == 'D4':
        tau = (lat_params['a'] ** 2) / (lat_params['c'] ** 2)
        if sig_type == 'specific':
            [nu, mu] = int_man.rat(tau)
            kmax = 4 * mu * nu
        if sig_type == 'common':
            kmax = []

    if cryst_ptgrp == 'D6':
        tau = (lat_params['a'] ** 2) / (lat_params['c'] ** 2)
        if sig_type == 'specific':
            [nu, mu] = int_man.rat(tau)
            if np.remainder(nu, 2) == 0:
                if np.remainder(nu, 4) == 0:
                    kmax = 3 * mu * nu
                else:
                    kmax = 6 * mu * nu
            else:
                kmax = 12 * mu * nu
        if sig_type == 'common':
            kmax = []

    if cryst_ptgrp == 'O':
        tau = 1
        kmax = []

    return tau, kmax
# -----------------------------------------------------------------------------------------------------------


def mesh_muvw(cryst_ptgrp, sigma, sig_type, *args):
    """
    Compute max allowed values of [m,U,V,W] and generates an array
    of integer quadruples

    Parameters
    ----------------
    cryst_ptgrp: string
        Proper point group in Schoenflies notation

    sigma: integer
        Sigma number

    sig_type: {'common', 'specific'}

    args[0]: dictionary
    keys: 'nu', 'mu', 'kmax'

    Returns
    -----------
    Integer quadruple numpy array
    """
    if sig_type == 'common':
        if cryst_ptgrp == 'D3':
            tu1 = np.ceil(2 * np.sqrt(sigma))
            m_max = tu1
            u_max = tu1
            v_max = tu1
            w_max = tu1
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [-v_max, v_max]
            wlims = [0, w_max]

        if cryst_ptgrp == 'D6':
            tu1 = np.ceil(np.sqrt(sigma / 3.0))
            tu2 = np.ceil(np.sqrt(sigma))
            m_max = tu1
            u_max = tu2
            v_max = tu2
            w_max = tu2
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [0, v_max]
            wlims = [0, w_max]

        if cryst_ptgrp == 'D4' or cryst_ptgrp == 'O':
            t1 = np.ceil(np.sqrt(sigma))
            m_max = t1
            u_max = t1
            v_max = t1
            w_max = t1
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [0, v_max]
            wlims = [0, w_max]

    elif sig_type == 'specific':
        mu = args[0]['mu']
        nu = args[0]['nu']
        kmax = args[0]['kmax']

        if cryst_ptgrp == 'D3':
            t1 = np.ceil(np.sqrt(sigma * kmax / (mu)))
            t2 = np.ceil(np.sqrt(sigma * kmax / (mu - 2 * nu)))
            m_max = t1
            u_max = t2
            v_max = t2
            w_max = t2
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [-v_max, v_max]
            wlims = [-w_max, w_max]

        if cryst_ptgrp == 'D6':
            m_max = np.ceil(np.sqrt(sigma * kmax / (3.0 * mu)))
            u_max = np.ceil(np.sqrt(sigma * kmax / (nu)))
            v_max = np.ceil(np.sqrt(sigma * kmax / (nu)))
            w_max = np.ceil(np.sqrt(sigma * kmax / (mu)))
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [0, v_max]
            wlims = [0, w_max]

        if cryst_ptgrp == 'D4':
            t1 = np.sqrt(sigma * kmax)
            m_max = np.ceil(t1 / np.sqrt(mu))
            u_max = np.ceil(t1 / np.sqrt(nu))
            v_max = np.ceil(t1 / np.sqrt(nu))
            w_max = np.ceil(t1 / np.sqrt(mu))
            mlims = [0, m_max]
            ulims = [0, u_max]
            vlims = [0, v_max]
            wlims = [0, w_max]
    else:
        raise Exception('sig_type: wrong input type')

    m_var = np.arange(mlims[0], mlims[1] + 1, 1)
    u_var = np.arange(ulims[0], ulims[1] + 1, 1)
    v_var = np.arange(vlims[0], vlims[1] + 1, 1)
    w_var = np.arange(wlims[0], wlims[1] + 1, 1)

    [x1, x2, x3, x4] = np.meshgrid(m_var, u_var, v_var, w_var)

    x1 = x1.ravel()
    x2 = x2.ravel()
    x3 = x3.ravel()
    x4 = x4.ravel()

    return np.vstack((x1, x2, x3, x4)).astype(int)
# -----------------------------------------------------------------------------------------------------------


def mesh_muvw_fz(quad_int, cryst_ptgrp, sig_type, *args):
    """
    For given integer quadruples, the set belonging to the corresponding
    fundamental zone are separated out and retruned.

    Parameters
    ----------------
    quad_int: numpy array
        Integer quadruples

    cryst_ptgrp: string
        Proper point group in Schoenflies notation

    sig_type: {'common', 'specific'}

    args[0]: dictionary
    keys: 'nu', 'mu', 'kmax'

    Returns
    -----------
    Integer quadruple numpy array belonging to the fundamental zone
    of the corresponding crystallographic point group
    """
    m = quad_int[0, :]
    u = quad_int[1, :]
    v = quad_int[2, :]
    w = quad_int[3, :]

    if sig_type == 'specific':
        if cryst_ptgrp == 'D3':
            mu = args[0]['mu']
            nu = args[0]['nu']
            tau = float(nu)/float(mu)

            cond0 = u + v + w >= 0
            cond1 = (u >= w) & (v >= w)
            condfin = cond0 & cond1
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

            cond0 = 2*m >= np.sqrt(1 - 3*tau)*(u + w)
            cond1 = m >= u + v + w
            condfin = cond0 & cond1
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

        if cryst_ptgrp == 'D4':
            mu = args[0]['mu']
            nu = args[0]['nu']
            tau = float(nu)/float(mu)
            cond0 = (u >= v)
            cond1 = (m >= (np.sqrt(2) + 1) * w)
            condfin = (cond0 & cond1)
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

            cond0 = (m >= np.sqrt(tau) * u)
            cond1 = (m >= np.sqrt(tau / 2) * (u + v))
            condfin = (cond0 & cond1)
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

        if cryst_ptgrp == 'D6':
            mu = args[0]['mu']
            nu = args[0]['nu']
            cond0 = (u >= 2 * v)
            cond1 = (m >= (2 / np.sqrt(3) + 1) * w)
            condfin = (cond0 & cond1)
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

            condfin = (m >= (np.sqrt(nu / (4 * mu)) * u))
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

            condfin = (m >= (np.sqrt(nu / (12 * mu)) * (u - 2 * v)))
            m = m[condfin]
            u = u[condfin]
            v = v[condfin]
            w = w[condfin]

        return np.vstack((m, u, v, w))
# -----------------------------------------------------------------------------------------------------------

def check_fsig_int(quad_int, cryst_ptgrp, sigma, *args):
    """
    For specific sigma rotations, a function of m, U, V, W (fsig) is computed.
    The ratio of fsig and sigma should be a divisor of kmax. This
    condition is checked and those integer quadruples that satisfy
    this condition are returned

    Parameters
    ----------------
    quad_int: numpy array
        Integer quadruples

    cryst_ptgrp: string
        Proper point group in Schoenflies notation

    sigma: float
        sigma number

    args[0]: dictionary
    keys: 'nu', 'mu', 'kmax'

    Returns
    -----------
    quad_int: numpy array
    Integer quadruple array that satisfy the above mentioned condition
    """
    mu = args[0]['mu']
    nu = args[0]['nu']
    kmax = args[0]['kmax']

    m = quad_int[0, :]
    u = quad_int[1, :]
    v = quad_int[2, :]
    w = quad_int[3, :]

    sigma = float(sigma)

    if cryst_ptgrp == 'D3':
        # $\frac{F}{$\Sigma$}$ should be a divisor of kmax
        # $\in (12\mu\nu, 6\mu\nu, 3\mu\nu)$
        # Keep only those quadruples for which the above condition is met
        fsig = ((mu * (m ** 2) + (mu - 2 * nu) * (u ** 2 + v ** 2 + w ** 2)
                 + 2 * nu * (u * v + v * w + w * u)) / sigma)
        cond1 = np.where(abs(fsig - np.round(fsig)) < 1e-06)[0]
        cond2 = np.where(np.remainder(kmax, fsig[cond1]) == 0)[0]
        quad_int = quad_int[:, cond1[cond2]]

    if cryst_ptgrp == 'D4':
        # $\frac{F}{$\Sigma$}$ should be a divisor of kmax
        # $\in (12\mu\nu, 6\mu\nu, 3\mu\nu)$
        # Keep only those quadruples for which the above condition is met
        fsig = (mu * (m ** 2 + w ** 2) + nu * (u ** 2 + v ** 2)) / sigma
        cond1 = np.where(abs(fsig - np.round(fsig)) < 1e-06)[0]
        cond2 = np.where(np.remainder(kmax, fsig[cond1]) == 0)[0]
        quad_int = quad_int[:, cond1[cond2]]

    if cryst_ptgrp == 'D6':
        # $\frac{F}{$\Sigma$}$ should be a divisor of kmax
        # $\in (12\mu\nu, 6\mu\nu, 3\mu\nu)$
        # Keep only those quadruples for which the above condition is met
        fsig = ((mu * (3 * (m ** 2) + w ** 2) +
                 nu * (u ** 2 - u * v + v ** 2)) / sigma)
        cond1 = np.where(abs(fsig - np.round(fsig)) < 1e-06)[0]
        cond2 = np.where(np.remainder(kmax, fsig[cond1]) == 0)[0]
        quad_int = quad_int[:, cond1[cond2]]

    return quad_int
# -----------------------------------------------------------------------------------------------------------

def eliminate_idrots(quad_int):
    """
    Eliminate the roations that belong to the identity matrix and return the
    integer quadruples
    """
    m = quad_int[0, :]
    u = quad_int[1, :]
    v = quad_int[2, :]
    w = quad_int[3, :]

    cond0 = (u == 0) & (v == 0) & (w == 0)
    cond1 = (m == 0) & (v == 0) & (w == 0)
    cond2 = (u == 0) & (m == 0) & (w == 0)
    cond3 = (u == 0) & (v == 0) & (m == 0)

    condfin = (cond0 | cond1 | cond2 | cond3)
    quad_int = np.delete(quad_int, np.where(condfin), axis=1)

    return quad_int
# -----------------------------------------------------------------------------------------------------------

def sigtype_muvw(quad_int, cryst_ptgrp, sig_type):
    """
    The type of integer quadruples are different for common and specific sigma
    rotations. For example, for D4 point group, common rotations satisfy the
    condition u = 0 and v = 0 or m = 0 and w = 0. The specific rotations belong
    to the complimentary set. Depending on the sig_type (common, specific), the
    appropriate set of the integer quadruples are returned.

    Parameters
    ----------------
    quad_int: numpy array
        Integer quadruples.

    cryst_ptgrp: string
        Proper point group in Schoenflies notation.

    sig_type: {'common', 'specific'}

    Returns
    -----------
    quad_int: numpy array
        Integer quadruple array that satisfy the above mentioned condition.
    """

    m = quad_int[0, :]
    u = quad_int[1, :]
    v = quad_int[2, :]
    w = quad_int[3, :]

    if cryst_ptgrp == 'D3':
        cond0 = (m == 0) & (u + v + w == 0)
        cond1 = (u == v) & (v == w)
        condfin = cond0 | cond1
    if cryst_ptgrp == 'D4' or cryst_ptgrp == 'D6':
        cond0 = (u == 0) & (v == 0)
        cond1 = (m == 0) & (w == 0)
        condfin = cond0 | cond1
    if cryst_ptgrp == 'O':
        cond0 = m >= u
        cond1 = u >= v
        cond2 = v >= w
        condfin = cond0 & cond1 & cond2

    if sig_type == 'specific':
        condfin = ~condfin

    return quad_int[:, condfin]
# -----------------------------------------------------------------------------------------------------------

def eliminate_mults(quad_int):
    """
    Divide all the integer quadruples by their corresponding least common
    multiples and return the unique set of integer quadruples
    """
    quad_gcd = int_man.gcd_array(quad_int.astype(int), 'columns')
    quad_gcd = np.tile(quad_gcd, (4, 1))

    a = quad_int / quad_gcd
    a = a.transpose()
    b = np.ascontiguousarray(a).view(np.dtype((np.void,
                                               a.dtype.itemsize * a.shape[1])))
    quad_int = np.unique(b).view(a.dtype).reshape(-1, a.shape[1])
    quad_int = quad_int.transpose()

    return quad_int
# -----------------------------------------------------------------------------------------------------------

def check_sigma(quad_int, sigma, cryst_ptgrp, sig_type, *args):
    """
    The integer quadruples that correspond to a sigma rotation satisfy
    certain conditions. These conditions are checked and all the
    quadruples that do not meet these requirements are filtered
    out. These conditions depend on the rotation type (common or
    specific) and the lattice type (crystallogrphic point group and mu, nu)

    Parameters
    ----------------
    quad_int: numpy array
        Integer quadruples.

    cryst_ptgrp: string
        Proper point group in Schoenflies notation.

    sig_type: {'common', 'specific'}

    args[0]: dictionary
    keys: 'nu', 'mu', 'kmax'

    Returns
    -----------
    quad_int: numpy array
        Integer quadruple array that satisfy the above mentioned condition.

    See Also
    -----------
    check_fsig_int

    """
    m = quad_int[0, :]
    u = quad_int[1, :]
    v = quad_int[2, :]
    w = quad_int[3, :]

    tol = 1e-10

    if sig_type == 'common':
        if cryst_ptgrp == 'D3':
            ind1 = np.where((u == v) & (v == w))[0]
            cond1 = ((np.remainder(m[ind1], 2) == 0)
                     & (np.remainder(u[ind1], 2) == 0))
            sig_inds = (m[ind1] ** 2 + 3 * (u[ind1] ** 2)) / 4.
            sig_inds[cond1] = 4 * sig_inds[cond1]
            ind2 = ind1[np.where(abs(sig_inds - sigma) < tol)[0]]

            ind3 = np.where((m == 0) & (u + v + w == 0))[0]
            sig_inds1 = u[ind3] ** 2 + u[ind3] * v[ind3] + v[ind3] ** 2
            ind4 = ind3[np.where(abs(sig_inds1 - sigma) < tol)[0]]

            inds = np.sort(np.concatenate((ind2, ind4)))
            return quad_int[:, inds]

        if cryst_ptgrp == 'D4':
            ind1 = np.where((u == 0) & (v == 0))[0]
            t1 = m[ind1] ** 2 + w[ind1] ** 2
            sig_inds = largest_odd_factor(t1)
            ind2 = ind1[np.where(abs(sig_inds - sigma) < tol)[0]]
            ind3 = np.where((m == 0) & (w == 0))[0]
            t2 = u[ind3] ** 2 + v[ind3] ** 2
            sig_inds1 = largest_odd_factor(t2)
            ind4 = ind3[np.where(abs(sig_inds1 - sigma) < tol)[0]]
            inds = np.sort(np.concatenate((ind2, ind4)))
            return quad_int[:, inds]

        if cryst_ptgrp == 'D6':
            ind1 = np.where((u == 0) & (v == 0))[0]

            t1 = gcd1d_arr((3 * np.ones(np.shape(w[ind1])), w[ind1]))
            t2 = gcd1d_arr((2 * np.ones(np.shape(w[ind1])), m[ind1] + w[ind1]))
            sig_inds = (3 * (m[ind1] ** 2) + w[ind1] ** 2) / (t1 * (t2 ** 2))
            ind2 = ind1[np.where(abs(sig_inds - sigma) < tol)[0]]

            ind3 = np.where((m == 0) & (w == 0))[0]
            t3 = gcd1d_arr((3 * np.ones(np.shape(u[ind3])), u[ind3] + v[ind3]))
            sig_inds1 = (u[ind3] ** 2 - u[ind3] * v[ind3] + v[ind3] ** 2) / t3
            ind4 = ind3[np.where(abs(sig_inds1 - sigma) < tol)[0]]

            inds = np.sort(np.concatenate((ind2, ind4)))
            return quad_int[:, inds]

        if cryst_ptgrp == 'O':
            sig_ind = m ** 2 + u ** 2 + v ** 2 + w ** 2
            return quad_int[:, np.where(sig_ind == sigma)[0]]

    elif sig_type == 'specific':
        nu = args[0]['nu']
        mu = args[0]['mu']
        if cryst_ptgrp == 'D3':
            t1s = np.ones(np.shape(m))
            twos = 2 * t1s
            f1 = gcd1d_arr(([twos, m + u + v + w]))
            f2 = gcd1d_arr(([twos, m + u + v + w, u - v, v - w]))

            t1 = 2 * (u - v) / f1
            cond1 = abs(t1 - np.round(t1)) < 1e-06
            t1 = 2 * (v - w) / f1
            cond2 = abs(t1 - np.round(t1)) < 1e-06
            t1 = 2 * m / f2
            cond3 = abs(t1 - np.round(t1)) < 1e-06
            condfin = (cond1 & cond2 & cond3)

            quad_int = quad_int[:, condfin]
            if np.shape(quad_int)[1] == 0:
                quad_int_out = []
                return quad_int_out

            m = quad_int[0, :]
            u = quad_int[1, :]
            v = quad_int[2, :]
            w = quad_int[3, :]

            t1s = np.ones(np.shape(m))
            twos = 2 * t1s
            f3 = gcd1d_arr(([twos, 2 * (u - v) / f1, 2 * (v - w) / f1]))

            t1 = 2*(mu*v + nu*(u-2*v+w))/f3
            cond1 = abs(t1 - np.round(t1)) < 1e-06
            condfin = cond1

            quad_int = quad_int[:, condfin]
            if np.shape(quad_int)[1] == 0:
                quad_int_out = []
                return quad_int_out

            m = quad_int[0, :]
            u = quad_int[1, :]
            v = quad_int[2, :]
            w = quad_int[3, :]

            f = mu*(m**2) + (mu - 2*nu)*(u**2 + v**2 + w**2) + \
                2*nu*(v*w + w*u + u*v)
            f1 = gcd1d_arr(([twos, m+u+v+w]))
            f2 = gcd1d_arr(([twos, m+u+v+w, u-v, v-w]))
            f3 = gcd1d_arr(([mu*t1s, 2*(u-v)/f1, 2*(v-w)/f1]))
            f4 = gcd1d_arr(([mu*t1s-3*nu*t1s, 2*m/f2, m+u+v+w,
                             2*(mu*v + nu*(u-2*v+w))/f3]))
            sig1 = f/(f1*f2*f3*f4)

            quad_int_out = quad_int[:, sig1 == sigma]
            return quad_int_out

        if cryst_ptgrp == 'D4':
            t1s = np.ones(np.shape(m))
            twos = 2*t1s
            threes = 3*t1s

            f2 = gcd1d_arr((twos, mu*t1s, u + v))
            f3 = gcd1d_arr((twos, nu*t1s, m+w))
            f4 = gcd1d_arr((twos, mu*t1s, u, v, m+w))
            f5 = gcd1d_arr((twos, nu*t1s, m, w, u+v))

            t1 = (mu*t1s)/f2
            cond1 = abs(t1 - np.round(t1)) < 1e-06
            t1 = (nu*t1s)/f3
            cond2 = abs(t1 - np.round(t1)) < 1e-06
            t1 = (u+v)/f4
            cond3 = abs(t1 - np.round(t1)) < 1e-06
            t1 = (m+w)/f5
            cond4 = abs(t1 - np.round(t1)) < 1e-06
            condfin = cond1 & cond2 & cond3 & cond4

            quad_int = quad_int[:, condfin]
            if np.shape(quad_int)[1] == 0:
                quad_int_out = []
                return quad_int_out

            m = quad_int[0, :]
            u = quad_int[1, :]
            v = quad_int[2, :]
            w = quad_int[3, :]

            t1s = np.ones(np.shape(m))
            twos = 2*t1s
            fours = 4*t1s
            f = mu*(m**2 + w**2) + nu*(u**2 + v**2)
            f1 = gcd1d_arr((fours, ((mu+nu)**2)*t1s,
                            m**2 + u**2 + v**2 + w**2))
            f2 = gcd1d_arr((twos, mu*t1s, u+v))
            f3 = gcd1d_arr((twos, nu*t1s, m+w))
            f4 = gcd1d_arr((twos, mu*t1s, u, v, m+w))
            f5 = gcd1d_arr((twos, nu*t1s, m, w, u+v))
            f6 = gcd1d_arr(((mu*t1s)/f2, u, v, (u+v)/f4))
            f7 = gcd1d_arr(((nu*t1s)/f3, m, w, (m+w)/f5))

            sig1 = f/(f1*f2*f3*f4*f5*f6*f7)

            quad_int_out = quad_int[:, sig1 == sigma]
            return quad_int_out

        if cryst_ptgrp == 'D6':
            twos = 2*np.ones(np.shape(m))
            threes = 3*np.ones(np.shape(m))

            f1 = gcd1d_arr(((u, v, m+w)))
            f2 = gcd1d_arr(((threes, u+v, w)))

            cond1 = abs(twos / f1 - np.round(twos / f1)) < 1e-06
            cond2 = abs(2*w / (f1*f2) - np.round(2*w / (f1*f2))) < 1e-06
            cond3 = abs(3*u / (f1*f2) - np.round(3*u / (f1*f2))) < 1e-06
            cond4 = abs((u+v) / f1 - np.round((u+v) / f1)) < 1e-06
            condfin = cond1 & cond2 & cond3 & cond4

            quad_int = quad_int[:, condfin]
            if np.shape(quad_int)[1] == 0:
                quad_int_out = []
                return quad_int_out

            m = quad_int[0, :]
            u = quad_int[1, :]
            v = quad_int[2, :]
            w = quad_int[3, :]

            twos = 2*np.ones(np.shape(m))
            nus = nu*np.ones(np.shape(m))
            f1 = gcd1d_arr(((u, v, m+w)))

            f3 = gcd1d_arr(((twos / f1, nus, m+w)))
            cond0 = abs(nus / f3 - np.round(nus / f3)) < 1e-06

            quad_int = quad_int[:, cond0]
            if np.shape(quad_int)[1] == 0:
                quad_int_out = []
                return quad_int_out

            m = quad_int[0, :]
            u = quad_int[1, :]
            v = quad_int[2, :]
            w = quad_int[3, :]

            twos = 2*np.ones(np.shape(m))
            threes = 3*np.ones(np.shape(m))
            nus = nu*np.ones(np.shape(m))
            mus = mu*np.ones(np.shape(m))

            f = mu*(3*(m**2) + w**2) + nu*(u**2 - u*v + v**2)
            f1 = gcd1d_arr((u, v, m+w))
            f2 = gcd1d_arr((threes, u+v, w))
            f3 = gcd1d_arr((twos / f1, nus, m+w))
            f4 = gcd1d_arr((nus / f3, 2*w / (f1*f2), m+w))
            f5 = gcd1d_arr((mus, 3*u / (f1*f2), (u+v) / f1))
            sig1 = f / (f1*f2*f3*f4*f5)

            quad_int_out = quad_int[:, sig1 == sigma]
            return quad_int_out
    else:
        raise Exception('sig_type: wrong input type')
# -----------------------------------------------------------------------------------------------------------

def gcd1d_arr(arr_tup):
    """
    A tuple of one-D arrays are passed with equal size and the gcd of
    their rows is computed

    Parameters
    ---------------
    arr_tup: tuple
        one-D arrays of integers of equal size.

    Returns
    -----------
    GCD of rows of 1D arrays of integers
    """
    asz = len(arr_tup)
    gc1 = arr_tup[0].astype(int)
    for ct1 in range(asz - 1):
        gc1 = np.copy(np.column_stack((gc1,
                                       arr_tup[ct1+1].astype(int))))

    if gc1.ndim == 1:
        gc1 = np.reshape(gc1, (1, np.size(gc1)))

    t1 = int_man.gcd_array(gc1, 'rows')
    return np.reshape(t1, (np.size(t1), ))
# -----------------------------------------------------------------------------------------------------------

def compute_tmat(quad_int, tau, lat_type):
    """
    The transformation matrix (r_g1tog2_g1) corresponding to the integer
    quadruple is computed. The matrix elements depend on m, U, V, W and the
    crystallographic point group and tau = (nu/mu)

    Parameters
    ------------------
    quad_int: numpy array
        Integer quadruples.

    tau: float
        :math:`\\frac{\\nu}{\\mu}`

    lat_type: Lattice class
        Attributes of the underlying lattice

    Returns
    ----------
    g : numpy array
        dimension = 3, n x 3 x  3 transformation matrices
    """
    m = quad_int[0, :].astype(float)
    u = quad_int[1, :].astype(float)
    v = quad_int[2, :].astype(float)
    w = quad_int[3, :].astype(float)

    sz = np.size(m)
    g = np.zeros((sz, 3, 3))

    cryst_ptgrp = proper_ptgrp(lat_type.cryst_ptgrp)

    if cryst_ptgrp == 'D3':
        k1 = (1 - 2 * tau)
        s = (m ** 2 + k1 * (u ** 2 + v ** 2 + w ** 2) +
             (2 * tau) * (v * w + w * u + u * v))

        g[:, 0, 0] = (m ** 2 + k1 * (u ** 2 - v ** 2 - w ** 2) +
                      (2 * tau) * (m * v - m * w - v * w)) / s
        g[:, 1, 1] = (m ** 2 + k1 * (-u ** 2 + v ** 2 - w ** 2)
                      + (2 * tau) * (m * w - m * u - w * u)) / s
        g[:, 2, 2] = (m ** 2 + k1 * (-u ** 2 - v ** 2 + w ** 2)
                      + (2 * tau) * (m * u - m * v - u * v)) / s

        g[:, 0, 1] = 2*(tau*u*(w + u - m) - (1-tau)*(m*w) + k1*(u*v)) / s
        g[:, 1, 0] = 2*(tau*v*(v + w + m) + (1-tau)*(m*w) + k1*(u*v)) / s

        g[:, 0, 2] = 2*(tau*u*(u + v + m) + (1-tau)*(m*v) + k1*(w*u)) / s
        g[:, 2, 0] = 2*(tau*w*(v + w - m) - (1-tau)*(m*v) + k1*(w*u)) / s

        g[:, 1, 2] = 2*(tau*v*(u + v - m) - (1-tau)*(m*u) + k1*(v*w)) / s
        g[:, 2, 1] = 2*(tau*w*(w + u + m) + (1-tau)*(m*u) + k1*(v*w)) / s

    if cryst_ptgrp == 'D6':
        s = (3*m ** 2 + w ** 2 + tau*(u ** 2 - u*v + v ** 2))/3
        g[:, 0, 0] = (3*m ** 2 + 2*m*w - w ** 2 + tau*(u ** 2 - v ** 2))/(3*s)
        g[:, 1, 1] = (3*m ** 2 - 2*m*w - w ** 2 - tau*(u ** 2 - v ** 2))/(3*s)
        g[:, 2, 2] = (3*m ** 2 + w ** 2 - tau*(u ** 2 - u*v + v ** 2)) / (3*s)

        g[:, 0, 1] = (tau*u*(2*v - u) - 4*m*w) / (3*s)
        g[:, 1, 0] = (tau*v*(2*u - v) + 4*m*w) / (3*s)

        g[:, 0, 2] = (2*(u*w + m*(2*v - u))) / (3*s)
        g[:, 2, 0] = (tau*(w*(2*u-v) - 3*m*v)) / (3*s)

        g[:, 1, 2] = (2*(v*w - m*(2*u - v))) / (3*s)
        g[:, 2, 1] = (tau*(w*(2*v-u) + 3*m*u)) / (3*s)

    if cryst_ptgrp == 'D4':

        s = m ** 2 + w ** 2 + tau*(u ** 2 + v ** 2)

        g[:, 0, 0] = (m ** 2 - w ** 2 + tau*(u ** 2 - v ** 2)) / s
        g[:, 1, 1] = (m ** 2 - w ** 2 - tau*(u ** 2 - v ** 2)) / s
        g[:, 2, 2] = (m ** 2 + w ** 2 - tau*(u ** 2 + v ** 2)) / s

        g[:, 0, 1] = (2*(tau*u*v - m*w)) / s
        g[:, 1, 0] = (2*(tau*u*v + m*w)) / s

        g[:, 0, 2] = (2*(u*w + m*v)) / s
        g[:, 2, 0] = tau*(2*(u*w - m*v)) / s

        g[:, 1, 2] = (2*(v*w - m*u)) / s
        g[:, 2, 1] = tau*(2*(v*w + m*u)) / s

    if cryst_ptgrp == 'O':

        s = m ** 2 + w ** 2 + u ** 2 + v ** 2

        g[:, 0, 0] = (m ** 2 - w ** 2 + u ** 2 - v ** 2) / s
        g[:, 1, 1] = (m ** 2 - w ** 2 - u ** 2 + v ** 2) / s
        g[:, 2, 2] = (m ** 2 + w ** 2 - u ** 2 - v ** 2) / s

        g[:, 0, 1] = 2*(u*v - m*w) / s
        g[:, 1, 0] = 2*(u*v + m*w) / s

        g[:, 0, 2] = 2*(u*w + m*v) / s
        g[:, 2, 0] = 2*(u*w - m*v) / s

        g[:, 1, 2] = 2*(v*w - m*u) / s
        g[:, 2, 1] = 2*(v*w + m*u) / s

        if lat_type.pearson == 'cF' or lat_type.pearson == 'cI':
            l_p_po = lat_type.l_p_po
            l_go_g = np.linalg.inv(l_p_po)
            for i in range(np.shape(g)[0]):
                g[i, :, :] = np.dot(np.dot(l_go_g, g[i, :, :]), l_p_po)

    if sz == 1:
        g = g.reshape((3, 3))

    return g
# -----------------------------------------------------------------------------------------------------------

def disorient_sigmarots(r_g1tog2_g1, l_p_po, cryst_ptgrp):
    """
    The disorientation corresponding to each rotation matrix is computed
    and the unique set is returned

    Parameters
    ----------------
    r_g1tog2_g1: numpy array (n x 3 x 3)
        Transformation matrices in g1 reference frame

    l_p_po: numpy array
        The primitive basis vectors of the underlying lattice in the orthogonal
        reference frame.

    cryst_ptgrp: string
        Proper point group in Schoenflies notation

    Returns
    ----------
    rots_g1tog2_g1: numpy array (n x 3 x 3)
        Transformation matrices in g1 reference frame in the fundamental zone
    """
    if r_g1tog2_g1.ndim == 2:
        r_g1tog2_g1 = np.reshape(r_g1tog2_g1, (1, 3, 3))

    l_go_g = np.linalg.inv(l_p_po)
    msz = np.shape(r_g1tog2_g1)[0]

    r_go1togo2_go1 = np.zeros(np.shape(r_g1tog2_g1))
    r_go1togo2_go1[:] = np.NaN

    for i in range(msz):
        r_go1togo2_go1[i, :, :] = np.dot(
            np.dot(l_p_po, r_g1tog2_g1[i, :, :]), l_go_g)

    q_go1togo2_go1 = trans.mat2quat(r_go1togo2_go1)
    qfz_go1togo2_go1 = mis_fz.misorient_fz(q_go1togo2_go1, cryst_ptgrp)

    # qt1 = quat.double(qfz_go1togo2_go1)
    qt1 = np.array(qfz_go1togo2_go1)
    qt1 = qt1.transpose()
    [t1, ia] = trans.unique_rows_tol(qt1, 1e-06, True)

    # Change rotations to the fundamental zone
    rots_g1tog2_g1 = np.zeros((np.size(ia), 3, 3))
    rots_g1tog2_g1[:] = np.NaN
    ct2 = 0
    for ct1 in range(np.size(ia)):
        if abs(abs(qfz_go1togo2_go1[:, ia[ct1]][0])-1) > 1e-10:
            mat1 = trans.quat2mat(qfz_go1togo2_go1[:, ia[ct1]])
            rots_g1tog2_g1[ct2, :, :] = np.dot(np.dot(l_go_g, mat1), l_p_po)
            ct2 += 1

    if ct2 == 0:
        return np.zeros(0)
    else:
        return rots_g1tog2_g1
# -----------------------------------------------------------------------------------------------------------

def check_sigma_rots(r_g1tog2_g1, sigma):
    """
    The sigma transformation matrix has the property that sigma is the
    smallest integer such that sigma*T is an integer matrix. This condition
    is checked and the numerator and denominatr(sigma) matrices are
    returned

    Parameters
    ----------------
    r_g1tog2_g1: numpy array (n x 3 x 3)
        Transformation matrices in g1 reference frame

    sigma: float
        sigma number

    Returns
    ----------
    {'N': rots_n, 'D': rots_d}: dictionary
    rots_n: numpy array
        numerator matrices n x 3 x3
    rots_d: numpy array
        denominator matrices n x 3 x3
    """
    rots_n = np.zeros(np.shape(r_g1tog2_g1))
    rots_d = np.zeros(np.shape(r_g1tog2_g1))

    msz = np.shape(r_g1tog2_g1)[0]
    for i in range(msz):
        tmat = r_g1tog2_g1[i, :, :]
        mult1 = int_man.int_mult(tmat, 1e-06)
        if abs(mult1[0]-sigma) > 1e-06:
            t_check = 0
            tol_arr = np.array([1e-05, 1e-06, 1e-07,
                                1e-08, 1e-09, 1e-10, 1e-4])
            for tols in tol_arr:
                mult1 = int_man.int_mult(tmat, tols)
                if abs(mult1[0]-sigma) < 1e-04:
                    t_check = 1
                    break
            if t_check == 0:
                raise Exception('Not a sigma rotation')

        rots_n[i, :, :] = mult1[1]
        rots_d[i, :, :] = int(np.round(mult1[0]))*np.ones(np.shape(tmat))

    if np.size(np.where(rots_d != sigma)[0]) > 0:
        raise Exception('Not a sigma rotation')

    return {'N': rots_n, 'D': rots_d}
# -----------------------------------------------------------------------------------------------------------

def csl_rotations(sigma, sig_type, lat_type):
    """
    The function computes the CSL rotation matrices r_g1tog2_g1 corresponding
    to a give sigma and lattice

    Parameters
    ----------
    sigma : int
        Sigma corresponding to the transformation matrix

    sig_type: {'common', 'specific'}
        If the sigma generating function depends on the lattice type, then
        sig_type is 'specific', otherwise it is 'common'

    lat_type: Lattice class
        Attributes of the underlying lattice

    Returns
    -------
    sig_rots: dictionary
    keys: 'N', 'D'
    sig_rots['N'], sig_rots['D']: Numerator and Integer matrices
        The transformation matrix is N/D in the g1 reference frame
        (i.e. r_g1tog2_g1)

    Notes
    -------
    The following steps are considered to obtain the sigma rotation:

    * compute_inp_params: computes tau and kmax that fixes the range of
      integer qudruples sampled
    * mesh_muvw: Creates the integer quadruples that depend on sigma,
      tau, kmax, crystallographic point group
    * eliminate_idrots: Eliminates Identity rotations
    * If specific rotations are desired:
        -  mesh_muvw_fz: Restricts quadruples to fundamental zone
        -  check_fsig_int: Filters out quadruples that do not meet the
           condition specified in this function
    * sigtype_muvw: Filters out quadruple combinations depending on
      the type of sigma rotation
    * eliminate_mults: Eliminates integer quadruples that are same except for
      a scaling factor
    * check_sigma: Returns integer quadruples that result in the sigma rotation
    * compute_tmat: Computes the transformation matrix from the integer
      quadruple
    * disorient_sigmarots: Converts all the transformations to the fundamental
      zone of the corresponding crystallogrphic point group
    * check_sigma_rots: Checks that the transformation matrix is a sigma
      rotation and returns them as numerator and denominator matrices
    """

    cryst_ptgrp = proper_ptgrp(lat_type.cryst_ptgrp)
    lat_type.cryst_ptgrp = cryst_ptgrp

    if cryst_ptgrp in ['T', 'Td', 'Th', 'O', 'Oh']:
        if np.remainder(sigma, 2) == 0:
            return {'N': np.empty(0), 'D': np.empty(0)}

    # Define Parameters
    [tau, kmax] = compute_inp_params(lat_type, sig_type)
    [nu, mu] = int_man.rat(tau)

    if sig_type == 'specific':
        lat_args = {}
        lat_args['mu'] = mu[0][0]
        lat_args['nu'] = nu[0][0]
        lat_args['kmax'] = kmax[0][0]
        # Create Integer Quadruples
        quad_int = mesh_muvw(cryst_ptgrp, sigma, sig_type, lat_args)
        # Restrict to Fundamental Zone
        quad_int = mesh_muvw_fz(quad_int, cryst_ptgrp, sig_type, lat_args)
        # Eliminate Identity Rotations
        quad_int = eliminate_idrots(quad_int)
        # Check $\frac{F}{\Sigma}$ is an iteger and a divisor of kmax
        quad_int = check_fsig_int(quad_int, cryst_ptgrp, sigma, lat_args)
    else:
        quad_int = mesh_muvw(cryst_ptgrp, sigma, sig_type)
        # Eliminate Identity Rotations
        quad_int = eliminate_idrots(quad_int)

    # Keep only 'specific' or 'common' rotations
    quad_int = sigtype_muvw(quad_int, cryst_ptgrp, sig_type)

    # Keep only those quadruples such that gcd(m, U, V, W) = 1
    quad_int = eliminate_mults(quad_int)

    # Compute $\Sigma$ and check with input $\Sigma$
    if sig_type == 'common':
        quad_int = check_sigma(quad_int, sigma, cryst_ptgrp, sig_type)
    if sig_type == 'specific':
        quad_int = check_sigma(quad_int, sigma, cryst_ptgrp,
                               sig_type, lat_args)

    if np.size(quad_int) == 0:
        return {'N': np.empty(0), 'D': np.empty(0)}

    # Compute rotation matrices in G1 lattice
    r_g1tog2_g1 = compute_tmat(quad_int, tau, lat_type)
    l_p_po = lat_type.l_p_po

    # Convert to disorientations to keep the unique rotations
    r_g1tog2_g1 = disorient_sigmarots(r_g1tog2_g1, l_p_po, cryst_ptgrp)

    if np.size(r_g1tog2_g1) == 0:
        return {'N': np.empty(0), 'D': np.empty(0)}
    else:
        # Check that r_g1tog2_g1 are rational with lcm of denominator matrices
        # equal to $\Sigman$
        sig_rots = check_sigma_rots(r_g1tog2_g1, sigma)
        return sig_rots
# -----------------------------------------------------------------------------------------------------------
