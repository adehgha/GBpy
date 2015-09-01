# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
from fractions import gcd
import integer_manipulations as int_man
import find_csl_dsc as fcd
import sys
import os
from tools import Col, lll_reduction, smith_nf, extgcd
# -----------------------------------------------------------------------------------------------------------


def check_int_mats(l1, l2):
    """
    The function checks whether or not the elements of
    $A_{left}^{-1}$*B are integers

    Parameters
    ----------
    l1, l2 : numpy arrays of same dimensions

    Returns
    -------
    cond: numpy boolean array
    """

    l1_linv = np.dot(np.linalg.inv(np.dot(l1.T, l1)), l1.T)
    cond = int_man.int_check(np.dot(l1_linv, l2), 10)
    return cond
# -----------------------------------------------------------------------------------------------------------


def check_2d_csl(l_pl1_g1, l_pl2_g1, l_csl_g1):
    """
    The function checks whether or not the CSL basis may be expressed
    as a linear integer combination of the plane bases of planes 1 and 2

    Parameters
    ----------------
    l_pl1_g1, l_pl2_g1: numpy arrays of basis vectors for plane 1 and 2
    in the g1 reference frame
    """

    ct = 1
    str1 = 'Check ' + str(ct) + ' : --> '
    if check_int_mats(l_pl1_g1, l_csl_g1).all():
        str2 = 'YES\n'
        color = 'yel'
    else:
        raise Exception('The 2D CSL does not contain base1.')

    print str1
    txt = Col()
    txt.c_prnt(str2, color)
    # --------------------------
    ct += 1
    str1 = 'Check ' + str(ct) + ' : --> '
    if check_int_mats(l_pl2_g1, l_csl_g1).all():
        str2 = 'YES\n'
        color = 'yel'
    else:
        raise Exception('The 2D CSL does not contain base2.')

    print str1
    txt = Col()
    txt.c_prnt(str2, color)
# -----------------------------------------------------------------------------------------------------------


def lbi_dioph_soln(a, b, c):
    """
    Computes the diophantaine solution for the equation ax + by = c
    """
    k = abs(gcd(a, b))
    if np.fmod(c, k) != 0.0:
        raise Exception('No Solutions Exist')

    mult = c / k

    # b1 = extended_euclid(a, b) ####<----------------- must be verified
    b1 = np.zeros(2)
    b1[0], b1[1], _ = extgcd(a, b)

    x1 = mult * b1[0]
    x2 = mult * b1[1]

    if a*x1 + b*x2 - c != 0:
        raise Exception('Something wrong with the Algorithm')

    tk = -(b-a)/(a**2 + b**2)
    tk_mat = np.array([[np.ceil(tk)], [np.floor(tk)]])
    x_sol = np.array([[x1], [x1]]) + tk_mat * b
    y_sol = np.array([[x2], [x2]]) - tk_mat * a
    sol_mag = np.power(x_sol, 2) + np.power(y_sol, 2)
    ind = np.where(sol_mag == min(sol_mag))[0]
    int_soln = [x_sol[ind[0]], y_sol[ind[0]]]
    return int_soln
# -----------------------------------------------------------------------------------------------------------


def compute_basis_vec(d_eq):
    """
    The function computes y1, y2, y3 such that h\*y1 + k\*y2 + l\*y3 = 0
    and modulus of y1 is a minimum

    Parameters
    -----------------
    d_eq: numpy array or list of size 3 and dimension 1
        h = d_eq[0], k = d_eq[1], l = d_eq[2]

    Returns
    ------------
    np.array([y1, y2, y3])

    """
    hp = d_eq[0]
    kp = d_eq[1]
    lp = d_eq[2]
    # Find the minimum y1 such that y2 and y3 are solutions of the equation.
    # kp*y2 + lp*y3 = -hp*y1 (Diaphontane Equation).
    # Solutions exist if gcd(kp,lp) is a multiple of hp*y1
    cond = 0
    y1 = 1
    while cond == 0:
        if np.fmod(hp * y1, gcd(kp, lp)) == 0:
            cond = 1
        else:
            y1 += 1

    # Diophantine Equation: ax + by = c
    # To solve: f = kp*x + lp*y + hp*m = 0;
    avar = kp
    bvar = lp
    cvar = -hp * y1
    int_soln = lbi_dioph_soln(avar, bvar, cvar)
    y2 = int_soln[0]
    y3 = int_soln[1]
    if (kp*y2 + lp*y3 + hp*y1) != 0:
        raise Exception('Error with Diophantine solution')

    return np.array([y1, y2, y3])
# -----------------------------------------------------------------------------------------------------------


def bp_basis(miller_ind):
    """
    The function computes the primitve basis of the plane if the
    boundary plane indices are specified

    Parameters
    ---------------
    miller_ind: numpy array
        Miller indices of the plane (h k l)

    Returns
    -----------
    l_pl_g1: numpy array
        The primitive basis of the plane in 'g1' reference frame
    """
    miller_ind = int_man.int_finder(miller_ind)
    h = miller_ind[0]
    k = miller_ind[1]
    l = miller_ind[2]

    if h == 0 and k == 0 and l == 0:
        raise Exception('hkl indices cannot all be zero')
    else:
        if h != 0 and k != 0 and l != 0:
            gc_f1_p = gcd(k, l)
            bv1_g1 = np.array([[0], [-l / gc_f1_p], [k / gc_f1_p]])
            bv2_g1 = compute_basis_vec([h, k, l])
        else:
                if h == 0:
                    if k == 0:
                        bv1_g1 = np.array([[1], [0], [0]])
                        bv2_g1 = np.array([[0], [1], [0]])
                    elif l == 0:
                        bv1_g1 = np.array([[0], [0], [1]])
                        bv2_g1 = np.array([[1], [0], [0]])
                    else:
                        gc_f1_p = gcd(k, l)
                        bv1_g1 = np.array([[0], [-l / gc_f1_p],
                                           [k / gc_f1_p]])
                        bv2_g1 = np.array([[1], [-l / gc_f1_p],
                                           [k / gc_f1_p]])
                else:
                    if k == 0:
                        if l == 0:
                            bv1_g1 = np.array([[0], [1], [0]])
                            bv2_g1 = np.array([[0], [0], [1]])
                        else:
                            gc_f1_p = gcd(h, l)
                            bv1_g1 = np.array([[-l / gc_f1_p], [0],
                                               [h / gc_f1_p]])
                            bv2_g1 = np.array([[-l / gc_f1_p], [1],
                                               [h / gc_f1_p]])
                    else:
                        if l == 0:
                            gc_f1_p = gcd(h, k)
                            bv1_g1 = np.array([[-k / gc_f1_p],
                                               [h / gc_f1_p], [0]])
                            bv2_g1 = np.array([[-k / gc_f1_p],
                                               [h / gc_f1_p], [1]])

    #  The reduced basis vectors for the plane
    l_pl_g1 = lll_reduction(np.column_stack([bv1_g1, bv2_g1]))
    return l_pl_g1
# -----------------------------------------------------------------------------------------------------------


def pl_density(l_pl_g1, l_g1_go1):
    """
    For a given two-dimensional plane basis, the planar density is
    computed

    Parameters
    ---------------
    l_pl_g1: numpy array

    l_g1_go1: numpy array
        Basis vectors of the underlying lattice with respect to the
        orthogonal reference frame 'go1'

    Returns
    ----------
    pd: float
        Planar density = (1/area covered by plane basis)
    """
    l_pl_go1 = np.dot(l_g1_go1, l_pl_g1)
    planar_basis_area = np.linalg.norm(np.cross(l_pl_go1[:, 0],
                                                l_pl_go1[:, 1]))
    pd = 1.0/planar_basis_area
    return pd
# -----------------------------------------------------------------------------------------------------------


def csl_finder_2d(l_pl1_g1, l_pl2_g1):
    """
    Given two plane bases, the 2D CSL bases are obtined by utilizing the
    smith normal form of the transformation between the two bases

    Parameters
    ---------------
    l_pl1_g1, l_pl2_g1: numpy array
        Basis vectors of planes 1 and 2 expressed in g1 reference frame

    Returns
    ---------------
    l_2d_csl_g1: numpy array
        The basis vectors of the 2D CSL expressed in g1 reference frame
    """

    l1 = np.array(l_pl1_g1)
    l2 = np.array(l_pl2_g1)

    # L2 = L1*T (T := T_Pl1ToPl2_Pl1)
    l1_linv = np.dot(np.linalg.inv(np.dot(l1.T, l1)), l1.T)
    t = np.dot(l1_linv, l2)

    # # Obtain N1 and N2
    n1, a = int_man.int_mult(t)

    if int_man.int_check(a, 12).all:
        a = np.around(a)

    # A = U*E*V'; E is a diagonal matrix
    u, e, v = smith_nf(a)

    l2_n = np.dot(l2, np.linalg.inv(v.T))
    l_2d_csl_g1 = fcd.csl_elem_div_thm_l2(e/n1, l2_n)

    if int_man.int_check(l_2d_csl_g1, 10).all():
        l_2d_csl_g1 = np.around(l_2d_csl_g1)
    else:
        raise Exception('Wrong CSL computation')

    l_2d_csl_g1 = lll_reduction(l_2d_csl_g1[:, 0:2])

    return l_2d_csl_g1
# -----------------------------------------------------------------------------------------------------------


def gb_2d_csl(inds, t_mat, l_g_go,
              inds_type='miller_index', mat_ref='g1'):
    """
    For a given boundary plane normal 'bp1_g1' and the misorientation
    matrix 't_g1tog2_g1', the two-dimensional CSL lattice is computed

    Parameters
    ------------------
    inds: numpy array
        The boundary plane indices

    inds_type: string
        {'miller_index', 'normal_go', 'normal_g'}

    t_mat: numpy array
        Transformation matrix from g1 to g2 in 'mat_ref' reference frame

    mat_ref: string
        {'go1', 'g1'}

    lattice: Lattice class
    Attributes of the underlying lattice

    Returns
    -----------
    l_2d_csl_g1, l_pl1_g1, l_pl2_g1: numpy arrays
        ``l_2d_csl_g1`` is the 2d CSL in g1 ref frame.\v
        ``l_pl1_g1`` is the plane 1 basis in g1 ref frame.\v
        ``l_pl2_g1`` is the plane 2 basis in g1 ref frame.\v
    """

    l_g1_go1 = l_g_go
    l_go1_g1 = np.linalg.inv(l_g1_go1)
    l_rg1_go1 = fcd.reciprocal_mat(l_g1_go1)
    l_go1_rg1 = np.linalg.inv(l_rg1_go1)

    if inds_type == 'normal_go':
        bp1_go1 = inds
        miller1_ind = int_man.int_finder(np.dot(l_go1_rg1, bp1_go1))
    elif inds_type == 'miller_index':
        miller1_ind = inds
    elif inds_type == 'normal_g':
        bp1_g1 = inds
        l_g1_rg1 = np.dot(l_go1_rg1, l_g1_go1)
        miller1_ind = int_man.int_finder(np.dot(l_g1_rg1, bp1_g1))
    else:
        raise Exception('Wrong index type')

    if mat_ref == 'go1':
        t_g1tog2_g1 = np.dot(l_go1_g1, np.dot(t_mat, l_g1_go1))
    elif mat_ref == 'g1':
        t_g1tog2_g1 = t_mat
    else:
        raise Exception('Wrong reference axis type')

    bp1_go1 = int_man.int_finder(np.dot(l_rg1_go1, miller1_ind))
    l_g2_g1 = t_g1tog2_g1
    l_g2_go1 = np.dot(l_g1_go1, l_g2_g1)
    l_rg2_go1 = fcd.reciprocal_mat(l_g2_go1)
    l_go1_rg2 = np.linalg.inv(l_rg2_go1)
    # bp2_g2 = int_man.int_finder(np.dot(-l_go1_g2, bp1_go1))
    miller2_ind = int_man.int_finder(np.dot(-l_go1_rg2, bp1_go1))

    l_pl1_g1 = bp_basis(miller1_ind)
    l_pl1_g1 = reduce_go1_mat(l_pl1_g1, l_g1_go1)
    l_pl2_g2 = bp_basis(miller2_ind)
    l_pl2_g2 = reduce_go1_mat(l_pl2_g2, l_g1_go1)
    l_pl2_g1 = np.dot(l_g2_g1, l_pl2_g2)

    l_2d_csl_g1 = csl_finder_2d(l_pl1_g1, l_pl2_g1)
    l_2d_csl_g1 = reduce_go1_mat(l_2d_csl_g1, l_g1_go1)

    return l_2d_csl_g1, l_pl1_g1, l_pl2_g1
# -----------------------------------------------------------------------------------------------------------

def reduce_go1_mat(l_pl1_g1, l_g1_go1):
    l_go1_g1 = np.linalg.inv(l_g1_go1)
    l_pl1_go1 = np.dot(l_g1_go1, l_pl1_g1);

    tmat1  = int_man.int_finder(l_pl1_go1)
    t_ind = np.where(abs(tmat1) == abs(tmat1).max())
    t_ind_1 = t_ind[0][0]; t_ind_2 = t_ind[1][0];
    Dnum = tmat1[t_ind_1, t_ind_2] / l_pl1_go1[t_ind_1, t_ind_2]

    l_pl1_go1 = lll_reduction(tmat1)/Dnum
    l_pl1_g1 = np.dot(l_go1_g1, l_pl1_go1)
    return l_pl1_g1

def bicryst_planar_den(inds, t_mat, l_g_go, inds_type='miller_index',
                       mat_ref='go1'):
    """
    The function computes the planar densities of the planes
    1 and 2 and the two-dimensional CSL

    Parameters
    ---------------
    inds: numpy array
        The boundary plane indices.

    inds_type: string
        {'miller_index', 'normal_go', 'normal_g'}

    t_mat: numpy array
        Transformation matrix from g1 to g2 in go1 (or g1) reference frame.

    mat_ref: string
        {'go1', 'g1'}

    lattice: Lattice class
        Attributes of the underlying lattice.

    Returns
    -----------
    pl_den_pl1, pl_den_pl2: numpy array
        The planar density of planes 1 and 2.

    pl_den_csl: numpy array
        The planare density of the two-dimensional CSL.
    """
    l_g1_go1 = l_g_go
    l_rg1_go1 = fcd.reciprocal_mat(l_g1_go1)
    l_go1_rg1 = np.linalg.inv(l_rg1_go1)

    if inds_type == 'normal_go':
        bp1_go1 = inds
        miller1_inds = int_man.int_finder(np.dot(l_go1_rg1, bp1_go1))
    elif inds_type == 'miller_index':
        miller1_inds = inds
    elif inds_type == 'normal_g':
        bp1_g1 = inds
        l_g1_rg1 = np.dot(l_go1_rg1, l_g1_go1)
        miller1_inds = int_man.int_finder(np.dot(l_g1_rg1, bp1_g1))
    else:
        raise Exception('Wrong index type')

    if mat_ref == 'go1':
        l_2d_csl_g1, l_pl1_g1, l_pl2_g1 = gb_2d_csl(miller1_inds,
                                                    t_mat, l_g_go,
                                                    'miller_index', 'go1')
    elif mat_ref == 'g1':
        l_2d_csl_g1, l_pl1_g1, l_pl2_g1 = gb_2d_csl(miller1_inds,
                                                    t_mat, l_g_go,
                                                    'miller_index', 'g1')
    else:
        raise Exception('Wrong reference axis type')

    check_2d_csl(l_pl1_g1, l_pl2_g1, l_2d_csl_g1)

    pl_den_pl1 = pl_density(l_pl1_g1, l_g1_go1)
    pl_den_pl2 = pl_density(l_pl2_g1, l_g1_go1)
    pl_den_csl = pl_density(l_2d_csl_g1, l_g1_go1)

    return pl_den_pl1, pl_den_pl2, pl_den_csl
# -----------------------------------------------------------------------------------------------------------

