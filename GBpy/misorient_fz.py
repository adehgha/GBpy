# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import sys
import pickle
import os
import quaternion as quat
# -----------------------------------------------------------------------------------------------------------


def check_cond(g, cryst_ptgrp, tol):
    """
    Parameters
    ----------------
    g: quaternion object
        Misorientation

    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    tol: float
        Tolerance for the misorientation to belong in the fundamental zone

    Returns
    ------------
    True or False: Boolean
        Depending on whether or not the misorientation is a disorientation
    """

    q0 = quat.getq0(g)
    q1 = quat.getq1(g)
    q2 = quat.getq2(g)
    q3 = quat.getq3(g)

    if cryst_ptgrp == 'D3' or cryst_ptgrp == 'D3d':
        cond1 = q3 > -tol
        cond2 = q3 - q0/np.sqrt(3) < tol
        cond3 = q1 - q0 < tol
        cond4 = q1 + np.sqrt(3)*q2 > -tol
        cond5 = q1 - np.sqrt(3)*q2 > -tol
        cond6 = q1 > -tol
        cond7 = True
        cond8 = True

        if q0 == q1 and q3 >= 0:
            cond7 = (q2 >= 0)
        if q3 == 0:
            cond8 = (q2 >= 0)

        if (cond1 and cond2 and cond3 and cond4 and cond5 and cond6
            and cond7 and cond8):
            return True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if cryst_ptgrp == 'D4' or cryst_ptgrp == 'D4h':
        cond1 = q3 > -tol
        cond2 = (q3 - (np.sqrt(2)-1)*q0) < tol
        cond3 = (q1 - q0) < tol
        cond4 = (q1 - q2) > -tol
        cond5 = q2 > -tol
        cond6 = (q1 + q2 - np.sqrt(2)*q0) < tol
        cond7 = True
        cond8 = True
        cond9 = True
        if q0 == q1:
            cond7 = q2 >= q3
        if np.sqrt(2)*q0 == q1 + q2:
            cond8 = q1 - q2 >= np.sqrt(2)*q3
        if q0 == (np.sqrt(2) + 1)*q3:
            cond9 = q1 >= (np.sqrt(2) + 1)*q2
        if (cond1 and cond2 and cond3 and cond4 and cond5 and cond6 and cond7
            and cond8 and cond9):
            return True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if cryst_ptgrp == 'D6' or cryst_ptgrp == 'D6h':
        cond1 = q3 > -tol
        cond2 = (q3 - (2-np.sqrt(3))*q0) < tol
        cond3 = (q1 - q0) < tol
        cond4 = (q1 - (np.sqrt(3))*q2) > -tol
        cond5 = q2 > -tol
        cond6 = (np.sqrt(3)*q1 + q2 - 2*q0) < tol
        cond7 = True
        cond8 = True
        cond9 = True
        if q0 == q1:
            cond7 = (q2 >= q3)
        if 2*q0 - (np.sqrt(3)*q1 + q2) == 0:
            cond8 = q1 >= (2 + np.sqrt(3))*q2
        if q0 - (2 + np.sqrt(3))*q3 == 0:
            cond9 = q1 - np.sqrt(3)*q2 >= q3
        if (cond1 and cond2 and cond3 and cond4 and cond5 and cond6 and cond7
            and cond8 and cond9):
            return True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if cryst_ptgrp == 'O' or cryst_ptgrp == 'Oh':
        cond1 = q1 - q2 > -tol
        cond2 = q2 - q3 > -tol
        cond3 = q3 > -tol
        cond4 = q0 - (np.sqrt(2)+1)*q1 > -tol
        cond5 = q0 - (q1 + q2 + q3) > -tol
        cond6 = True
        if abs(q0 - (np.sqrt(2)+1)*q1) <= 0:
            cond6 = q2 - (np.sqrt(2)+1)*q3 <= tol
        if cond1 and cond2 and cond3 and cond4 and cond5 and cond6:
            return True
# ------------------------------------------------------------------------------------------------------


def misorient_fz(misquats, cryst_ptgrp, tol=1e-12):
    """
    The function takes as input the misorientations and the corresponding
    crystallographic point group. It converts them using symmetry operations
    and returns the disorientations

    Parameters
    ----------
    misquats: Quaternion class
        Quaternion misorientations

    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    tol: float
        Tolerance for the disorientation to belong in the fundamental zone

    Returns
    -------
    disquats: quaternion class
        Disorientations for the given misorientations
    """

    file_dir = os.path.dirname(os.path.realpath(__file__))
    fil2 = file_dir+'/pkl_files/symm_quats_'+cryst_ptgrp+'.pkl'
    pkl_fil2 = open(fil2, 'rb')
    symm_quat = pickle.load(pkl_fil2)

    if misquats.ndim == 1:
        misquats = np.reshape(misquats, (5, 1))

    disquats = quat.Quaternion(np.zeros(np.shape(misquats)))
    disquats[:] = np.NaN

    msz = quat.get_size(misquats)
    symm_sz = quat.get_size(symm_quat)
    for ct1 in range(msz):
        tq1 = misquats[:, ct1]
        tcheck = 0
        for j in range(symm_sz):
            tsymm_q1 = quat.mtimes(tq1, symm_quat[:, j])
            for k in range(symm_sz):
                tsymm_q2 = quat.mtimes(symm_quat[:, k], tsymm_q1)

                tqn1 = quat.Quaternion(np.copy(tsymm_q2))
                tqn2 = quat.antipodal(tsymm_q2)
                tqn3 = quat.inverse(tsymm_q2)
                tqn4 = quat.inverse(quat.antipodal(tsymm_q2))

                if check_cond(tqn1, cryst_ptgrp, tol):
                    tcheck = 1
                    disquats[:, ct1] = tqn1
                    break
                elif check_cond(tqn2, cryst_ptgrp, tol):
                    tcheck = 1
                    disquats[:, ct1] = tqn2
                    break
                elif check_cond(tqn3, cryst_ptgrp, tol):
                    tcheck = 1
                    disquats[:, ct1] = tqn3
                    break
                elif check_cond(tqn4, cryst_ptgrp, tol):
                    tcheck = 1
                    disquats[:, ct1] = tqn4
                    break

            if tcheck == 1:
                break
        if tcheck == 0:
            raise Exception('FZ Quaternion not found')
    return disquats
# ------------------------------------------------------------------------------------------------------
