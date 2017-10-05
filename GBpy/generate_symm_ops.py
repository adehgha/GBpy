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
import shutil
import quaternion as quat
import tools as trans
# -----------------------------------------------------------------------------------------------------------

def generate_symm_mats(cryst_ptgrp, tol=1e-10):
    """
    Give crystallographic point group, this function generates all the symmetry
    operations (as matrices) that belong to the point group using 'generators'

    Parameters
    -----------------
    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    tol: float
        The tolerance used to check if two matrices are the same

    Returns
    ------------
    symm_mat: numpy array
        Size: n x 3 x3 \v
        Symmetry operations as matrices for the corresponding point group
    """

    prop_grps = ['C1', 'C2', 'C3', 'C4', 'C6', 'D2', 'D3', 'D4', 'D6', 'D8'
                 'T', 'O']
    laue_grps = ['Ci', 'C2h', 'C3i', 'C4h', 'C6h', 'D2h', 'D3d', 'D4h', 'D6h', 'D8h',
                 'Th', 'Oh']
    noncentro_grps = ['Cs', 'S4', 'S6', 'C2v', 'C3v', 'C4v', 'C6v', 'D2d', 'D3h', 'Td']

    if cryst_ptgrp in prop_grps:

        if cryst_ptgrp == 'C1':
            n = 1
            gsz = 1
            order_ptgrp = n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C2':
            n = 2
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C3':
            n = 3
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C4':
            n = 4
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C6':
            n = 6
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D2':
            n = 2
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D3':
            n = 3
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D4':
            n = 4
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D6':
            n = 6
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D8':
            n = 8
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'T':
            gsz = 3
            order_ptgrp = 12
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/2, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 1, 1, 2*np.pi/3, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'O':
            gsz = 3
            order_ptgrp = 24
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, np.pi/2, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi/2, 1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

    if cryst_ptgrp in laue_grps:

        if cryst_ptgrp == 'Ci':
            n = 1
            gsz = 2
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C2h':
            n = 2
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C3i':
            n = 3
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C4h':
            n = 4
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'C6h':
            n = 6
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D2h':
            n = 2
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D3d':
            n = 3
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D4h':
            n = 4
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D6h':
            n = 6
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'D8h':
            n = 8
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'Td':
            gsz = 4
            order_ptgrp = 2*12
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 2*np.pi/2, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 1, 1, 2*np.pi/3, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

        elif cryst_ptgrp == 'Oh':
            gsz = 4
            order_ptgrp = 2*24
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, np.pi/2, 1]))
            generators[2, :, :] = trans.vrrotvec2mat(
                np.array([1, 0, 0, np.pi/2, 1]))
            generators[3, :, :] = trans.vrrotvec2mat(
                np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

    if cryst_ptgrp in noncentro_grps:
        if cryst_ptgrp == 'Cs':
            n = 1
            gsz = 2
            order_ptgrp = 2*n
            # Generators
            generators = np.zeros((gsz, 3, 3))
            generators[0, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, 1]))
            generators[1, :, :] = trans.vrrotvec2mat(np.array([0, 0, 1, 0, -1]))
            generators[1, :, :] = np.dot(generators[1, :, :], trans.vrrotvec2mat(
                np.array([0, 0, 1, np.pi, 1])))

            # Symmetry Operations
            symm_mat = np.zeros([order_ptgrp, 3, 3])
            symm_mat[:gsz, :, :] = generators

    count1 = 1
    numops = gsz-1

    while numops < order_ptgrp-1:
        initsize = numops
        for ct1 in np.arange(count1, initsize+1):
            t1 = np.copy(symm_mat[ct1, :, :])
            t2 = np.copy(symm_mat[count1, :, :])
            tM1 = np.dot(t1, t2)
            tcheck = 0
            # for ct2 in np.arange(0, numops+1):
            ct2 = np.arange(0, numops+1)
            if trans.eq(tM1, symm_mat[ct2, :, :], tol):
                tcheck = 1
            if tcheck == 0:
                symm_mat[numops+1, :, :] = np.copy(tM1)
                numops = numops + 1
        if numops == initsize:
            count1 = count1+1


    # print symm_mat
    return symm_mat
# -----------------------------------------------------------------------------------------------------------


def generate_symm_quats(cryst_ptgrp, tol=1e-10):
    """
    Give crystallographic point group, this function generates all the symmetry
    operations (as quaternions) that belong to the point group
    using 'generators'

    Parameters
    -----------------
    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    tol: float
        The tolerance used to check if two matrices are the same

    Returns
    ------------
    symm_quat: quaternion array
        Size: n x 5 \v
        Symmetry operations as matrices for the corresponding point group
    """

    prop_grps = ['C1', 'C2', 'C3', 'C4', 'C6', 'D2', 'D3', 'D4', 'D6',
                 'T', 'O']
    laue_grps = ['Ci', 'C2h', 'C3i', 'C4h', 'C6h', 'D2h', 'D3d', 'D4h', 'D6h', 'D8h',
                 'Th', 'Oh']
    noncentro_grps = ['Cs', 'S4', 'S6', 'C2v', 'C3v', 'C4v', 'C6v', 'D2d', 'D3h', 'Td']

    # if cryst_ptgrp in prop_grps:
    #     rot_type = 'proper'
    # elif cryst_ptgrp in laue_grps:
    #     rot_type = 'improper'
    # else:
    #     raise Exception('Wrong Input for crystal point group')

    if cryst_ptgrp in prop_grps:

        if cryst_ptgrp == 'C1':
            n = 1
            gsz = 1
            order_ptgrp = n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C2':
            n = 2
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C3':
            n = 3
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C4':
            n = 4
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C6':
            n = 6
            gsz = 2
            order_ptgrp = n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D2':
            n = 2
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D3':
            n = 3
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D4':
            n = 4
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D6':
            n = 6
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'T':
            gsz = 3
            order_ptgrp = 12
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/2, 1]))
            generators[:, 2] = trans.axang2quat(
                np.array([1, 1, 1, 2*np.pi/3, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'O':
            gsz = 3
            order_ptgrp = 24
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, np.pi/2, 1]))
            generators[:, 2] = trans.axang2quat(
                np.array([1, 0, 0, np.pi/2, 1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

    if cryst_ptgrp in laue_grps:

        if cryst_ptgrp == 'Ci':
            n = 1
            gsz = 2
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C2h':
            n = 2
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C3i':
            n = 3
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C4h':
            n = 4
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'C6h':
            n = 6
            gsz = 3
            order_ptgrp = 2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D2h':
            n = 2
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D3d':
            n = 3
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D4h':
            n = 4
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D6h':
            n = 6
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'D8h':
            n = 8
            gsz = 4
            order_ptgrp = 2*2*n
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/n, 1]))
            generators[:, 2] = trans.axang2quat(np.array([1, 0, 0, np.pi, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'Td':
            gsz = 4
            order_ptgrp = 2*12
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, 2*np.pi/2, 1]))
            generators[:, 2] = trans.axang2quat(
                np.array([1, 1, 1, 2*np.pi/3, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

        elif cryst_ptgrp == 'Oh':
            gsz = 4
            order_ptgrp = 2*24
            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            generators[:, 1] = trans.axang2quat(
                np.array([0, 0, 1, np.pi/2, 1]))
            generators[:, 2] = trans.axang2quat(
                np.array([1, 0, 0, np.pi/2, 1]))
            generators[:, 3] = trans.axang2quat(np.array([0, 0, 1, 0, -1]))

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators


    if cryst_ptgrp in noncentro_grps:
        if cryst_ptgrp == 'Cs':
            n = 1
            gsz = 2
            order_ptgrp = 2*n
            # Generators

            # Generators
            generators = quat.Quaternion(np.zeros((5, gsz)))
            generators[:, 0] = trans.axang2quat(np.array([0, 0, 1, 0, 1]))
            t_q1 = trans.axang2quat(np.array([0, 0, 1, 0, -1]))
            t_q2 = trans.axang2quat(np.array([0, 0, 1, np.pi, 1]))
            t_mat = quat.mtimes(t_q1, t_q2)
            generators[:, 1] = t_mat

            # Symmetry Operations
            symm_quat = quat.Quaternion(np.zeros([5, order_ptgrp]))
            symm_quat[:, :gsz] = generators

    count1 = 1
    numops = gsz-1

    while numops < order_ptgrp-1:
        initsize = numops
        for ct1 in np.arange(count1, initsize+1):
            t1 = symm_quat[:, ct1]
            t2 = symm_quat[:, count1]
            tM1 = quat.mtimes(t1, t2)
            tcheck = 0
            ct2 = np.arange(0, numops+1)
            if np.any(quat.eq(tM1, symm_quat[:, ct2], tol)):
                tcheck = 1

            if tcheck == 0:
                symm_quat[:, numops+1] = tM1
                numops += 1
        if numops == initsize:
            count1 += 1

    symm_quat = quat.antipodal(symm_quat)

    # print quat.display(symm_quat)
    return symm_quat
# -----------------------------------------------------------------------------------------------------------


def save_symm_pkl(cryst_ptgrp, op_type):
    """
    A pkl file with the symmetry operations of op_type (matrices or
    quaternions) are created and stored in the 'pkl_files' directory

    Parameters
    -----------------
    cryst_ptgrp: string
        Crystallogrphic point group in Schoenflies notation

    op_type: {'matrices', 'quats'}
        Creates matrices or quaternion symmetry operations depending on op_type
    """

    tol = 1e-10
    if op_type == 'matrices':
        symm_ops = generate_symm_mats(cryst_ptgrp, tol)
        fstr = 'mats'
    elif op_type == 'quats':
        symm_ops = generate_symm_quats(cryst_ptgrp, tol)
        fstr = 'quats'

    pkl_file = 'symm_' + fstr + '_' + cryst_ptgrp + '.pkl'
    jar = open(pkl_file, 'wb')
    pickle.dump(symm_ops, jar)
    jar.close()

    # Move to pkl_files
    shutil.move(pkl_file, 'pkl_files/')
# -----------------------------------------------------------------------------------------------------------

# save_symm_pkl('D3', 'quats')
# save_symm_pkl('D4', 'matrices')
# save_symm_pkl('D6', 'quats')
# save_symm_pkl('D6h', 'quats')
# save_symm_pkl('O', 'quats')
# save_symm_pkl('Oh', 'quats')
# save_symm_pkl('C2h', 'quats')
# save_symm_pkl('D8h', 'matrices')

