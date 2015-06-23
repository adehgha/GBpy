# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np
import sys

import pickle
import os

import scipy.io

file_dir = os.path.dirname(os.path.realpath(__file__))

# Load Rotations Transformations , Matrix operations,
# Unique Rows With Tolerance, Modules
path_dir = file_dir + '/../../'
sys.path.append(path_dir)
# import transformations as trans
import tools as tl

# Load Quaternion Module
# path_dir = file_dir + '/../../geometry/'
# sys.path.append(path_dir)
import quaternion as quat

# # Load Misorientation Fundamental Zones Module
# path_dir = file_dir + '/../../symmetry/'
# sys.path.append(path_dir)
import misorient_fz as mis_fz


def test_cubic_cslmats(l1):
        """
        Testing the CSL matrices generated using MATLAB with those
        generated using python code

        """
        # Define "ideal" cubic fcc lattice
        l_g_go = l1.l_g_go
        l_go_g = np.linalg.inv(l_g_go)

        cryst_ptgrp = l1.cryst_ptgrp

        # Load matlab generated CSL matrices
        mat_file = '../matlab_csls/432_CommonCSL.mat'
        matf = scipy.io.loadmat(mat_file)
        sig_rots = matf['Sigma_Rots']
        sig_type = 'common'

        # Load Python genererated CSL matrices
        pkl_file = ('../python_csls/' + l1.elem_type +
                    '_csl_' + sig_type + '_rotations' + '.pkl')
        jar1 = open(pkl_file, 'rb')
        sig_rots_p = pickle.load(jar1)

        sig_rots_m = {}
        for ct1 in range(np.size(sig_rots)):
                sig_num = sig_rots[0][ct1][0][0][0][0]
                rot_n = sig_rots[0][ct1][0][1]
                rot_d = sig_rots[0][ct1][0][2]
                sig_rots_m[str(sig_num)] = {}
                if rot_n.ndim == 3:
                        msz = np.shape(rot_n)[2]
                        n_mat = np.zeros((msz, 3, 3))
                        d_mat = np.zeros((msz, 3, 3))
                        for ct2 in range(msz):
                            n_mat[ct2, :, :] = rot_n[:, :, ct2]
                            d_mat[ct2, :, :] = rot_d[:, :, ct2]
                        sig_rots_m[str(sig_num)]['N'] = n_mat
                        sig_rots_m[str(sig_num)]['D'] = d_mat
                elif rot_n.ndim == 2:
                        msz = 1
                        n_mat = np.zeros((msz, 3, 3))
                        d_mat = np.zeros((msz, 3, 3))
                        n_mat[0, :, :] = rot_n
                        d_mat[0, :, :] = rot_d
                        sig_rots_m[str(sig_num)]['N'] = n_mat
                        sig_rots_m[str(sig_num)]['D'] = d_mat
                else:
                        raise Exception('Wrong Dimensions')

                sig_rots_m[str(sig_num)]['N'] = n_mat
                sig_rots_m[str(sig_num)]['D'] = d_mat

        sig_vals = sig_rots_m.keys()
        for ct1 in sig_vals:
                print ct1
                rot_np = sig_rots_p[str(ct1)]['N']
                rot_dp = sig_rots_p[str(ct1)]['D']

                rot_nm = sig_rots_m[str(ct1)]['N']
                rot_dm = sig_rots_m[str(ct1)]['D']

                if rot_nm.ndim == 3:
                        msz1 = np.shape(rot_np)[0]
                        msz2 = np.shape(rot_nm)[0]
                        if msz1 != msz2:
                                raise Exception('No Good')
                        
                        inds = []
                        for ct2 in range(msz1):
                                tn1 = rot_np[ct2, :, :]
                                td1 = rot_dp[ct2, :, :]
                                matp1 = tn1.astype(float) / td1.astype(float)
                                matp2 = np.dot(np.dot(l_g_go, matp1), l_go_g)
                                quat_p2 = tl.mat2quat(matp2)
                                disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)
                                tcheck = 0

                                for ct3 in range(msz2):
                                        tn2 = rot_nm[ct3, :, :]
                                        td2 = rot_dm[ct3, :, :]
                                        mat_m = tn2.astype(float) / td2.astype(float)
                                        quat_m = tl.mat2quat(mat_m)
                                        disquat_m= mis_fz.misorient_fz(quat_m, cryst_ptgrp)
                                        if quat.eq(disquat_p2, disquat_m, 1e-10):
                                                tcheck = 1
                                                inds.append(ct3)
                                                break
                                if (tcheck == 0):
                                        raise Exception('No Good')

                elif rot_nm.ndim == 2:
                        tn1 = rot_np[0, :, :]
                        td1 = rot_dp[0, :, :]
                        matp1 = tn1.astype(float) / td1.astype(float)
                        matp2 = np.dot(np.dot(l_g_go, matp1), l_go_g)
                        quat_p2 = tl.mat2quat(matp2)
                        disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)

                        tn2 = rot_nm[0, :, :]
                        td2 = rot_dm[0, :, :]
                        mat_m = tn2.astype(float) / td2.astype(float)
                        quat_m = tl.mat2quat(mat_m)
                        disquat_m = mis_fz.misorient_fz(quat_m, cryst_ptgrp)

                        # if mat_ops.eq(mat_m, matp2, 1e-10):
                        if quat.eq(disquat_p2, disquat_m, 1e-10):
                                print 'matp1 exists in mat_m'
                        # else:
                        #         raise Exception('No Good')
                else:
                        raise Exception('Wrong Dimensions')
