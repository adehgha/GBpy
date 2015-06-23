# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np
import sys

import pickle
import os

import string

file_dir = os.path.dirname(os.path.realpath(__file__))
# Load Rotations Transformations , Matrix operations,
# Unique Rows With Tolerance, Modules
path_dir = file_dir + '/../../tools/'
sys.path.append(path_dir)
import transformations as trans

# # Load lattice module
# path_dir = file_dir + '/../../lattice/'
# sys.path.append(path_dir)
# import lattice as lat

# Load Quaternion Module
path_dir = file_dir + '/../../geometry/'
sys.path.append(path_dir)
import quaternion as quat

# Load Misorientation Fundamental Zones Module
path_dir = file_dir + '/../../symmetry/'
sys.path.append(path_dir)
import misorient_fz as mis_fz

def test_tet_common_cslmats(l1):
        """
        """
        # Define "ideal" cubic fcc lattice
        l_g_go = l1.l_g_go
        l_go_g = np.linalg.inv(l_g_go)

        cryst_ptgrp = l1.cryst_ptgrp
        sig_type = 'common'

        # Load prior lit cls
        mat_file = '../prior_lit_csls/'+ l1.elem_type + \
                   '_lit_csl_' + sig_type + '_rotations' + '.pkl'
        jar1 = open(mat_file, 'rb')
        csl_rots = pickle.load(jar1)
        
        if csl_rots['type'] == 'matrices':
                csl_rots_m = csl_rots['rots']
                k1 = csl_rots_m.keys()
                sig_inds = {}
                for ct1 in k1:                        
                        if ct1[-1] in string.ascii_lowercase:
                                tstr1 = ct1[0:-1]
                        else:
                                tstr1 = ct1
                        if tstr1 in sig_inds.keys():
                                sig_inds[tstr1] += 1
                        else:
                                sig_inds[tstr1] = 1

                sig_rots_m = {}
                sig_inds1 = {}
                for ct1 in sig_inds.keys():
                        sig_rots_m[ct1] = np.zeros((sig_inds[ct1], 3, 3))
                for ct1 in k1:
                        if ct1[-1] in string.ascii_lowercase:
                                tstr1 = ct1[0:-1]
                        else:
                                tstr1 = ct1

                        if tstr1 in sig_inds1.keys():
                                sig_inds1[tstr1] += 1
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = csl_rots_m[ct1]
                        else:
                                sig_inds1[tstr1] = 0
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = csl_rots_m[ct1]

                sig_rots_m1 = {}
                for ct1 in sig_rots_m.keys():
                        sig_rots_m1[ct1] ={}
                        sig_rots_m1[ct1]['N'] = sig_rots_m[ct1]
                        sig_rots_m1[ct1]['D'] = \
                                                float(ct1)*np.ones(np.shape(sig_rots_m[ct1]))
                        
                        
        # Load Python genererated CSL matrices
        pkl_file = '../python_csls/'+l1.elem_type + \
                   '_csl_' + sig_type + '_rotations' + '.pkl'
        jar2 = open(pkl_file, 'rb')
        sig_rots_p = pickle.load(jar2)

        for ct1 in sig_rots_m1.keys():
                rot_np = sig_rots_p[ct1]['N']
                rot_dp = sig_rots_p[ct1]['D']

                rot_nm = sig_rots_m1[ct1]['N']
                rot_dm = sig_rots_m1[ct1]['D']

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
                                quat_p2 = trans.mat2quat(matp2)
                                disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)
                                tcheck = 0

                                for ct3 in range(msz2):
                                        tn2 = rot_nm[ct3, :, :]
                                        td2 = rot_dm[ct3, :, :]
                                        mat_m = tn2.astype(float) / td2.astype(float)
                                        quat_m = trans.mat2quat(mat_m)
                                        disquat_m= mis_fz.misorient_fz(quat_m, cryst_ptgrp)
                                        if quat.eq(disquat_p2, disquat_m, 1e-10):
                                                tcheck = 1
                                                inds.append(ct3)
                                                break
                                if (tcheck == 0):
                                        raise Exception('No Good')
                                else:
                                        print 'matp1 exists in mat_m'

                elif rot_nm.ndim == 2:
                        tn1 = rot_np[0, :, :]
                        td1 = rot_dp[0, :, :]
                        matp1 = tn1.astype(float) / td1.astype(float)
                        matp2 = np.dot(np.dot(l_g_go, matp1), l_go_g)
                        quat_p2 = trans.mat2quat(matp2)
                        disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)

                        tn2 = rot_nm[0, :, :]
                        td2 = rot_dm[0, :, :]
                        mat_m = tn2.astype(float) / td2.astype(float)
                        quat_m = trans.mat2quat(mat_m)
                        disquat_m = mis_fz.misorient_fz(quat_m, cryst_ptgrp)

                        # if mat_ops.eq(mat_m, matp2, 1e-10):
                        if quat.eq(disquat_p2, disquat_m, 1e-10):
                                print 'matp1 exists in mat_m'
                        else:
                                raise Exception('No Good')
                else:
                        raise Exception('Wrong Dimensions')



# test_case = 1;
# if test_case == 1:
#     sig_type = 'common'
#     l1 = lat.Lattice()
#     test_cubic_cslmats(l1)
# elif test_case == 2:
#     sig_type = 'common'
#     lat_type = 'tP_Id'
#     l1 = lat.Lattice(lat_type)
# elif test_case == 3:
#     sig_type = 'common'
#     lat_type = 'tP_Id'
#     l1 = lat.Lattice(lat_type)
# elif test_case == 4:
#     sig_type = 'common'
#     lat_type = 'tP_Id'
#     l1 = lat.Lattice(lat_type)
# elif test_case == 5:
#     sig_type = 'specific'
#     ca_rat = 3
#     lat_type = 'tP_ca'
#     l1 = lat.Lattice(lat_type, ca_rat)
