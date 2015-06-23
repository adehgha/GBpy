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
path_dir = file_dir + '/../../'
sys.path.append(path_dir)
import tools as tl
# import transformations as trans

# # Load lattice module
# path_dir = file_dir + '/../../lattice/'
# sys.path.append(path_dir)
# import lattice as lat

# # Load Quaternion Module
# path_dir = file_dir + '/../../geometry/'
# sys.path.append(path_dir)
import quaternion as quat

# # Load Misorientation Fundamental Zones Module
# path_dir = file_dir + '/../../symmetry/'
# sys.path.append(path_dir)
import misorient_fz as mis_fz


# path_dir = file_dir + '/../../'
# sys.path.append(path_dir)
import integer_manipulations as int_man

path_dir = file_dir + '/../'
sys.path.append(path_dir)
import csl_utility_functions as csl_util


####
def compute_sig_rots(csl_rots, l1, sig_type):
        [tau, kmax] = csl_util.compute_inp_params(l1, sig_type)
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
                sig_rots_m[ct1][:] = np.NAN

        for ct1 in k1:
                if ct1[-1] in string.ascii_lowercase:
                        tstr1 = ct1[0:-1]
                else:
                        tstr1 = ct1

                if tstr1 in sig_inds1.keys():
                        sig_inds1[tstr1] += 1
                        if csl_rots['type'] == 'quads':
                                if l1.pearson[0:2]=='hR':
                                        m = csl_rots_m[ct1][0]
                                        u = csl_rots_m[ct1][1]
                                        v = csl_rots_m[ct1][2]
                                        w = csl_rots_m[ct1][3]
                                        U = (2*u + v + w)/3.0
                                        V = (-u + v + w)/3.0
                                        W = (w - u - 2*v)/3.0
                                else:
                                        m = csl_rots_m[ct1][0]
                                        U = csl_rots_m[ct1][1]
                                        V = csl_rots_m[ct1][2]
                                        W = csl_rots_m[ct1][3]
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = \
                                                                    csl_util.compute_tmat(np.array([m,U,V,W]).reshape(4, 1), tau, l1)
                        elif csl_rots['type'] == 'matrices':
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = csl_rots_m[ct1]
                else:
                        sig_inds1[tstr1] = 0
                        if csl_rots['type'] == 'quads':
                                if l1.pearson[0:2]=='hR':
                                        m = csl_rots_m[ct1][0]
                                        u = csl_rots_m[ct1][1]
                                        v = csl_rots_m[ct1][2]
                                        w = csl_rots_m[ct1][3]
                                        U = (2*u + v + w)/3.0
                                        V = (-u + v + w)/3.0
                                        W = (w - u - 2*v)/3.0
                                else:
                                        m = csl_rots_m[ct1][0]
                                        U = csl_rots_m[ct1][1]
                                        V = csl_rots_m[ct1][2]
                                        W = csl_rots_m[ct1][3]
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = \
                                                                    csl_util.compute_tmat(np.array([m,U,V,W]).reshape(4, 1), tau, l1)
                        elif csl_rots['type'] == 'matrices':
                                sig_rots_m[tstr1][sig_inds1[tstr1], :, :] = csl_rots_m[ct1]

        sig_rots_m1 = {}
        for ct1 in sig_rots_m.keys():
                sig_rots_m1[ct1] ={}
                if csl_rots['type'] == 'quads':
                        mats1 = sig_rots_m[ct1]
                        msz1 = np.shape(mats1)[0]
                        sig_rots_m1[ct1]['N'] = np.zeros((msz1, 3, 3))
                        sig_rots_m1[ct1]['D'] = np.zeros((msz1, 3, 3))
                        for ct2 in range(msz1):
                                [mult, num_mat] = int_man.int_mult(mats1[ct2, :, :], 1e-06)
                                if int_man.int_check(mult):
                                        mult = np.round(mult)
                                        if int(mult) != int(ct1):
                                                raise Exception('Not a sigma rotation')
                                else:
                                        raise Exception('Not a sigma rotation')

                                sig_rots_m1[ct1]['N'][ct2, :, :] = num_mat
                                sig_rots_m1[ct1]['D'][ct2, :, :] = mult*np.ones(np.shape(num_mat))
                elif csl_rots['type'] == 'matrices':
                        sig_rots_m1[ct1]['N'] = sig_rots_m[ct1]
                        sig_rots_m1[ct1]['D'] = \
                                        float(ct1)*np.ones(np.shape(sig_rots_m[ct1]))

        return sig_rots_m1

def check_sigma_thm(r_g1tog2_g1, sigma):
        """
        """
        if r_g1tog2_g1.ndim == 2:
                r_g1tog2_g1 = np.reshape(r_g1tog2_g1, (1, 3, 3))

        msz = np.shape(r_g1tog2_g1)[0]
        check_mats = np.zeros((msz, 3, 3))

        ct1 = 0
        for i in range(msz):        
                tmat1 = r_g1tog2_g1[i, :, :]
                tmat2 = np.linalg.inv(r_g1tog2_g1[i, :, :])
                if (int_man.int_check(sigma*tmat1).all() and 
                    int_man.int_check(sigma*np.linalg.det(tmat1)*tmat2).all()):
                        check_mats[ct1, :, :] = r_g1tog2_g1[i, :, :]
                        ct1 += 1
                else:
                        print('Not a Sigma Rotation')

        tarr1 = np.arange(ct1, msz)
        return np.delete(check_mats, tarr1, 0)


def disorients(sig_rots_m1, l1):
        """
        """
        l_g_go = l1.l_g_go
        cryst_ptgrp = csl_util.proper_ptgrp(l1.cryst_ptgrp)

        sig_rots_m2 = {}
        for ct1 in sig_rots_m1.keys():
                n_mats = sig_rots_m1[ct1]['N']
                d_mats = sig_rots_m1[ct1]['D']
                msz = np.shape(n_mats)[0]
                r_g1tog2_g1 = np.zeros((msz, 3, 3))
                for ct2 in range(msz):
                        r_g1tog2_g1[ct2, :, :] = n_mats[ct2, :, :]/d_mats[ct2, :, :]
                r_g1tog2_g1 = check_sigma_thm(r_g1tog2_g1, int(ct1))
                r_g1tog2_g1 = csl_util.disorient_sigmarots(r_g1tog2_g1, l_g_go, cryst_ptgrp)
                sig_rots_m2[ct1] = csl_util.check_sigma_rots(r_g1tog2_g1, int(ct1))

        return sig_rots_m2


def compare_sig_rots(sig_rots_m1, sig_rots_p, l1):
        """
        """
        l_g_go = l1.l_g_go
        l_go_g = np.linalg.inv(l_g_go)

        cryst_ptgrp = csl_util.proper_ptgrp(l1.cryst_ptgrp)
        l1.cryst_ptgrp = cryst_ptgrp

        for ct1 in sig_rots_m1.keys():
                rot_np = sig_rots_p[ct1]['N']
                rot_dp = sig_rots_p[ct1]['D']

                rot_nm = sig_rots_m1[ct1]['N']
                rot_dm = sig_rots_m1[ct1]['D']

                if rot_nm.ndim == 3:
                        msz1 = np.shape(rot_np)[0]
                        msz2 = np.shape(rot_nm)[0]

                        if msz1 < msz2:
                                raise Exception('No Good')
                        if msz1 > msz2:
                                print('msz1= %d \t msz2 = %d \n' %(msz1, msz2))

                        inds = []
                        for ct3 in range(msz2):
                                tn2 = rot_nm[ct3, :, :]
                                td2 = rot_dm[ct3, :, :]
                                mat_m1 = tn2.astype(float) / td2.astype(float)
                                mat_m2 = np.dot(np.dot(l_g_go, mat_m1), l_go_g)
                                quat_m = tl.mat2quat(mat_m2)
                                disquat_m= mis_fz.misorient_fz(quat_m, cryst_ptgrp)

                                tcheck = 0

                                for ct2 in range(msz1):
                                        tn1 = rot_np[ct2, :, :]
                                        td1 = rot_dp[ct2, :, :]
                                        matp1 = tn1.astype(float) / td1.astype(float)
                                        matp2 = np.dot(np.dot(l_g_go, matp1), l_go_g)
                                        quat_p2 = tl.mat2quat(matp2)
                                        disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)

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
                        quat_p2 = tl.mat2quat(matp2)
                        disquat_p2 = mis_fz.misorient_fz(quat_p2, cryst_ptgrp)

                        tn2 = rot_nm[0, :, :]
                        td2 = rot_dm[0, :, :]
                        mat_m1 = tn2.astype(float) / td2.astype(float)
                        mat_m2 = np.dot(np.dot(l_g_go, mat_m1), l_go_g)
                        quat_m = tl.mat2quat(mat_m2)                        
                        disquat_m = mis_fz.misorient_fz(quat_m, cryst_ptgrp)

                        # if mat_ops.eq(mat_m, matp2, 1e-10):
                        if quat.eq(disquat_p2, disquat_m, 1e-10):
                                print 'matp1 exists in mat_m'
                        else:
                                raise Exception('No Good')
                else:
                        raise Exception('Wrong Dimensions')

####
def test_lit_common_cslmats(elem_type):
        """
        """
        sig_type = 'common'

        # Load prior lit cls
        mat_file = '../prior_lit_csls/'+ elem_type + '_lit_csl_' + sig_type + '_rotations' + '.pkl'
        jar1 = open(mat_file, 'rb')
        sig_rots = pickle.load(jar1)
        
        csl_rots = sig_rots['common']
        l1 = csl_rots['lattice']
        sig_rots_m1 = compute_sig_rots(csl_rots, l1, sig_type)
        sig_rots_m2 = disorients(sig_rots_m1, l1)

        for ct1 in sig_rots_m2.keys():
                sig_rots_m3={}
                sig_rots_p={}
                sig_rots_p[ct1] = csl_util.csl_rotations(int(ct1), sig_type, l1)
                sig_rots_m3[ct1] = sig_rots_m2[ct1]
                compare_sig_rots(sig_rots_m3, sig_rots_p, l1)

####
def test_lit_specific_cslmats(elem_type, sig_type):
        """
        """
        # sig_type = 'specific'
        # sig_type = 'common'

        # Load prior lit cls
        mat_file = '../prior_lit_csls/'+ elem_type + '_lit_csl_' + sig_type + '_rotations' + '.pkl'
        jar1 = open(mat_file, 'rb')
        sig_rots = pickle.load(jar1)

        for ct1 in sig_rots.keys():
                print ct1
                csl_rots = sig_rots[ct1]
                l1 = csl_rots['lattice']
                sig_rots_m1 = compute_sig_rots(csl_rots, l1, sig_type)
                sig_rots_m2 = disorients(sig_rots_m1, l1)
                for ct2 in sig_rots_m2.keys():
                        sig_rots_m3={}
                        sig_rots_p={}
                        sig_rots_p[ct2] = csl_util.csl_rotations(int(ct2), sig_type, l1)
                        sig_rots_m3[ct2] = sig_rots_m2[ct2]
                        compare_sig_rots(sig_rots_m3, sig_rots_p, l1)

# test_lit_specific_cslmats('hP_Id', 'common')
# test_lit_specific_cslmats('hR_Id', 'common')
# test_lit_specific_cslmats('tP_Id', 'common')
# test_lit_specific_cslmats('hP_ca', 'specific')
# test_lit_specific_cslmats('hR_ca', 'specific')
test_lit_specific_cslmats('tP_ca', 'specific')

