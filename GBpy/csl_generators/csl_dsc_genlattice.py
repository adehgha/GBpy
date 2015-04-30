# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
#  Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446

import numpy as np
from numpy import dot
from integer_manipulations import rat
from integer_manipulations import lcm_array
from integer_manipulations import int_check
from find_csl_dsc import find_csl_dsc
from sys import exit
import scipy.io as sio
from Col import Col
from csl_finder_smith import check_csl_finder_smith
from dsc_finder import check_dsc_finder
import os
import unittest

""""
# p = mfilename('fullpath');
# Pfd = p(1:length(p)-length(mfilename))
# PathTree=genpath(Pfd) # Returns the pwd (present working directory) and any folder within that file name
# addpath(PathTree) # Adds PathTree to the top of the search path so any folder within it can be easily accessed by the following code
"""


def csl_dsc_genlattice():
    check = True
    # ### Vectors for the basis of a primitive FCC cell
    b1x = np.array([[0.], [1.], [1.]]) / 2
    b1y = np.array([[1.], [0.], [1.]]) / 2
    b1z = np.array([[1.], [1.], [0.]]) / 2
    L_G1_GO1 = np.concatenate((b1x, b1y, b1z), axis=1)
    L_GO1_G1 = np.linalg.inv(L_G1_GO1)

    ### Cubic Systems
    ### Common
    mat_contents = sio.loadmat('432_CommonCSL.mat')
    for cnt in range(mat_contents['Sigma_Rots'].shape[1]):
        print 'case: ', cnt + 1
        Sigma_CSL = mat_contents['Sigma_Rots'][0, cnt][0][0]
        R_N_tmp = mat_contents['Sigma_Rots'][0, cnt][0][1]
        R_D_tmp = mat_contents['Sigma_Rots'][0, cnt][0][2]
        R_N = R_N_tmp.astype(float, copy=True)
        R_D = R_D_tmp.astype(float, copy=True)
        # % The rotation matrix is obtained (numerator matrix)/denominator matrix
        if len(R_N.shape) == 2:
            R_G1toG2_GO1 = R_N / R_D
        else:
            ct2 = 0
            R_G1toG2_GO1 = R_N[:, :, ct2] / R_D[:, :, ct2]

        R_G1toG2_G1 = dot(dot(L_GO1_G1, R_G1toG2_GO1), L_G1_GO1)

        ##
        # Find_CSL_DSC computes the CSL and DSC lattice using the lattice unit
        # cell and rotation matrix as parameters
        [L_CSL_G1, L_DSC_G1] = find_csl_dsc(L_G1_GO1, R_G1toG2_G1)
        ###########printing and checking #############
        txt = Col()
        txt.c_prnt('CSL', 'yel')
        print L_CSL_G1
        if check == True :check_csl_finder_smith(R_G1toG2_G1, Sigma_CSL, L_G1_GO1, L_CSL_G1)
        txt.c_prnt('DSC', 'dgrn')
        print L_DSC_G1
        if check == True: check_dsc_finder(R_G1toG2_G1, Sigma_CSL, L_G1_GO1, L_DSC_G1, L_CSL_G1)
        print '\n-------------------'
        ##############################################

    return L_CSL_G1, L_DSC_G1

# -----------------------------------------------------------------------------------------------------------



