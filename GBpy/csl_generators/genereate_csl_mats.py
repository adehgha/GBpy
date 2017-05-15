# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
#  Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446

import numpy as np
import sys
import pickle
import os
import GBpy

file_dir = os.path.dirname(os.path.realpath(__file__))
path_dir3 = file_dir + '/../'
sys.path.append(path_dir3)
# Load Integer Manipulations Module
import GBpy.integer_manipulations as int_man
# Load Lattice Module
import GBpy.lattice as lat
# Load CSL Utility function
import GBpy.csl_utility_functions as csl_util

# Test Cases
# 1: Common rotations for cubic lattices
# 2: Common rotations for primitive tetrahedral lattices
# 3: Common rotations for primitive hexagonal lattices
# 4: Common rotations for primitive rhombohedral lattices

test_case = 6
# Input parameters for pkl files
if test_case == 1:
    sig_type = 'common'
    l1 = lat.Lattice()
elif test_case == 2:
    sig_type = 'common'
    lat_type = 'tP_Id'
    l1 = lat.Lattice(lat_type)
elif test_case == 3:
    sig_type = 'common'
    lat_type = 'hP_Id'
    l1 = lat.Lattice(lat_type)
elif test_case == 4:
    sig_type = 'common'
    lat_type = 'hR_Id'
    l1 = lat.Lattice(lat_type)
elif test_case == 5:
    sig_type = 'specific'
    ca_rat = 3
    lat_type = 'tP_ca'
    l1 = lat.Lattice(lat_type, ca_rat)
elif test_case == 6:
    sig_type = 'common'
    lat_type = 'cP_Id'
    l1 = lat.Lattice(lat_type)

sig_rots = {}
sig_num = np.arange(100) + 1

for i in sig_num:
    print i
    sig_rots[str(i)] = csl_util.csl_rotations(i, sig_type, l1)

if sig_type == 'common':
    pkl_file = l1.elem_type + '_csl_' + sig_type + '_rotations' + '.pkl'
    jar = open(pkl_file, 'wb')
    pickle.dump(sig_rots, jar)
    jar.close()

if sig_type == 'specific':
    lat_tau = l1.lat_params['a']**2/l1.lat_params['c']**2
    [N, D] = int_man.rat(lat_tau)
    pkl_file = l1.elem_type + '_csl_' + sig_type + \
               '_tau_'+str(N[0][0])+'_'+str(D[0][0])+'_rotations' + '.pkl'
    jar = open(pkl_file, 'wb')
    pickle.dump(sig_rots, jar)
    jar.close()
