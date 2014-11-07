# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import sys
import os

file_dir = os.path.dirname(os.path.realpath(__file__))

# Load lattice module
path_dir = file_dir + '/../../lattice/'
sys.path.append(path_dir)
import lattice as lat
# Load  testing python files
import test_cubic_cslmats as t_cubic
import test_lit_cslmats as t_lit

# Test Cases
# 1: Compare python vs lit cubic CSLs
# 2: Compare python vs lit tetrahedral CSLs
# 3: Compare python vs lit hexagonal CSLs
# 4: Compare python vs lit rhombohedral CSLs
test_case = 7

if test_case == 1:
    sig_type = 'common'
    lat_type = 'cF_Id'
    l1 = lat.Lattice()
    t_cubic.test_cubic_cslmats(lat_type)
elif test_case == 2:
    sig_type = 'common'
    lat_type = 'tP_Id'
    l1 = lat.Lattice(lat_type)
    t_lit.test_lit_common_cslmats(lat_type)
elif test_case == 3:
    sig_type = 'common'
    lat_type = 'hP_Id'
    l1 = lat.Lattice(lat_type)
    t_lit.test_lit_common_cslmats(lat_type)
elif test_case == 4:
    sig_type = 'common'
    lat_type = 'hR_Id'
    l1 = lat.Lattice(lat_type)
    t_lit.test_lit_common_cslmats(lat_type)
elif test_case == 5:
    sig_type = 'specific'
    ca_rat = 3
    lat_type = 'tP_ca'
    l1 = lat.Lattice(lat_type, ca_rat)
    t_lit.test_lit_specific_cslmats(lat_type)
elif test_case == 6:
    sig_type = 'specific'
    elem_type = 'hP_ca'
    t_lit.test_lit_specific_cslmats(elem_type)
elif test_case == 7:
    sig_type = 'specific'
    elem_type = 'hR_ca'
    t_lit.test_lit_specific_cslmats(elem_type)
