# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
#  Journal  of Applied Crystallography 48, 585 - 588 (2015). doi:10.1107/S1600576715004446
import numpy as np
import pickle

import os
file_dir = os.path.dirname(os.path.realpath(__file__))

import sys
# Load Lattice Module
path_lat = file_dir + '/../../lattice/'
sys.path.append(path_lat)
import lattice as lat

path_dir = file_dir + '/../../'
sys.path.append(path_dir)
import integer_manipulations as int_man

## Test Cases
# 1: Common rotations for cubic lattices
# 2: Common rotations for primitive tetrahedral lattices
# 3: Common rotations for primitive hexagonal lattices
# 4: Common rotations for primitive rhombohedral lattices

test_case = 7

def save_csl_rots(sig_rots, sig_type, lat1):
    if sig_type == 'common':
        pkl_file = lat1.elem_type + '_lit_csl_' + sig_type + '_rotations' + '.pkl'
        jar = open(pkl_file, 'wb')
        pickle.dump(sig_rots, jar)
        jar.close()
    if sig_type == 'specific':
        pkl_file = lat1.elem_type + '_lit_csl_' + sig_type + '_rotations' + '.pkl'
        jar = open(pkl_file, 'wb')
        pickle.dump(sig_rots, jar)
        jar.close()

def lit_csl_rots (lat1, sig_type):
    """
    """

    csl_rots = {}
    csl_rots['lattice'] = lat1
    csl_rots['rots'] = {}
    csl_rots['sig_type'] = sig_type

    lat_elem = lat1.elem_type
    if lat_elem == 'tP_Id':
        if sig_type == 'common':
            csl_rots['type'] = 'matrices'
            csl_rots['rots']['5'] = np.array([[4, -3, 0], [3, 4, 0], [0, 0, 5]])
            csl_rots['rots']['13'] = np.array([[12, -5, 0], [5, 12, 0], [0, 0, 13]])
            csl_rots['rots']['17'] = np.array([[15, -8, 0], [8, 15, 0], [0, 0, 17]])
            csl_rots['rots']['25'] = np.array([[24, -7, 0], [7, 24, 0], [0, 0, 25]])
            csl_rots['rots']['29'] = np.array([[21, -20, 0], [20, 21, 0], [0, 0, 29]])
            csl_rots['rots']['37'] = np.array([[35, -12, 0], [12, 35, 0], [0, 0, 37]])
            csl_rots['rots']['41'] = np.array([[40, -9, 0], [9, 40, 0], [0, 0, 41]])

    if lat_elem == 'hP_Id':
        if sig_type == 'common':
            csl_rots['type'] = 'quads'
            csl_rots['rots']['7'] = np.array([3, 0, 0, 1])
            csl_rots['rots']['13'] = np.array([7, 0, 0, 3])
            csl_rots['rots']['19'] = np.array([5, 0, 0, 1])
            csl_rots['rots']['31'] = np.array([11, 0, 0, 3])
            csl_rots['rots']['37'] = np.array([7, 0, 0, 1])
            csl_rots['rots']['43'] = np.array([13, 0, 0, 3])
            csl_rots['rots']['49'] = np.array([4, 0, 0, 1])

    if lat_elem == 'hR_Id':
        if sig_type == 'common':
            csl_rots['type'] = 'quads'
            csl_rots['rots']['3'] = np.array([3, 0, 0, 3])
            csl_rots['rots']['7'] = np.array([5, 0, 0, 3])
            csl_rots['rots']['13'] = np.array([7, 0, 0, 3])
            csl_rots['rots']['19'] = np.array([4, 0, 0, 3])
            csl_rots['rots']['21'] = np.array([9, 0, 0, 3])
            csl_rots['rots']['31'] = np.array([11, 0, 0, 3])
            csl_rots['rots']['37'] = np.array([11, 0, 0, 9])
            csl_rots['rots']['39'] = np.array([6, 0, 0, 3])
            csl_rots['rots']['43'] = np.array([13, 0, 0, 3])
            csl_rots['rots']['49'] = np.array([13, 0, 0, 9])
            csl_rots['rots']['57'] = np.array([15, 0, 0, 3])
    
    if lat_elem == 'tP_ca':
        if sig_type == 'specific':
            lat_tau = lat1.lat_params['a']**2/lat1.lat_params['c']**2
            csl_rots['type'] = 'matrices'

            if lat_tau == 1./9.0:
                csl_rots['rots']['3'] = np.array([[3, 0, 0], [0, 0, -9], [0, 1, 0]])
                csl_rots['rots']['5a'] = np.array([[5, 0, 0], [0, 4, -9], [0, 1, 4]])
                csl_rots['rots']['5b'] = np.array([[0, -4, 9], [5, 0, 0], [0, 1, 4]])
                csl_rots['rots']['7'] = np.array([[6, 2, 9], [2, 3, -18], [-1, 2, 2]])
                csl_rots['rots']['9a'] = np.array([[6, 3, 18], [3, 6, -18], [-2, 2, 3]])
                csl_rots['rots']['9b'] = np.array([[6, -6, 9], [6, 3, -18], [1, 2, 6]])
                csl_rots['rots']['9c'] = np.array([[6, -3, 18], [6, 6, -9], [-1, 2, 6]])
                csl_rots['rots']['11a'] = np.array([[9, 2, 18], [2, 9, -18], [-2, -2, 7]])
                csl_rots['rots']['11b'] = np.array([[6, 2, 27], [7, 6, -18], [-2, 3, 2]])
                csl_rots['rots']['13a'] = np.array([[13, 0, 0], [0, 5, -36], [0, 4, 5]])
                csl_rots['rots']['13b'] = np.array([[4, -3, 36], [12, 4, -9], [-1, 4, 4]])
                csl_rots['rots']['13c'] = np.array([[12, -4, 9], [4, 3, -36], [1, 4, 4]])
                csl_rots['rots']['15a'] = np.array([[15, 0, 0], [0, 9, -36], [0, 4, 9]])
                csl_rots['rots']['15b'] = np.array([[12, 0, 27], [9, 0, -36], [0, 5, 0]])
                csl_rots['rots']['17a'] = np.array([[17, 0, 0], [0, 8, -45], [0, 5, 8]])
                csl_rots['rots']['17b'] = np.array([[9, 8, 36], [8, 9, -36], [-4, 4, 1]])
                csl_rots['rots']['17c'] = np.array([[8, -12, 27], [12, -1, -36], [3, 4, 8]])
                csl_rots['rots']['17d'] = np.array([[12, 1, 36], [8, 12, -27], [-3, 4, 8]])
                csl_rots['rots']['19a'] = np.array([[18, 1, 18], [1, 18, -18], [-2, 2, 17]])
                csl_rots['rots']['19b'] = np.array([[17, 6, 18], [6, 1, -54], [-2, 6, -1]])
                csl_rots['rots']['19c'] = np.array([[10, -6, 45], [15, 10, -18], [-2, 5, 10]])
                csl_rots['rots']['19d'] = np.array([[17, -6, 18], [6, -1, -54], [2, 6, 1]])
                csl_rots['rots']['19e'] = np.array([[15, -10, 18], [10, 6, -45], [2, 5, 10]])
                csl_rots['rots']['21a'] = np.array([[18, -6, 27], [9, 18, -18], [-2, 3, 18]])
                csl_rots['rots']['21b'] = np.array([[9, -6, 54], [18, 9, -18], [-2, 6, 9]])
                csl_rots['rots']['21c'] = np.array([[18, -9, 18], [9, 6, -54], [2, 6, 9]])
                csl_rots['rots']['23a'] = np.array([[14, -18, 9], [18, 13, -18], [1, 2, 22]])
                csl_rots['rots']['23b'] = np.array([[18, -13, 18], [14, 18, -9], [-1, 2, 22]])
                csl_rots['rots']['23c'] = np.array([[22, 3, 18], [3, 14, -54], [-2, 6, 13]])
                csl_rots['rots']['23d'] = np.array([[22, -6, 9], [6, 13, -54], [1, 6, 14]])
                csl_rots['rots']['23e'] = np.array([[6, -13, 54], [22, 6, -9], [-1, 6, 14]])
                csl_rots['rots']['25a'] = np.array([[25, 0, 0], [0, 7, -72], [0, 8, 7]])
                csl_rots['rots']['25b'] = np.array([[15, -16, 36], [20, 12, -27], [0, 5, 20]])
                csl_rots['rots']['25c'] = np.array([[20, -12, 27], [15, 16, -36], [0, 5, 20]])
                csl_rots['rots']['25d'] = np.array([[9, -20, 36], [20, 0, -45], [4, 5, 16]])
                csl_rots['rots']['25e'] = np.array([[20, 0, 45], [9, 20, -36], [-4, 5, 16]])
                csl_rots['rots']['27a'] = np.array([[24, 3, 36], [3, 24, -36], [-4, 4, 21]])
                csl_rots['rots']['27b'] = np.array([[12, -21, 36], [24, 12, -9], [-1, 4, 24]])
                csl_rots['rots']['27c'] = np.array([[21, 12, 36], [12, 3, -72], [-4, 8, -3]])
                csl_rots['rots']['27d'] = np.array([[24, -12, 9], [12, 21, -36], [1, 4, 24]])
                csl_rots['rots']['27e'] = np.array([[3, -24, 36], [24, -3, -36], [4, 4, 21]])
                csl_rots['rots']['27f'] = np.array([[21, -12, 36], [12, -3, -72], [4, 8, 3]])
                csl_rots['rots']['27g'] = np.array([[12, -12, 63], [24, 3, -36], [1, 8, 12]])
                csl_rots['rots']['27h'] = np.array([[12, 3, 72], [21, 12, -36], [-4, 8, 3]])
                csl_rots['rots']['27i'] = np.array([[24, -3, 36], [12, 12, -63], [-1, 8, 12]])
                csl_rots['rots']['29a'] = np.array([[29, 0, 0], [0, 20, -63], [0, 7, 20]])
                csl_rots['rots']['29b'] = np.array([[16, -12, 63], [24, 11, -36], [-1, 8, 16]])
                csl_rots['rots']['29c'] = np.array([[24, -11, 36], [16, 12, -63], [1, 8, 16]])
                csl_rots['rots']['29d'] = np.array([[21, -16, 36], [16, 3, -72], [4, 8, 11]])
                csl_rots['rots']['29e'] = np.array([[16, -3, 72], [21, 16, -36], [-4, 8, 11]])
                csl_rots['rots']['31a'] = np.array([[21, -14, 54], [22, 6, -63], [2, 9, 14]])
                csl_rots['rots']['31b'] = np.array([[30, 5, 18], [5, 6,-90 ], [-2, 10, 5]])
                csl_rots['rots']['31c'] = np.array([[22, -6, 63], [21, 14, -54], [-2, 9, 14]])
                csl_rots['rots']['31d'] = np.array([[27, -14, 18], [14, 18, -63], [2, 7, 22]])
                csl_rots['rots']['31e'] = np.array([[14, -18, 63], [27, 14, -18], [-2, 7, 22]])
                csl_rots['rots']['33a'] = np.array([[18, -18, 63], [27, 6, -54], [2, 9, 18]])
                csl_rots['rots']['33b'] = np.array([[18, -21, 54], [27, 18, -18], [-2, 6, 27]])
                csl_rots['rots']['33c'] = np.array([[27, -6, 54], [18, 18, -63], [-2, 9, 18]])
                csl_rots['rots']['35a'] = np.array([[17, -6, 90], [30, 10, -45], [-2, 11, 10]])
                csl_rots['rots']['35b'] = np.array([[26, -18, 45], [18, 1, -90], [5, 10, 10]])
                csl_rots['rots']['35c'] = np.array([[26, 15, 54], [15, 10, -90], [-6, 10, 1]])
                csl_rots['rots']['35d'] = np.array([[18, -1, 90], [26, 18, -45], [-5, 10, 10]])
                csl_rots['rots']['35e'] = np.array([[30, -10, 45], [17, 6, -90], [2, 11, 10]])
                csl_rots['rots']['35f'] = np.array([[33, -10, 18], [10, 15, -90], [2, 10, 17]])
                csl_rots['rots']['35g'] = np.array([[10, -30, 45], [30, 1, -54], [5, 6, 26]])
                csl_rots['rots']['35h'] = np.array([[30, -1, 54], [10, 30, -45], [-5, 6, 26]])
                csl_rots['rots']['35i'] = np.array([[15, 10, 90], [26, 15, -54], [-6, 10, -1]])
                csl_rots['rots']['35j'] = np.array([[10, -15, 90], [33, 10, -18], [-2, 10, 17]])
                csl_rots['rots']['37a'] = np.array([[37, 0, 0], [0, 35, -36], [0, 4, 35]])
                csl_rots['rots']['37b'] = np.array([[0, -35, 36], [37, 0, 0], [0, 4, 35]])
                csl_rots['rots']['37c'] = np.array([[28, -12, 63], [21, 28, -36], [-4, 7, 28]])
                csl_rots['rots']['37d'] = np.array([[24, -8, 81], [28, 3, -72], [1, 12, 8]])
                csl_rots['rots']['37e'] = np.array([[28, -3, 72], [24, 8, -81], [-1, 12, 8]])
                csl_rots['rots']['37f'] = np.array([[36, -8, 9], [8, 27, -72], [1, 8, 28]])
                csl_rots['rots']['37g'] = np.array([[21, -28, 36], [28, 12, -63], [4, 7, 28]])
                csl_rots['rots']['39a'] = np.array([[39, 0, 0], [0, 36, -45], [0, 5, 36]])
                csl_rots['rots']['39b'] = np.array([[12,-36 ,27 ], [36, 9, -36], [3, 4, 36]])
                csl_rots['rots']['39c'] = np.array([[36, 9, 36], [9, 12, -108], [-4, 12, 9]])
                csl_rots['rots']['39d'] = np.array([[36, -9, 36], [12, 36, -27], [-3, 4, 36]])
                csl_rots['rots']['39e'] = np.array([[0, -36, 45], [39, 0, 0], [0, 5, 36]])
                csl_rots['rots']['41a'] = np.array([[41, 0, 0], [0, 40, -27], [0, 3, 40]])
                csl_rots['rots']['41b'] = np.array([[32, 9, 72], [9, 32, -72], [-8, 8, 23]])
                csl_rots['rots']['41c'] = np.array([[39, 4, 36], [4, 33, -72], [-4, 8, 31]])
                csl_rots['rots']['41d'] = np.array([[32, -24, 27], [24, 23, -72], [3, 8, 32]])
                csl_rots['rots']['41e'] = np.array([[0, -40, 27], [41, 0, 0], [0, 3, 40]])
                csl_rots['rots']['41f'] = np.array([[4, -33, 72], [39, -4, 36], [4, 8, 31]])
                csl_rots['rots']['41g'] = np.array([[31, -12, 72], [24, -4, -99], [4, 13, 4]])
                csl_rots['rots']['41h'] = np.array([[24, 4, 99], [31, 12, -72], [-4, 13, 4]])
                csl_rots['rots']['41i'] = np.array([[24, -23, 72], [32, 24, -27], [-3, 8, 32]])
                csl_rots['rots']['43a'] = np.array([[25, 18, 90], [18, 25, -90], [-10, 10, 7]])
                csl_rots['rots']['43b'] = np.array([[42, 2, 27], [2, 39, -54], [-3, 6, 38]])
                csl_rots['rots']['43c'] = np.array([[25, -30, 54], [30, 7, -90], [6, 10, 25]])
                csl_rots['rots']['43d'] = np.array([[38, -9, 54], [18, -2, -117], [3, 14, 2]])
                csl_rots['rots']['43e'] = np.array([[30, -7, 90], [25, 30, -54], [-6, 10, 25]])
                csl_rots['rots']['45a'] = np.array([[30, -30, 45], [33, 30, -18], [-2, 5, 42]])
                csl_rots['rots']['45b'] = np.array([[42, 6, 45], [6, 33, -90], [-5, 10, 30]])
                csl_rots['rots']['45c'] = np.array([[15, -30, 90], [42, 6, -45], [2, 11, 30]])
                csl_rots['rots']['45d'] = np.array([[33, -30, 18], [30, 30, -45], [2, 5, 42]])
                csl_rots['rots']['45e'] = np.array([[15, -30, 90], [42, 15, -18], [-2, 10, 33]])
                csl_rots['rots']['45f'] = np.array([[42, -15, 18], [15, 30, -90], [2, 10, 33]])
                csl_rots['rots']['45g'] = np.array([[33, -6, 90], [30, 15, -90], [-2, 14, 15]])
                csl_rots['rots']['45h'] = np.array([[42, -6, 45], [15, 30, -90], [-2, 11, 30]])
                csl_rots['rots']['45i'] = np.array([[30, -15, 90], [33, 6, -90], [2, 14, 15]])
                csl_rots['rots']['45j'] = np.array([[30, 6, 99], [30, 15, -90], [-5, 14, 6]])
                csl_rots['rots']['47a'] = np.array([[27, -34, 54], [38, 27, -18], [-2, 6, 43]])
                csl_rots['rots']['47b'] = np.array([[38, 27, 18], [27, 34, -54], [2, 6, 43]])
                csl_rots['rots']['47c'] = np.array([[38, 18, 63], [18, 11, -126], [7, 14, 2]])
                csl_rots['rots']['47d'] = np.array([[11, -42, 54], [42, 2, -63], [6, 7, 38]])
                csl_rots['rots']['47e'] = np.array([[43, 6, 54], [6, 38, -81], [-6, 9, 34]])
                csl_rots['rots']['47f'] = np.array([[6, -38, 81], [43, -6, -54], [6, 9, 34]])
                csl_rots['rots']['47g'] = np.array([[38, -21, 54], [21, 2, -126], [6, 14, 11]])
                csl_rots['rots']['47h'] = np.array([[18, 11, 126], [38, 18, -63], [-7, 14, -2]])
                csl_rots['rots']['47i'] = np.array([[42, -2, 63], [11, 42, -54], [-6, 7, 38]])
                csl_rots['rots']['47j'] = np.array([[21, -2, 126], [38, 21, -54], [-6, 14, 11]])
                csl_rots['rots']['49a'] = np.array([[36, -23, 72], [31, 36, -36], [-4, 8, 41]])
                csl_rots['rots']['49b'] = np.array([[41, 12, 72], [12, 31, -108], [-8, 12, 23]])
                csl_rots['rots']['49c'] = np.array([[48, 4, 27], [4, 33, -108], [-3, 12, 32]])
                csl_rots['rots']['49d'] = np.array([[40, -15, 72], [24, 40, -45], [-5, 8, 40]])
                csl_rots['rots']['49e'] = np.array([[41, -24, 36], [24, 23, -108], [4, 12, 31]])
                csl_rots['rots']['49f'] = np.array([[24, -40, 45], [40, 15, -72], [5, 8, 40]])
                csl_rots['rots']['49g'] = np.array([[24, -23, 108], [41, 24, -36], [-4, 12, 31]])


    if lat_elem == 'hP_ca':
        if sig_type == 'specific':
            lat_tau = lat1.lat_params['a']**2/lat1.lat_params['c']**2
            ### Hexagonal Specific
            csl_rots = {}
            csl_rots['lattice'] = lat1
            csl_rots['rots'] = {}
            tol1 = 1e-10

            if abs(lat_tau - 3.0/7.0 ) < tol1:
                ### c/a  = 1.528
                csl_rots['type'] = 'quads'                
                csl_rots['rots']['10'] = np.array([1, 2, 1, 0])
                csl_rots['rots']['16'] = np.array([3, 7, 0, 0])
                csl_rots['rots']['19'] = np.array([6, 14, 7, 0])

            if abs(lat_tau - 8.0/19.0 ) < tol1:
                ### c/a  = 1.541
                csl_rots['type'] = 'quads'
                csl_rots['rots']['21'] = np.array([2, 2, 1, 0])

            if abs(lat_tau - 5.0/12.0 ) < tol1:
                ### c/a  = 1.549
                csl_rots['type'] = 'quads'
                csl_rots['rots']['9'] = np.array([5, 12, 0, 0])
                csl_rots['rots']['12'] = np.array([5, 6, 0, 0])
                csl_rots['rots']['16'] = np.array([5, 12, 6, 0])
                csl_rots['rots']['17'] = np.array([1, 2, 1, 0])
                csl_rots['rots']['19'] = np.array([5, 8, 4, 0])
                csl_rots['rots']['21a'] = np.array([5, 3, 0, 0])
                csl_rots['rots']['21b'] = np.array([2, 3, 0, 0])

            if abs(lat_tau - 16.0/39.0 ) < tol1:
                ### c/a  = 1.561
                csl_rots['type'] = 'quads'
                csl_rots['rots']['17'] = np.array([2, 3, 0, 0])

            if abs(lat_tau - 11.0/27.0 ) < tol1:
                ### c/a  = 1.567
                csl_rots['type'] = 'quads'
                csl_rots['rots']['20'] = np.array([11, 27, 0, 0])

            if abs(lat_tau - 2.0/5.0 ) < tol1:
                ### c/a  = 1.581
                csl_rots['type'] = 'quads'
                csl_rots['rots']['7'] = np.array([1, 2, 1, 0])
                csl_rots['rots']['11a'] = np.array([2, 2, 1, 0])
                csl_rots['rots']['11b'] = np.array([2, 5, 0, 0])
                csl_rots['rots']['13a'] = np.array([2, 3, 0, 0])
                csl_rots['rots']['13b'] = np.array([4, 10, 5, 0])
                csl_rots['rots']['17a'] = np.array([1, 1, 0, 0])
                csl_rots['rots']['17b'] = np.array([2, 5, 1, 0])
                csl_rots['rots']['19a'] = np.array([3, 5, 0, 1])
                csl_rots['rots']['19b'] = np.array([2, 6, 3, 0])

            if abs(lat_tau - 7.0/18.0 ) < tol1:
                ### c/a  = 1.604
                csl_rots['type'] = 'quads'
                csl_rots['rots']['13'] = np.array([7, 18, 0, 0])
                csl_rots['rots']['17'] = np.array([7, 9, 0, 0])

            if abs(lat_tau - 5.0/13.0 ) < tol1:
                ### c/a  = 1.612
                csl_rots['type'] = 'quads'
                csl_rots['rots']['18'] = np.array([1, 2, 1, 0])

            if abs(lat_tau - 8.0/21.0 ) < tol1:
                ### c/a  = 1.620
                csl_rots['type'] = 'quads'
                csl_rots['rots']['9'] = np.array([2, 3, 0, 0])
                csl_rots['rots']['13'] = np.array([2, 6, 3, 0])
                csl_rots['rots']['15a'] = np.array([4, 3, 0, 0])
                csl_rots['rots']['15b'] = np.array([8, 21, 0, 0])
                csl_rots['rots']['17'] = np.array([4, 6, 3, 0])
                csl_rots['rots']['21a'] = np.array([4, 9, 3, 0])
                csl_rots['rots']['21b'] = np.array([6, 14, 7, 2])

            if abs(lat_tau - 3.0/8.0 ) < tol1:
                ### c/a  = 1.633
                csl_rots['type'] = 'quads'
                csl_rots['rots']['10'] = np.array([3, 8, 4, 0])
                csl_rots['rots']['11'] = np.array([1, 2, 1, 0])
                csl_rots['rots']['14'] = np.array([3, 4, 2, 0])
                csl_rots['rots']['17'] = np.array([3, 8, 0, 0])
                csl_rots['rots']['18'] = np.array([1, 2, 0, 0])

            if abs(lat_tau - 10.0/27.0 ) < tol1:
                ### c/a  = 1.643
                csl_rots['type'] = 'quads'
                csl_rots['rots']['19'] = np.array([10, 27, 0, 0])
                csl_rots['rots']['21'] = np.array([5, 9, 0, 0])

            if abs(lat_tau - 11.0/30.0 ) < tol1:
                ### c/a  = 1.651
                csl_rots['type'] = 'quads'
                csl_rots['rots']['21'] = np.array([11, 30, 0, 0])

            if abs(lat_tau - 4.0/11.0 ) < tol1:
                ### c/a  = 1.658
                csl_rots['type'] = 'quads'
                csl_rots['rots']['12'] = np.array([2, 2, 1, 0])
                csl_rots['rots']['14'] = np.array([2, 3, 0, 0])
                csl_rots['rots']['15'] = np.array([1, 2, 1, 0])
                csl_rots['rots']['18'] = np.array([2, 5, 1, 0])
                csl_rots['rots']['20'] = np.array([2, 6, 3, 0])

            if abs(lat_tau - 5.0/14.0 ) < tol1:
                ### c/a  = 1.673
                csl_rots['type'] = 'quads'
                csl_rots['rots']['17'] = np.array([5, 14, 7, 0])
                csl_rots['rots']['19'] = np.array([1, 2, 1, 0])


            if abs(lat_tau - 16.0/45.0 ) < tol1:
                ### c/a  = 1.677
                csl_rots['type'] = 'quads'
                csl_rots['rots']['16'] = np.array([4, 9, 0, 0])
                csl_rots['rots']['17'] = np.array([8, 15, 0, 0])
                csl_rots['rots']['19'] = np.array([2, 3, 0, 0])

    if lat_elem == 'hR_ca':
        if sig_type == 'specific':
            
            talpha = lat1.lat_params['alpha']
            lat_tau = np.cos(talpha)/(1 + 2*np.cos(talpha))

            ### Rhombohedral Specific
            csl_rots = {}
            csl_rots['lattice'] = lat1
            csl_rots['type'] = 'matrices'
            csl_rots['rots'] = {}
            tol1 = 1e-10

            if abs(lat_tau - 4.0/15.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.739
                ### x.y reads as x subscript y
                csl_rots['rots']['7a'] = np.array([1, 0, 2, 1])
                csl_rots['rots']['7b'] = np.array([3, 6, 0, 3])
                csl_rots['rots']['11a'] = np.array([1, 1, 0, 1])
                csl_rots['rots']['11b'] = np.array([3, 0, 3, 3])
                csl_rots['rots']['11c'] = np.array([1, 0, 3, 0])
                csl_rots['rots']['13a'] = np.array([2, 0, 3, 0])
                csl_rots['rots']['13b'] = np.array([6, 0, 15, 6])
                csl_rots['rots']['13c'] = np.array([2, 5, 0, 2])
                csl_rots['rots']['17a'] = np.array([3, 2, 2, 3])
                csl_rots['rots']['17b'] = np.array([1, 2, 1, 1])
                csl_rots['rots']['17c'] = np.array([3, 3, 6, 3])
                csl_rots['rots']['19a'] = np.array([3, 0, 5, 1])
                csl_rots['rots']['19b'] = np.array([9, 15, 0, 3])
                csl_rots['rots']['19c'] = np.array([2, 3, 3, 0])
                csl_rots['rots']['21'] = np.array([1, 1, 1, 0])
                csl_rots['rots']['23a'] = np.array([6, 5, 5, 0])
                csl_rots['rots']['23b'] = np.array([2, 0, 5, 1])
                csl_rots['rots']['23c'] = np.array([6, 15, 0, 3])
                csl_rots['rots']['23d'] = np.array([3, 4, 4, 3])

            if abs(lat_tau - 9.0/34.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.699
                csl_rots['rots']['8'] = np.array([1, 0, 2, 1])
                csl_rots['rots']['15'] = np.array([7, 0, 17, 7])
                csl_rots['rots']['22'] = np.array([21, 34, 0, 7])
                csl_rots['rots']['24'] = np.array([1, 1, 1, 0])
                csl_rots['rots']['25'] = np.array([1, 1, 0, 1])

            if abs(lat_tau - 31.0/117.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.704
                csl_rots['rots']['25'] = np.array([3, 0, 9, 3])

            if abs(lat_tau - 13.0/49.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.711
                csl_rots['rots']['14'] = np.array([5, 0, 7, 5])
                csl_rots['rots']['18'] = np.array([1, 1, 0, 1])
                csl_rots['rots']['21'] = np.array([5, 14, 0, 5])
                csl_rots['rots']['23'] = np.array([1, 0, 2, 1])
                csl_rots['rots']['28a'] = np.array([3, 7, 0, 1])
                csl_rots['rots']['28b'] = np.array([1, 2, 1, 1])

            if abs(lat_tau - 17.0/64.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.717
                csl_rots['rots']['14'] = np.array([13, 0, 32, 13])
                csl_rots['rots']['15'] = np.array([1, 0, 2, 1])
                csl_rots['rots']['20'] = np.array([13, 16, 0, 13])

            if abs(lat_tau - 25.0/94.0) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.724
                csl_rots['rots']['22'] = np.array([1, 0, 2, 1])

            if abs(lat_tau - 33.0/124.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.728
                csl_rots['rots']['27'] = np.array([25, 0, 62, 25])

            if abs(lat_tau - 31.0/116.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.750
                # csl_rots['rots']['25'] = np.array([25, 0, 58, 23])
                csl_rots['rots']['27'] = np.array([1, 0, 2, 1])

            if abs(lat_tau - 23.0/86.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.755
                csl_rots['rots']['20'] = np.array([1, 0, 2, 1])

            if abs(lat_tau - 19.0/71.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.758
                csl_rots['rots']['26'] = np.array([1, 1, 0, 1])

            if abs(lat_tau - 15.0/56.0 ) < tol1:
                csl_rots['type'] = 'quads'
                ### c/a  = 2.763
                csl_rots['rots']['12'] = np.array([11, 0, 28, 11])
                csl_rots['rots']['13'] = np.array([1, 0, 2, 1])
                csl_rots['rots']['17'] = np.array([11, 14, 0, 11])
                csl_rots['rots']['27'] = np.array([11, 8, 0, 11])

    return csl_rots


## Input parameters for pkl files
if test_case == 1:
    sig_type = 'common'
    l1 = lat.Lattice()
    csl_rots = lit_csl_rots(l1, sig_type)
    save_csl_rots(csl_rots, sig_type, l1)

elif test_case == 2:
    sig_type = 'common'
    lat_type = 'tP_Id'
    l1 = lat.Lattice(lat_type)
    csl_rots = lit_csl_rots(l1, sig_type)
    sig_rots = {}
    sig_rots[sig_type] = csl_rots
    save_csl_rots(sig_rots, sig_type, l1)

elif test_case == 3:
    sig_type = 'common'
    lat_type = 'hP_Id'
    l1 = lat.Lattice(lat_type)
    csl_rots = lit_csl_rots(l1, sig_type)
    sig_rots = {}
    sig_rots[sig_type] = csl_rots
    save_csl_rots(sig_rots, sig_type, l1)

elif test_case == 4:
    sig_type = 'common'
    lat_type = 'hR_Id'
    l1 = lat.Lattice(lat_type)
    csl_rots = lit_csl_rots(l1, sig_type)
    sig_rots = {}
    sig_rots[sig_type] = csl_rots
    save_csl_rots(sig_rots, sig_type, l1)

elif test_case == 5:
    lat_tau = []
    lat_tau.append(1.0/9.0)
    ca_rat = []
    sig_type = 'specific'
    sig_rots = {}
    for tau_vals in lat_tau:
        ca_rat = np.sqrt(1/tau_vals)
        lat_type = 'tP_ca'
        l1 = lat.Lattice(lat_type, ca_rat)
        csl_rots = lit_csl_rots(l1, sig_type)
        [n1, d1] = int_man.rat(tau_vals)
        nu = n1[0][0]
        mu = d1[0][0]
        tau_str = str(nu) + '_' + str(mu)
        sig_rots[tau_str] = csl_rots
        save_csl_rots(sig_rots, sig_type, l1)

elif test_case == 6:
    lat_tau = []
    lat_tau.append(3.0/7.0)
    lat_tau.append(8.0/19.0)
    lat_tau.append(5.0/12.0)
    lat_tau.append(16.0/39.0)
    lat_tau.append(11.0/27.0)
    lat_tau.append(2.0/5.0)
    lat_tau.append(7.0/18.0)
    lat_tau.append(5.0/13.0)
    lat_tau.append(8.0/21.0)
    lat_tau.append(3.0/8.0)
    lat_tau.append(10.0/27.0)
    lat_tau.append(11.0/30.0)
    lat_tau.append(4.0/11.0)
    lat_tau.append(5.0/14.0)
    lat_tau.append(16.0/45.0)
    ca_rat = []
    sig_type = 'specific'
    sig_rots = {}
    for tau_vals in lat_tau:
        ca_rat = np.sqrt(1/tau_vals)
        lat_type = 'hP_ca'
        l1 = lat.Lattice(lat_type, ca_rat)
        csl_rots = lit_csl_rots(l1, sig_type)

        [n1, d1] = int_man.rat(tau_vals)
        nu = n1[0][0]
        mu = d1[0][0]
        tau_str = str(nu) + '_' + str(mu)
        sig_rots[tau_str] = csl_rots

    save_csl_rots(sig_rots, sig_type, l1)

elif test_case == 7:
    lat_tau = []
    lat_tau.append(4.0/15.0)
    lat_tau.append(9.0/34.0 )
    lat_tau.append(31.0/117.0 )
    lat_tau.append(13.0/49.0 )
    lat_tau.append(17.0/64.0 )
    lat_tau.append(25.0/94.0)
    lat_tau.append(33.0/124.0 )
    lat_tau.append(31.0/116.0 )
    lat_tau.append(23.0/86.0 )
    lat_tau.append(19.0/71.0 )
    lat_tau.append(15.0/56.0 )

    ca_rat = []
    sig_type = 'specific'
    sig_rots = {}
    for tau_vals in lat_tau:
        ca_rat = np.sqrt(3/(2 - 6*tau_vals))
        lat_type = 'hR_ca'
        l1 = lat.Lattice(lat_type, ca_rat)
        csl_rots = lit_csl_rots(l1, sig_type)

        [n1, d1] = int_man.rat(tau_vals)
        nu = n1[0][0]
        mu = d1[0][0]
        tau_str = str(nu) + '_' + str(mu)
        sig_rots[tau_str] = csl_rots

    save_csl_rots(sig_rots, sig_type, l1)



