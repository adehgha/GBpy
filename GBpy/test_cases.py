# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import integer_manipulations as int_man
import os
import sys
import find_csl_dsc as fcd
import lattice as lat
import bp_basis as bpb
from tools import smith_nf


def test_int_mult():
    """
    Test cases for the int_mult function in integer_manipulations
    """
    Mat = np.zeros(5, dtype=[('Matrix', '(3,3)float64')])
    Mat['Matrix'][0] = np.array(([1, 2.0e-5, 0], [-4, 5, 6], [7, 8, 15]))
    Mat['Matrix'][1] = np.array(([-1, 2.5, 3], [-4, 0, 6.5], [1./3, 0, 0]))
    Mat['Matrix'][2] = np.array(([1, 2./5, 0], [-4, 5, 1.0/6], [3, 0, -1]))
    Mat['Matrix'][3] = np.array(([1, 2, 3], [-4, 5.5, 6], [2, -2, 1]))
    Mat['Matrix'][4] = np.array(([1, 2, 25], [-4, 5, 6], [3, -2, 1e-7]))
    for i in range(Mat.shape[0]):
        a, b = int_man.int_mult(Mat['Matrix'][i])
        print (Mat['Matrix'][i], '\nIntegral Form: \n',
               b, '\nMultiplier:\n', a, '\n-------\n')
# -----------------------------------------------------------------------------------------------------------


def test_int_finder():
    """
    Test cases for the int_finder function in integer_manipulations
    """
    Mat = np.zeros(3, dtype=[('Matrix', '(3,3)float64'),
                             ('row', '(1,3)float64'), ('col', '(3,1)float64')])
    Mat['Matrix'][0] = np.array(([1.5, 2, 3.0e-7],
                                 [-4, 5, 6], [7, 2.0e-5, 1.0e-5]))
    Mat['Matrix'][1] = np.array(([1.5, 2, 3e-6],
                                 [-4, 0, 6], [3.5, -2e-6, 1e-6]))
    Mat['Matrix'][2] = np.array(([1.5, 2, 3e-7],
                                 [-4, 5, 6], [3.5, -2e-6, 1e-6]))
    Mat['row'][0] = np.array(([1.5e-7, 0, 3e-6]))
    Mat['row'][1] = np.array(([1.5e-7, 1, 3e-6]))
    Mat['row'][2] = np.array(([1.5, 4, 3e-7]))
    Mat['col'][0] = np.array(([1.5e-7], [0], [3e-6]))
    Mat['col'][1] = np.array(([1.5e-7], [1], [3e-6]))
    Mat['col'][2] = np.array(([1.5], [4], [3e-7]))
    order = np.zeros(2, 'a4')
    order[0] = 'rows'
    order[1] = 'col'
    # order[2] = 'cols'
    cnt = 0
    for j in Mat.dtype.names:
        for i in range(Mat.shape[0]):
            for k in range(len(order)):
                cnt += 1
                print 'case:', cnt, '\n'
                print Mat[j][i], '\n\n', 'order:', order[k], '\nanswer:\n'

                # a = int_man.int_finder(Mat[j][i], tolerance, order[k])
                a = int_man.int_finder(Mat[j][i])
                print a, '\n', '--'
                a = int_man.int_finder(Mat[j][i], 1.0e-5, 'rows', 1.0e-5)
                print a, '\n', '--'
                a = int_man.int_finder(Mat[j][i], 1e-5, 'rows')
                print a, '\n', '--'
                a = int_man.int_finder(Mat[j][i], 1e-5, 'col')
                print a, '\n', '--'
                a = int_man.int_finder(Mat[j][i], 1e-5, 'columns')
                print a, '\n', '--'
                a = int_man.int_finder(Mat[j][i], 1e-5, order[k], 1e-5)
                print a, '\n', '--'

                print '\n', '-----------------------------------------'

    print cnt, ' test cases have been tried.'
    if __name__ == '__main__':
        test_int_finder
        # unittest.main()
# -----------------------------------------------------------------------------------------------------------


def test_int_check():
    """
    Test cases for the int_check function in integer_manipulations
    """
    b = np.array([[2, 3, 5], [6, 6.000002, -2.000001], [-0.00002, 1.5, 4]])
    a = int_man.int_check(b)
    print a
    # ------------
    b = 2.5
    a = int_man.int_check(b)
    print a
    if __name__ == '__main__':
        test_int_check

# -----------------------------------------------------------------------------------------------------------


def test_csl_finder_smith():
    """
    Test cases for the csl_finder find_csl_dsc
    """

    Mat = np.zeros(5, dtype=[('Matrix', '(3,3)float64')])
    # Mat['Matrix'][0] = np.array(([1, 2, 0], [-4, 5, 6], [7, 8, 15]))
    Mat['Matrix'][0] = np.array(([7./37, 42./37, 12./37],
                                 [30./37, -5./37, -12./37],
                                 [-42./37, -30./37, -35./37]))
    Mat['Matrix'][1] = np.array(([-1, 2, 3], [-4, 0, 6], [3, 0, 0]))
    Mat['Matrix'][2] = np.array(([1, 2, 0], [-4, 5, 6], [3, 0, -1]))
    Mat['Matrix'][3] = np.array(([1, 2, 3], [-4, 5, 6], [2, -2, 1]))
    Mat['Matrix'][4] = np.array(([1, 2, 25], [-4, 5, 6], [3, -2, 0]))
    for i in range(Mat.shape[0]):
        a = fcd.csl_finder_smith(Mat['Matrix'][i])
        print Mat['Matrix'][i], '\n CSL: \n', a, '\n-------\n'
# -----------------------------------------------------------------------------------------------------------


def test_smith_nf():
    """
    Test for Smith Normal Form computations from smith_nf
    """
    # a, b ,g= extgcd(15, -5)
    # print a, b, g
    import numpy as np
    D = np.array([[7, 42, 12], [30, -5, -12], [-42, -30, -35]])
    # print D
    # index0, index1 = np.nonzero(D)
    # print (index1)
    # # print min(np.nonzero(D != 0))
    U, S, V = smith_nf(D)
    print '\n---\nU=\n', U, '\n---\nS=\n', S, '\n---\nV=\n', V
# -----------------------------------------------------------------------------------------------------------


def test_csl_elem_div_thm_l1():
    """
    Test the csl_elem function from find_csl_dsc
    """
    import numpy as np
    a = np.array([[2, 3, 0], [4, -2, 7], [0, 2, 8]])
    print a, '\n'
    b = fcd.csl_elem_div_thm_l1(a, 2.5)
    print b, '\n--------\n'
    a = np.array([[2, 3], [4, -2]])
    print a, '\n'
    b = fcd.csl_elem_div_thm_l1(a, 2.5)
    print b
# -----------------------------------------------------------------------------------------------------------


def test_gb_2d_csl():
    """
    Test the two-dimensional boundary plane basis in bp_basis function
    """
    Mat = np.zeros(6, dtype=[('t_g1tog2_go1',
                              '(3,3)float64'), ('bp1_go1', '(3,1)float64')])
    Mat['bp1_go1'][0] = np.array([[0.5], [0.], [-0.5]])
    Mat['bp1_go1'][1] = np.array([[1.5], [-0.5], [-1.0]])
    Mat['bp1_go1'][2] = np.array([[11.], [-4.0], [-4.0]])
    Mat['bp1_go1'][3] = np.array([[12.5], [0.5], [-1.0]])
    Mat['bp1_go1'][4] = np.array([[3.5], [2.0], [0.5]])
    Mat['bp1_go1'][5] = np.array([[3.], [-0.5], [-0.5]])
    Mat['t_g1tog2_go1'][0] = np.array([[2./3, -1./3, 2./3],
                                       [2./3, 2./3, -1./3],
                                       [-1./3, 2./3, 2./3]])
    Mat['t_g1tog2_go1'][1] = np.array([[2./3, -1./3, 2./3],
                                       [2./3, 2./3, -1./3],
                                       [-1./3, 2./3, 2./3]])
    Mat['t_g1tog2_go1'][2] = np.array([[1./3, 2./3, 2./3],
                                       [-2./3, 2./3, -1./3],
                                       [-2./3, -1./3, 2./3]])
    Mat['t_g1tog2_go1'][3] = np.array([[-2./3, -1./3, 2./3],
                                       [1./3, 2./3, 2./3],
                                       [-2./3, 2./3, -1./3]])
    Mat['t_g1tog2_go1'][4] = np.array([[-2./3, -2./3, 1./3],
                                       [-2./3, 1./3, -2./3],
                                       [1./3, -2./3, -2./3]])
    Mat['t_g1tog2_go1'][5] = np.array([[1./3, 2./3, 2./3],
                                       [-2./3, 2./3, -1./3],
                                       [-2./3, -1./3, 2./3]])

    for i in range(Mat.shape[0]):
        AL = lat.Lattice('Al')
        bp1_go1 = Mat['bp1_go1'][i]
        t_g1tog2_go1 = Mat['t_g1tog2_go1'][i]
        a, b, c = bpb.bicryst_planar_den(bp1_go1, t_g1tog2_go1, AL)
        print ('Pl Density 1=', a, '\nPl Density 2=',
               b, '\nPl Density_2D CSL=', c)
        print '\n------------\n'
test_gb_2d_csl()
# -----------------------------------------------------------------------------------------------------------


def test_dsc_finder():
    """
    contains a number of random test cases for the dsc_finder algorithm in
    find_csl_dsc
    """

    Mat = np.zeros(1, dtype=[('Matrix', '(3,3)float64')])
    # Mat['Matrix'][0] = np.array(([1, 2, 0], [-4, 5, 6], [7, 8, 15]))
    # Mat['Matrix'][0] = np.array(([2./3, -1./2, 2./3],
    # [2./3, 2./3, -1./3],
    # [-1./3, 2./3, 2./3]))

    Mat['Matrix'][0] = np.array(([7./37, 42./37, 12./37],
                                 [30./37, -5./37, -12./37],
                                 [-42./37, -30./37, -35./37]))

    # Mat['Matrix'][1] = np.array(([-1, 2, 3], [-4, 0, 6], [3, 0, 0]))
    # Mat['Matrix'][2] = np.array(([1, 2, 0], [-4, 5, 6], [3, 0, -1]))
    # Mat['Matrix'][3] = np.array(([1, 2, 3], [-4, 5, 6], [2, -2, 1]))
    # Mat['Matrix'][4] = np.array(([1, 2, 25], [-4, 5, 6], [3, -2, 0]))

    b1x = np.array([[0.], [1.], [1.]])/2
    b1y = np.array([[1.], [0.], [1.]])/2
    b1z = np.array([[1.], [1.], [0.]])/2
    L_G1_GO1 = np.concatenate((b1x, b1y, b1z), axis=1)

    for i in range(Mat.shape[0]):
        a = fcd.dsc_finder(Mat['Matrix'][i], L_G1_GO1)
        print ('\n R_G1ToG2_G1: \n',
               Mat['Matrix'][i], '\n DSC: \n', a, '\n-------\n')
        # txt=Col()
        # txt.c_prnt('R_G1ToG2_G1', 'yel')
        # print Mat['Matrix'][i]
        # txt.c_prnt('DSC', 'dgrn')
        # print a, '\n-------\n'
# -----------------------------------------------------------------------------------------------------------
