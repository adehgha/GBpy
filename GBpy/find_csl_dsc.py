# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import integer_manipulations as int_man
from tools import lll_reduction
from tools import smith_nf
from tools import message_display
import bp_basis as bpb
# -----------------------------------------------------------------------------------------------------------


def sigma_calc(t_g1tog2_g1):
    """
    Computes the sigma of the transformation matrix

    * if det(T) = det(T^{-1}) then sigma1 = sigma2 is returned
    * if det(T) \\neq det(T^{-1}) then max(sigma1, sigma2) is returned

    """

    R = np.array(t_g1tog2_g1)
    R2 = np.linalg.det(R)*np.linalg.inv(R)
    n1, d1 = int_man.rat(R)
    n2, d2 = int_man.rat(R2)
    # -----------------------------
    Sigma21 = int_man.lcm_array(d1[:])
    # -----------------------------
    Sigma22 = int_man.lcm_array(d2[:])
    Sigma = np.array([Sigma21, Sigma22]).max()
    return Sigma
# -----------------------------------------------------------------------------------------------------------


def reciprocal_mat(l_g_go):
    """
    The reciprocal matrix with reciprocal basis vectors is computed for the
    input matrix with primitve basis vectors

    Parameters
    ----------------
    l_g_go: numpy array
        The primitive basis vectors b1x, b1y, b1z

    Returns
    -----------
    rl_g_go: numpy array
        The primitve reciprocal basis vectors
    """
    InMat = np.array(l_g_go)

    L3 = np.cross(InMat[:, 0], InMat[:, 1]) / np.linalg.det(InMat)
    L1 = np.cross(InMat[:, 1], InMat[:, 2]) / np.linalg.det(InMat)
    L2 = np.cross(InMat[:, 2], InMat[:, 0]) / np.linalg.det(InMat)
    rl_g_go = np.vstack((L1, L2, L3)).T
    return rl_g_go
# -----------------------------------------------------------------------------------------------------------


def csl_elem_div_thm_l1(T0, l_g1n_g1):
    """
    The csl basis vectors are obtained from the diagonal matrix using the
    algorithm specified in doi:10.1107/S056773947601231X. There are two
    algorithms specified based on numerators or denominators of the T0 matrix.
    The numerators are used in this function.

    Parameters
    ---------------
    T0: numpy array
        The transformation matrix in G1n reference frame

    l_g1n_g1: numpy array
        The 'new' basis vectors of g1 lattice (g1n) in g1 reference frame

    Returns
    ------------
    l_csl_g1: numpy array
        The CSL basis vectors in g1 reference frame
    """
    T0 = np.array(T0)
    L1 = np.array(l_g1n_g1)

    if T0.shape[0] == 3:
        p1 = int_man.rat(np.array(T0[0, 0]), 1e-06)[0][0][0]
        p2 = int_man.rat(np.array(T0[1, 1]), 1e-06)[0][0][0]
        p3 = int_man.rat(np.array(T0[2, 2]), 1e-06)[0][0][0]
        l_csl_g1 = np.dot(L1, np.array([[p1, 0, 0], [0, p2, 0], [0, 0, p3]]))
    elif T0.shape[0] == 2:
        p1 = int_man.rat(np.array(T0[0, 0]), 1e-06)[0][0][0]
        p2 = int_man.rat(np.array(T0[1, 1]), 1e-06)[0][0][0]
        l_csl_g1 = np.dot(L1, np.array([[p1, 0, 0], [0, p2, 0]]))
    return l_csl_g1
# -----------------------------------------------------------------------------------------------------------


def csl_elem_div_thm_l2(t0, l_g2n_g2):
    """
    The csl basis vectors are obtained from the diagonal matrix using the
    algorithm specified in doi:10.1107/S056773947601231X. There are two
    algorithms specified based on numerators or denominators of the T0 matrix.
    The denominators are used in this function.

    Parameters
    ---------------
    T0: numpy array
        The transformation matrix in G1n reference frame

    l_g2n_g2: numpy array
        The 'new' basis vectors of g2 lattice (g2n) in g2 reference frame

    Returns
    ------------
    l_csl_g2: numpy array
        The CSL basis vectors in g2 reference frame
    """
    t0 = np.array(t0)
    l2 = np.array(l_g2n_g2)

    if t0.shape[0] == 3:
        q1 = int_man.rat(np.array(t0[0, 0]), 1e-06)[1][0][0]
        q2 = int_man.rat(np.array(t0[1, 1]), 1e-06)[1][0][0]
        q3 = int_man.rat(np.array(t0[2, 2]), 1e-06)[1][0][0]
        l_csl_g2 = np.dot(l2, np.array([[q1, 0, 0], [0, q2, 0], [0, 0, q3]]))
    elif t0.shape[0] == 2:
        q1 = int_man.rat(np.array(t0[0, 0]), 1e-06)[1][0][0]
        q2 = int_man.rat(np.array(t0[1, 1]), 1e-06)[1][0][0]
        l_csl_g2 = np.dot(l2, np.array([[q1, 0, 0], [0, q2, 0]]))
    return l_csl_g2
# -----------------------------------------------------------------------------------------------------------


def csl_finder_smith(r_g1tog2_g1):
    """
    This funciton extracts the CSL basis when transformation between the two
    lattices is given (r_g1tog2_g1). The algorithms used are based on the
    following article: doi:10.1107/S056773947601231X)

    Parameters
    ----------------
    r_g1tog2_g1: numpy array
        The 3x3 transformation matrix in g1 reference frame

    Returns
    -----------
    l_csl_g1: numpy array
        3 x 3 matrix with the csl basis vectors as columns

    Notes
    ---------
    The "Reduced" refer to the use of LLL algorithm to compute a
    basis that is as close to orthogonal as possible.
    (Refer to http://en.wikipedia.org/wiki/Lattice_reduction) for further
    detials on the concept of Lattice Reduction
    """
    R_G1ToG2_G1 = np.array(r_g1tog2_g1)
    L_G2_G1 = R_G1ToG2_G1

    # Obtain N1 and N2
    N1, _ = int_man.int_mult(L_G2_G1)
    A = np.dot(N1, L_G2_G1)

    # Check that A is an integer matrix
    if int_man.int_check(A, 12).all():
        A = np.around(A)
    else:
        raise Exception('L_G2_G1 is not a Sigma Transformation')
    # Smith Normal Form of A
    # A = U*E*V'
    U, E, _ = smith_nf(A)
    L_G1n_G1 = U
    # L_G1_G1n = np.linalg.inv(L_G1n_G1)

    # CSL Bases
    l_csl_g1 = csl_elem_div_thm_l1(E / N1, L_G1n_G1)

    if(int_man.int_check(l_csl_g1, 1e-08)).all():
        l_csl_g1 = np.around(l_csl_g1)
    else:
        raise Exception('l_csl_g1 is not defined in L_G1_G1 axis')

    # Reduced CSL bases using the LLL algorithm
    # Actually don't reduce yet because the matrix is in "g1" reference frame!
    # l_csl_g1 = lll_reduction((l_csl_g1))
    return l_csl_g1
# -----------------------------------------------------------------------------------------------------------


def check_csl_finder_smith(r_g1tog2_g1, Sigma, L_G1_GO1, L_CSL_G1):
    """
    This function checks the obtained CSL basis vectors are correct by
    using the following conditions:
    * The CSL basis vectors are integer combinations of basis vectors of
    lattice 1
    * The CSL basis vectors are integer combinations of basis vectors of
    lattice 2
    * The volume enclosed by the CSL is sigma times the volume of lattice 1
    """
    R_G1ToG2_G1 = np.array(r_g1tog2_g1)
    Sigma = np.array(Sigma)
    L_G1_GO1 = np.array(L_G1_GO1)
    L_CSL_G1 = np.array(L_CSL_G1)

    L_GO1_G1 = np.linalg.inv(L_G1_GO1)
    R_G1ToG2_GO1 = np.dot(np.dot(L_G1_GO1, R_G1ToG2_G1), L_GO1_G1)
    L_CSL_GO1 = np.dot(L_G1_GO1, L_CSL_G1)
    L_G2_GO1 = np.dot(R_G1ToG2_GO1, L_G1_GO1)

    print '*** CSL checks ***'
    # -----Check-1: L_CSL_GO1 is defined in the L_G1_GO1 lattice
    CheckBase1 = np.dot(L_GO1_G1, L_CSL_GO1)
    Precis = 10
    message_display(CheckBase1, 1,
                    'L_CSL_GO1 is defined in the L_G1_GO1 lattice', Precis)

    # -----Check-2: L_CSL_GO1 is defined in the L_G2_GO1 lattice
    CheckBase2 = np.dot(np.linalg.inv(L_G2_GO1), (L_CSL_GO1))
    message_display(CheckBase2, 2,
                    'L_CSL_GO1 is defined in the L_G2_GO1 lattice', Precis)

    # -----Check-3: Check that we have the right volume for L_CSL_GO1
    CheckBase3 = np.linalg.det(L_CSL_GO1) / (Sigma * np.linalg.det(L_G1_GO1))
    Disp_str = ('V(CSL_GO1)/V(G1_GO1) = Sigma =  ' + "%0.0f"
                % (np.linalg.det(L_CSL_GO1) / np.linalg.det(L_G1_GO1)))
    message_display(CheckBase3, 3, Disp_str, Precis)
# -----------------------------------------------------------------------------------------------------------


def dsc_finder(L_G2_G1, L_G1_GO1):
    """
    The DSC lattice is computed for the bi-crystal, if the transformation
    matrix l_g2_g1 is given and the basis vectors of the underlying crystal
    l_g_go (in the orthogonal reference go frame) are known. The following
    relationship is used: **The reciprocal of the coincidence site lattice of
    the reciprocal lattices is the DSC lattice**

    Parameters
    ----------------
    l_g2_g1: numpy array
        transformation matrix (r_g1tog2_g1)

    l_g1_go1: numpy array
        basis vectors (as columns) of the underlying lattice expressed in the
        orthogonal 'go' reference frame

    Returns
    ------------
    l_dsc_g1: numpy array
        The dsc lattice basis vectors (as columns) expressed in the g1 reference

    Notes
    ---------
    The "Reduced" refer to the use of LLL algorithm to compute a
    basis that is as close to orthogonal as possible.
    (Refer to http://en.wikipedia.org/wiki/Lattice_reduction) for further
    detials on the concept of Lattice Reduction
    """

    L_G2_G1 = np.array(L_G2_G1)
    L_G1_GO1 = np.array(L_G1_GO1)

    L_GO1_G1 = np.linalg.inv(L_G1_GO1)
    # % % Reciprocal lattice of G1
    # --------------------------------------------------------------
    L_rG1_GO1 = reciprocal_mat(L_G1_GO1)
    L_GO1_rG1 = np.linalg.inv(L_rG1_GO1)
    # L_rG1_G1 = np.dot(L_GO1_G1, L_rG1_GO1)
    # % % L_G1_rG1 = L_rG1_G1^(-1);
    # % % Reciprocal lattice of G2
    L_G2_GO1 = np.dot(L_G1_GO1, L_G2_G1)
    L_rG2_GO1 = reciprocal_mat(L_G2_GO1)

    # % % Transformation of the Reciprocal lattices
    # % % R_rG1TorG2_rG1 = L_rG2_G1*L_G1_rG1;
    L_rG2_rG1 = np.dot(L_GO1_rG1, L_rG2_GO1)
    Sigma_star = sigma_calc(L_rG2_rG1)

    # % % CSL of the reciprocal lattices
    L_rCSL_rG1 = csl_finder_smith(L_rG2_rG1)
    L_rCSL_GO1 = np.dot(L_rG1_GO1, L_rCSL_rG1)
    # % % Reciprocal of the CSL of the reciprocal lattices
    L_DSC_GO1 = reciprocal_mat(L_rCSL_GO1)
    L_DSC_G1 = np.dot(L_GO1_G1, L_DSC_GO1)
    # L_DSC_G1 = bpb.reduce_go1_mat(L_DSC_G1, L_G1_GO1)

    # # % % Reduction of the DSC lattice in G1 reference frame
    # DSC_Int = int_man.int_finder(L_DSC_G1, 1e-06)
    # t_ind = np.where(abs(DSC_Int) == abs(DSC_Int).max())
    # t_ind_1 = t_ind[0][0]
    # t_ind_2 = t_ind[1][0]
    # Mult1 = DSC_Int[t_ind_1, t_ind_2] / L_DSC_G1[t_ind_1, t_ind_2]
    # DSC_Reduced = lll_reduction(DSC_Int)
    # DSC_Reduced = DSC_Reduced / Mult1
    # L_DSC_G1 = DSC_Reduced

    # % % % Check this assertion: L_DSC_G1 = [Int_Matrix]/Sigma
    if int_man.int_check(Sigma_star*L_DSC_G1, 10).all():
        L_DSC_G1 = np.around(Sigma_star*L_DSC_G1) / Sigma_star
    else:
        raise Exception('L_DSC_G1 is not equal to [Int_Matrix]/Sigma')
    return L_DSC_G1
# -----------------------------------------------------------------------------------------------------------


def check_dsc_finder(R_G1ToG2_G1, Sigma, L_G1_GO1, L_DSC_G1, L_CSL_G1):
    """
    This function checks the obtained DSC basis vectors are correct by
    using the following conditions:
    * Lattice 1 basis vectors are integer combinations of basis vectors of
    the DSC lattice
    * Lattice 2 basis vectors are integer combinations of basis vectors of
    the DSC lattice
    * The volume enclosed by the DSC is 1/sigma times the volume of lattice 1
    """
    L_GO1_G1 = np.linalg.inv(L_G1_GO1)
    R_G1ToG2_GO1 = np.dot(np.dot(L_G1_GO1, R_G1ToG2_G1), L_GO1_G1)
    L_DSC_GO1 = np.dot(L_G1_GO1, L_DSC_G1)
    L_CSL_GO1 = np.dot(L_G1_GO1, L_CSL_G1)
    L_G2_GO1 = np.dot(R_G1ToG2_GO1, L_G1_GO1)

    print '*** DSC checks ***'
    # -----Check-1:
    Check1 = np.dot(np.linalg.inv(L_DSC_GO1), L_G1_GO1)
    Precis = 10
    message_display(Check1, 1,
                    'L_G1_GO1 is defined in the obtained DSC basis', Precis)

    # -----Check-2:
    Check2 = np.dot(np.linalg.inv(L_DSC_GO1), L_G2_GO1)
    message_display(Check2, 2,
                    'L_G2_GO1 is defined in the obtained DSC basis', Precis)

    # -----Check-3:
    Check3 = np.dot(np.linalg.inv(L_DSC_GO1), L_CSL_GO1)
    message_display(Check3, 3,
                    'L_CSL_GO1 is defined in the obtained DSC basis', Precis)

    # -----Check-4:
    CheckBase4 = np.linalg.det(L_DSC_GO1) * Sigma / np.linalg.det(L_G1_GO1)
    Disp_str = ('V(G1_GO1)/V(DSC_GO1) = Sigma =  ' + "%0.0f"
                % (np.linalg.det(L_G1_GO1) / np.linalg.det(L_DSC_GO1)))
    Precis = 6
    message_display(CheckBase4, 4, Disp_str, Precis)
# -----------------------------------------------------------------------------------------------------------


def find_csl_dsc(L_G1_GO1, R_G1ToG2_G1):
    """
    This function calls the csl_finder and dsc_finder and returns
    the CSL and DSC basis vectors in 'g1' reference frame.

    Parameters
    -----------------
    L_G1_GO1: numpy array
        The three basis vectors for the primitive unit cell
        (as columns) are given with respect to the GO1 reference
        frame.

    R_G1ToG2_G1: 3X3 numpy array
        The rotation matrix defining the
        transformation in 'G1' reference frame. The subscript 'G1' refers
        to the primitive unit cell of G lattice.

    Returns
    l_csl_g1, l_dsc_g1: numpy arrays
        The basis vectors of csl and dsc lattices in the g1 reference frame
    """

    R_G1ToG2_G1 = np.array(R_G1ToG2_G1)
    L_G1_GO1 = np.array(L_G1_GO1)

    L_CSL_G1 = csl_finder_smith(R_G1ToG2_G1)
    print np.dot(np.linalg.inv(R_G1ToG2_G1), L_CSL_G1)

    L_CSL_G1 = bpb.reduce_go1_mat(L_CSL_G1, L_G1_GO1)

    Sigma, _ =  int_man.int_mult(R_G1ToG2_G1)
    check_csl_finder_smith(R_G1ToG2_G1, Sigma, L_G1_GO1, L_CSL_G1)

    # Finding the DSC lattice from the obtained CSL.
    L_DSC_G1 = dsc_finder(R_G1ToG2_G1, L_G1_GO1)

    L_CSL_G1 = make_right_handed(L_CSL_G1, L_G1_GO1)
    L_DSC_G1 = make_right_handed(L_DSC_G1, L_G1_GO1)

    return L_CSL_G1, L_DSC_G1
# -----------------------------------------------------------------------------------------------------------

def make_right_handed(l_csl_p1, l_p_po):
    if (np.linalg.det(l_csl_p1) < 0):
        t1_array = l_csl_p1.copy()
        t1_array[:, 0] = l_csl_p1[:, 1]
        t1_array[:, 1] = l_csl_p1[:, 0]
        l_csl_p1 = t1_array.copy()
    return l_csl_p1
