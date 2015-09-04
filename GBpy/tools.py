# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import os, sys, inspect
import integer_manipulations as int_man
# -----------------------------------------------------------------------------------------------------------


class Col(object):
    """
    This class is defined to ouput a word or sentence in a different color
    to the standard shell.
    The colors available are:
    ``pink``, ``blue``, ``green``, ``dgrn``: dark green, ``yellow``, ``amber``
    """
    def __init__(self):
        self.pink = '\033[95m'
        self.blue = '\033[94m'
        self.green = '\033[92m'
        self.dgrn = '\033[1;32m'
        self.yellow = '\033[93m'
        self.amber = '\033[91m'
        self.ENDC = '\033[0m'

    def c_prnt(self, text, color):
        """
        Print a string in color,

        Parameters
        ----------
        Col: Col class instance
            an instance of the ``Col`` class

        text: string
            Text to be shown in color.

        color: string

        Returns
        --------
        N/A
        """
        if color is 'pink':
            a = self.pink
        elif color is 'blue':
            a = self.blue
        elif color is 'green':
            a = self.green
        elif color is 'dgrn':
            a = self.dgrn
        elif color is 'yel':
            a = self.yellow
        elif color is 'amber':
            a = self.amber
        else:
            raise Exception('The color you selected is not acceptable')
        print a + text + self.ENDC
# -----------------------------------------------------------------------------------------------------------


def lll_reduction(matrix, delta=0.75):
    """
    _________________________________________________________________________
    This function has been borrowed from pymatgen project with slight changes
    in input and output handling.
    Link to the project:
        http://pymatgen.org/
    Link to the source code:
        https://github.com/materialsproject/pymatgen/blob/92ee88ab6a6ec6e27b717150931e6d484d37a4e6/pymatgen/core/lattice.py
    _________________________________________________________________________
    Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
    c-reduced basis. This method returns a basis which is as "good" as
    possible, with "good" defined by orthongonality of the lattice vectors.

    Parameters
    ----------
    delta: float
        Reduction parameter. \v
        Default of 0.75 is usually fine.

    Returns
    -------
    a: numpy array
        Reduced lattice:
    """

    if int_man.int_check(matrix, 10).all():

        # Transpose the lattice matrix first so that basis vectors are columns.
        # Makes life easier.
        a = np.array(matrix)
        sz = a.shape

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coeffieicnts
        m = np.zeros(3)  # These are the norm squared of each vec.

        b[:, 0] = a[:, 0]
        m[0] = np.dot(b[:, 0], b[:, 0])
        for i in xrange(1, sz[1]):
            u[i, 0:i] = np.dot(a[:, i].T, b[:, 0:i]) / m[0:i]
            b[:, i] = a[:, i] - np.dot(b[:, 0:i], u[i, 0:i].T)
            m[i] = np.dot(b[:, i], b[:, i])

        k = 2

        while k <= sz[1]:
            # Size reduction.
            for i in xrange(k - 1, 0, -1):
                q = round(u[k - 1, i - 1])
                if q != 0:
                    # Reduce the k-th basis vector.
                    a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                    uu = list(u[i - 1, 0:(i - 1)])
                    uu.append(1)
                    # Update the GS coefficients.
                    u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

            # Check the Lovasz condition.
            if np.dot(b[:, k - 1], b[:, k - 1]) >=\
                    (delta - abs(u[k - 1, k - 2]) ** 2) *\
                    np.dot(b[:, (k - 2)], b[:, (k - 2)]):
                # Increment k if the Lovasz condition holds.
                k += 1
            else:
                # If the Lovasz condition fails,
                # swap the k-th and (k-1)-th basis vector
                v = a[:, k - 1].copy()
                a[:, k - 1] = a[:, k - 2].copy()
                a[:, k - 2] = v
                # Update the Gram-Schmidt coefficients
                for s in xrange(k - 1, k + 1):
                    u[s - 1, 0:(s - 1)] = np.dot(a[:, s - 1].T,
                                            b[:, 0:(s - 1)]) / m[0:(s - 1)]
                    b[:, s - 1] = a[:, s - 1] - np.dot(b[:, 0:(s - 1)],
                                                    u[s - 1, 0:(s - 1)].T)
                    m[s - 1] = np.dot(b[:, s - 1], b[:, s - 1])

                if k > (sz[1]-1):
                    k -= 1
                else:
                    # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                    p = np.dot(a[:, k:sz[1]].T, b[:, (k - (sz[1]-1)):k])
                    q = np.diag(m[(k - (sz[1]-1)):k])
                    result = np.linalg.lstsq(q.T, p.T)[0].T
                    u[k:sz[1], (k - (sz[1]-1)):k] = result

        # checking if the reduced matrix is right handed
        # if np.linalg.det(a) < 0:
        #     a[:, [1, 2]] = a[:, [2, 1]]
        return a
    else:
        raise Exception(
            'The input to the lll_algorithm is expected to be integral')
# -----------------------------------------------------------------------------------------------------------


def eq(m1, m2, tol):
    """
    Check if the two rotation matrices are the same
    """
    if m1.ndim == 2 and m2.ndim == 2:
        m = abs(m1 - m2)

        if np.amax(m) < tol:
            return True
        else:
            return False
    elif m1.ndim == 2:
        msz = np.shape(m2)[0]
        tmat1 = m1.reshape((1, 9))
        tmat2 = np.tile(tmat1, (msz, 1))
        tmat3 = tmat2.reshape(msz, 3, 3)

        m = abs(tmat3 - m2)
        max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
        if np.any(max1):
            return True
        else:
            return False

    elif m2.ndim == 2:
        msz = np.shape(m1)[0]
        tmat1 = m2.reshape(msz, (1, 9))
        tmat2 = np.tile(tmat1, (msz, 1))
        tmat3 = tmat2.reshape(msz, 3, 3)

        m = abs(m1 - tmat3)
        max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
        if np.any(max1):
            return True
        else:
            return False
    else:
        if np.shape(m1)[0] == np.shape(m2)[0]:
            m = abs(m1 - m2)
            max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
            return np.where(max1)
        else:
            raise Exception('Wrong Input Types')
# -----------------------------------------------------------------------------------------------------------


def message_display(CheckMatrix, Checknumber, Message, Precis):
    """
    This function displays a Message (passed as input) and gives and error
    in case the matrix passed to it is not integral.`
    """
    cond = int_man.int_check(CheckMatrix, Precis)
    print Checknumber, '.', Message, '-> ',
    txt = Col()
    if cond.all():
        txt.c_prnt('YES', 'yel')
    else:
        txt.c_prnt('<<<Error>>>', 'amber')
        raise Exception('Something wrong!!')
# -----------------------------------------------------------------------------------------------------------


def extgcd(x, y):
    """
    Return a tuple (u, v, d); they are the greatest common divisor d
    of two integers x and y and u, v such that d = x * u + y * v.
    """
    # Crandall & Pomerance "PRIME NUMBERS", Algorithm 2.1.4
    a, b, g, u, v, w = 1, 0, x, 0, 1, y
    while w:
        q, t = divmod(g, w)
        a, b, g, u, v, w = u, v, w, a-q*u, b-q*v, t
    if g >= 0:
        return a, b, g
    else:
        return -a, -b, -g
# -----------------------------------------------------------------------------------------------------------


def ehermite(a, b):
    """
    Elementary Hermite tranformation.
    For integers a and b, E = ehermite(a,b) returns
    an integer matrix with determinant 1 such that E * [a;b] = [g;0],
    where g is the gcd of a and b.
    E = ehermite(a,b)
    This function is in some ways analogous to GIVENS.
    John Gilbert, 415-812-4487, December 1993
    gilbert@parc.xerox.com
    Xerox Palo Alto Research Center

    Parameters
    ----------
    a, b: integers

    Returns
    -------
    E: numpy array 3x3
        integer matrix with determinant 1 such that E * [a;b] = [g;0],
        where g is the gcd of a and b.

    """
    [c, d, g] = extgcd(a, b)
    if g:
        E = np.array([[c, d], [-b/g, a/g]])
    else:
        E = np.array([[1, 0], [0, 1]])

    return E
# -----------------------------------------------------------------------------------------------------------


def left_matrix_division(X, Y):
    # solving the left matrix division X / Y
    # # ---------
    # # 1st alternative ---> Leaves unwanted decimals in some cases!
    tmp_solution = np.linalg.lstsq(Y.T, X.T)[0].T
    # # ---------
    # # 2nd alternative ---> Also Leaves unwanted decimals in some cases!
    # solution = np.dot(np.dot(X, Y.T), np.linalg.inv(np.dot(Y, Y.T)))
    # # ---------
    solution = (np.around(tmp_solution*1e10))/1e10
    return solution
# -----------------------------------------------------------------------------------------------------------


def smith_nf(matrix):
    """
    Smith normal form of an integer matrix.
    [U,S,V] = smith(A) returns integer matrices U, S, and V such that
    A = U*S*V',
    S is diagonal and nonnegative, S(i,i) divides S(i+1,i+1) for all i,
    det U =+-1, and det V =+-1.
    s = smith(A) just returns diag(S).
    Uses function ehermite.
    [U,S,V] = smith(A);

    This function is in some ways analogous to SVD.
    Originally implemented by: John Gilbert, 415-812-4487, December 1993
    gilbert@parc.xerox.com
    Xerox Palo Alto Research Center

    Parameters
    -----------
    matrix: numpy array

    Returns
    --------
    S: numpy array
        S is diagonal and nonnegative, S(i,i) divides S(i+1,i+1) for all i
    U: numpy array
        det(U) =+-1
    V: numpy array
        det(V) =+-1
    """

    A=np.copy(matrix)
    if (np.around(A) != A).any():
        raise Exception('This function requires integer input.')

    # This looks much like an SVD algorithm that first bidiagonalizes
    # A by Givens rotations and then chases zeros, except for
    # the construction of the 2 by 2 elementary transformation.

    m, n = A.shape

    S = A
    U = np.eye(m)
    V = np.eye(n)

    # Bidiagonalize S with elementary Hermite transforms.
    for j in range(min(m, n)):
        # Zero column j below the diagonal.
        for i in range(j+1, m):
            if S[i, j]:
                # Construct an elementary Hermite transformation E
                # to zero S(i,j) by combining rows i and j.
                E = ehermite(S[j, j], S[i, j])
                # Apply the transform to S and U.
                S[[j, i], :] = np.dot(E, S[[j, i], :])
                # U[:, [j, i]] = U[:, [j, i]] / E
                U[:, [j, i]] = left_matrix_division(U[:, [j, i]], E) # solving the left matrix division

        # % Zero row j after the superdiagonal.
        for i in range(j+2, n):
            if S[j, i]:
                # Construct an elementary Hermite transformation E
                # to zero S(j,i) by combining columns j+1 and i.
                E = ehermite(S[j, j+1], S[j, i])
                # Apply the transform to S and V.
                S[:, [j+1, i]] = np.dot(S[:, [j+1, i]], E.T)
                # V[:, [j+1, i]] = V[:, [j+1, i]] / E
                V[:, [j+1, i]] = left_matrix_division(V[:, [j+1, i]], E) # solving the left matrix division

    # Now S is upper bidiagonal.
    # Chase the superdiagonal nonzeros away.

    D = np.diag(S, 1)
    while any(D):
        b = min(np.where(D))[0]
        # Start chasing bulge at first nonzero superdiagonal element.
        # To guarantee reduction in S(b,b), first make S(b,b) positive
        # and make S(b,b+1) nonnegative and less than S(b,b).
        if S[b, b] < 0:
            S[b, :] = -S[b, :]
            U[:, b] = -U[:, b]

        q = np.floor(S[b, b+1] / S[b, b])
        E = np.array([[1, 0], [-q, 1]])
        S[:, [b, b+1]] = np.dot(S[:, [b, b+1]], E.T)
        # V[:, [b, b+1]] = V[:, [b, b+1]] / E
        V[:, [b, b+1]] = left_matrix_division(V[:, [b, b+1]], E) # solving the left matrix division

        if S[b, b+1]:
            # Zero the first nonzero superdiagonal element
            # using columns b and b+1, to start the bulge at S(b+1,b).
            E = ehermite(S[b, b], S[b, b+1])
            S[:, [b, b+1]] = np.dot(S[:, [b, b+1]], E.T)
            # V[:, [b, b+1]] = V[:, [b, b+1]] / E
            V[:, [b, b+1]] = left_matrix_division(V[:, [b, b+1]], E)

            for j in range(min(m, n)):
                if j+1 < m:
                    # Zero S(j+1,j) using rows j and j+1.
                    E = ehermite(S[j, j], S[j+1, j])
                    S[[j, j+1], :] = np.dot(E, S[[j, j+1], :])
                    # U[:, [j, j+1]] = U[:, [j, j+1]] / E
                    U[:, [j, j+1]] = left_matrix_division(U[:, [j, j+1]], E)
                if j+2 < n:
                    # Zero S(j,j+2) using columns j+1 and j+2.
                    E = ehermite(S[j, j+1], S[j, j+2])
                    S[:, [j+1, j+2]] = np.dot(S[:, [j+1, j+2]], E.T)
                    # V[:, [j+1, j+2]] = V[:, [j+1, j+2]] / E
                    V[:, [j+1, j+2]] = left_matrix_division(V[:, [j+1, j+2]], E)
        D = np.diag(S, 1)

    # Now S is diagonal. Make it nonnegative.

    for j in range(min(m, n)):
        if S[j, j] < 0:
            S[j, :] = -S[j, :]
            U[:, j] = -U[:, j]

    # Squeeze factors to lower right to enforce divisibility condition.

    for i in range(min(m, n)):
        for j in range(i+1, min(m, n)):
            # Replace S(i,i), S(j,j) by their gcd and lcm respectively.
            a = S[i, i]
            b = S[j, j]
            [c, d, g] = extgcd(a, b)
            E = np.array([[1, d], [-b/g, a*c/g]])
            F = np.array([[c, 1], [-b*d/g, a/g]])
            S[np.ix_([i, j], [i, j])] = np.dot(np.dot(E, S[:, [i, j]][[i, j], :]), F.T)
            # S[i, i] = tmp_arr[0, 0]
            # S[i, j] = tmp_arr[0, 1]
            # S[j, i] = tmp_arr[1, 0]
            # S[j, j] = tmp_arr[1, 1]
            U[:, [i, j]] = left_matrix_division(U[:, [i, j]], E)
            V[:, [i, j]] = left_matrix_division(V[:, [i, j]], F)

    U = np.around(U)
    V = np.around(V)
    return U, S, V
# -----------------------------------------------------------------------------------------------------------


def vrrotvec2mat(ax_ang):
    """
    Create a Rotation Matrix from Axis-Angle vector:

    Parameters
    ----------
    ``ax_ang``: numpy 5xn array
        The 3D rotation axis and angle (ax_ang) \v
        5 entries: \v
        First 3: axis \v
        4: angle \v
        5: 1 for proper and -1 for improper \v

    Returns
    -------
    mtx: nx3x3 numpy array
        3x3 rotation matrices

    See Also
    --------
    mat2quat, axang2quat, vrrotmat2vec
    """
    
    #file_dir = os.path.dirname(os.path.realpath(__file__))
    #path_dir2 = file_dir + '/../geometry/'
    #sys.path.append(path_dir2)
    
    if ax_ang.ndim == 1:
        if np.size(ax_ang) == 5:
            ax_ang = np.reshape(ax_ang, (5, 1))
            msz = 1
        elif np.size(ax_ang) == 4:
            ax_ang = np.reshape(np.hstack((ax_ang, np.array([1]))), (5, 1))
            msz = 1
        else:
            raise Exception('Wrong Input Type')
    elif ax_ang.ndim == 2:
        if np.shape(ax_ang)[0] == 5:
            msz = np.shape(ax_ang)[1]
        elif np.shape(ax_ang)[1] == 5:
            ax_ang = ax_ang.transpose()
            msz = np.shape(ax_ang)[1]
        else:
            raise Exception('Wrong Input Type')
    else:
        raise Exception('Wrong Input Type')

    direction = ax_ang[0:3, :]
    angle = ax_ang[3, :]

    d = np.array(direction, dtype=np.float64)
    d /= np.linalg.norm(d, axis=0)
    x = d[0, :]
    y = d[1, :]
    z = d[2, :]
    c = np.cos(angle)
    s = np.sin(angle)
    tc = 1 - c

    mt11 = tc*x*x + c
    mt12 = tc*x*y - s*z
    mt13 = tc*x*z + s*y

    mt21 = tc*x*y + s*z
    mt22 = tc*y*y + c
    mt23 = tc*y*z - s*x

    mt31 = tc*x*z - s*y
    mt32 = tc*y*z + s*x
    mt33 = tc*z*z + c

    mtx = np.column_stack((mt11, mt12, mt13, mt21, mt22, mt23, mt31, mt32, mt33))

    inds1 = np.where(ax_ang[4, :] == -1)
    mtx[inds1, :] = -mtx[inds1, :]

    if msz == 1:
        mtx = mtx.reshape(3, 3)
    else:
        mtx = mtx.reshape(msz, 3, 3)

    return mtx
# -----------------------------------------------------------------------------------------------------------


def vrrotmat2vec(mat1, rot_type='proper'):
    """
    Create an axis-angle np.array from Rotation Matrix:

    Parameters
    ----------
    mat1: nx3x3 numpy array
        The nx3x3 rotation matrices to convert
    rot_type: string ('proper' or 'improper')
        ``improper`` if there is a possibility of
        having improper matrices in the input,
        ``proper`` otherwise. \v
        Default: ``proper``

    Returns
    -------
    ``ax_ang``: numpy 5xn array
        The 3D rotation axis and angle (ax_ang) \v
        5 entries: \v
        First 3: axis \v
        4: angle \v
        5: 1 for proper and -1 for improper \v

    See Also
    --------
    mat2quat, axang2quat, vrrotvec2mat
    """
    mat = np.copy(mat1)
    if mat.ndim == 2:
        if np.shape(mat) == (3, 3):
            mat = np.copy(np.reshape(mat, (1, 3, 3)))
        else:
            raise Exception('Wrong Input Type')
    elif mat.ndim == 3:
        if np.shape(mat)[1:] != (3, 3):
            raise Exception('Wrong Input Type')
    else:
        raise Exception('Wrong Input Type')

    msz = np.shape(mat)[0]
    ax_ang = np.zeros((5, msz))

    epsilon = 1e-12
    if rot_type == 'proper':
        ax_ang[4, :] = np.ones(np.shape(ax_ang[4, :]))
    elif rot_type == 'improper':
        for i in range(msz):
            det1 = np.linalg.det(mat[i, :, :])
            if abs(det1 - 1) < epsilon:
                ax_ang[4, i] = 1
            elif abs(det1 + 1) < epsilon:
                ax_ang[4, i] = -1
                mat[i, :, :] = -mat[i, :, :]
            else:
                raise Exception('Matrix is not a rotation: |det| != 1')
    else:
        raise Exception('Wrong Input parameter for rot_type')



    mtrc = mat[:, 0, 0] + mat[:, 1, 1] + mat[:, 2, 2]


    ind1 = np.where(abs(mtrc - 3) <= epsilon)[0]
    ind1_sz = np.size(ind1)
    if np.size(ind1) > 0:
        ax_ang[:4, ind1] = np.tile(np.array([0, 1, 0, 0]), (ind1_sz, 1)).transpose()


    ind2 = np.where(abs(mtrc + 1) <= epsilon)[0]
    ind2_sz = np.size(ind2)
    if ind2_sz > 0:
        # phi = pi
        # This singularity requires elaborate sign ambiguity resolution

        # Compute axis of rotation, make sure all elements >= 0
        # real signs are obtained by flipping algorithm below
        diag_elems = np.concatenate((mat[ind2, 0, 0].reshape(ind2_sz, 1),
                                     mat[ind2, 1, 1].reshape(ind2_sz, 1),
                                     mat[ind2, 2, 2].reshape(ind2_sz, 1)), axis=1)
        axis = np.sqrt(np.maximum((diag_elems + 1)/2, np.zeros((ind2_sz, 3))))
        # axis elements that are <= epsilon are set to zero
        axis = axis*((axis > epsilon).astype(int))

        # Flipping
        #
        # The algorithm uses the elements above diagonal to determine the signs
        # of rotation axis coordinate in the singular case Phi = pi.
        # All valid combinations of 0, positive and negative values lead to
        # 3 different cases:
        # If (Sum(signs)) >= 0 ... leave all coordinates positive
        # If (Sum(signs)) == -1 and all values are non-zero
        #   ... flip the coordinate that is missing in the term that has + sign,
        #       e.g. if 2AyAz is positive, flip x
        # If (Sum(signs)) == -1 and 2 values are zero
        #   ... flip the coord next to the one with non-zero value
        #   ... ambiguous, we have chosen shift right

        # construct vector [M23 M13 M12] ~ [2AyAz 2AxAz 2AxAy]
        # (in the order to facilitate flipping):    ^
        #                                  [no_x  no_y  no_z ]

        m_upper = np.concatenate((mat[ind2, 1, 2].reshape(ind2_sz, 1),
                                  mat[ind2, 0, 2].reshape(ind2_sz, 1),
                                  mat[ind2, 0, 1].reshape(ind2_sz, 1)), axis=1)

        # elements with || smaller than epsilon are considered to be zero
        signs = np.sign(m_upper)*((abs(m_upper) > epsilon).astype(int))

        sum_signs = np.sum(signs, axis=1)
        t1 = np.zeros(ind2_sz,)
        tind1 = np.where(sum_signs >= 0)[0]
        t1[tind1] = np.ones(np.shape(tind1))

        tind2 = np.where(np.all(np.vstack(((np.any(signs == 0, axis=1) == False), t1 == 0)), axis=0))[0]
        t1[tind2] = 2*np.ones(np.shape(tind2))

        tind3 = np.where(t1 == 0)[0]
        flip = np.zeros((ind2_sz, 3))
        flip[tind1, :] = np.ones((np.shape(tind1)[0], 3))
        flip[tind2, :] = np.copy(-signs[tind2, :])

        t2 = np.copy(signs[tind3, :])

        shifted = np.column_stack((t2[:, 2], t2[:, 0], t2[:, 1]))
        flip[tind3, :] = np.copy(shifted + (shifted == 0).astype(int))

        axis = axis*flip
        ax_ang[:4, ind2] = np.vstack((axis.transpose(), np.pi*(np.ones((1, ind2_sz)))))

    ind3 = np.where(np.all(np.vstack((abs(mtrc + 1) > epsilon, abs(mtrc - 3) > epsilon)), axis=0))[0]
    ind3_sz = np.size(ind3)
    if ind3_sz > 0:
        phi = np.arccos((mtrc[ind3]-1)/2)
        den = 2*np.sin(phi)
        a1 = (mat[ind3, 2, 1]-mat[ind3, 1, 2])/den
        a2 = (mat[ind3, 0, 2]-mat[ind3, 2, 0])/den
        a3 = (mat[ind3, 1, 0]-mat[ind3, 0, 1])/den
        axis = np.column_stack((a1, a2, a3))
        ax_ang[:4, ind3] = np.vstack((axis.transpose(), phi.transpose()))

    return ax_ang
# -----------------------------------------------------------------------------------------------------------


def quat2mat(q):
    """
    Convert Quaternion Arrays to Rotation Matrix

        Parameters
    ----------
    q: numpy array (5 x 1)
        quaternion

    Returns
    ----------
    g: numpy array (3 x 3)
        rotation matrix

    See Also
    --------
    mat2quat, axang2quat
    """
    import quaternion as quat
    sz = quat.get_size(q)
    q0 = quat.getq0(q)
    q1 = quat.getq1(q)
    q2 = quat.getq2(q)
    q3 = quat.getq3(q)
    qt = quat.get_type(q)

    g = np.zeros((sz, 3, 3))
    g[:, 0, 0] = np.square(q0) + np.square(q1) - np.square(q2) - np.square(q3)
    g[:, 0, 1] = 2*(q1*q2 - q0*q3)
    g[:, 0, 2] = 2*(q3*q1 + q0*q2)
    g[:, 1, 0] = 2*(q1*q2 + q0*q3)
    g[:, 1, 1] = np.square(q0) - np.square(q1) + np.square(q2) - np.square(q3)
    g[:, 1, 2] = 2*(q2*q3 - q0*q1)
    g[:, 2, 0] = 2*(q3*q1 - q0*q2)
    g[:, 2, 1] = 2*(q2*q3 + q0*q1)
    g[:, 2, 2] = np.square(q0) - np.square(q1) - np.square(q2) + np.square(q3)

    if sz == 1:
        g = g.reshape((3, 3))
        if qt == -1:
            g = -g
    else:
        inds1 = np.where(qt == -1)
        g[inds1, :, :] = -g[inds1, :, :]

    return g
# -----------------------------------------------------------------------------------------------------------


def mat2quat(mat, rot_type='proper'):
    """
    Convert Rotation Matrices to Quaternions

    Parameters
    ----------
    mat: numpy array or a list of (3 x 3)
        rotation matrix

    rot_type: string ('proper' or 'improper')
        ``improper`` if there is a possibility of
        having improper matrices in the input,
        ``proper`` otherwise. \v
        Default: ``proper``

    Returns
    ----------
    quaternion_rep: numpy array (5 x 1)

    See Also
    --------
    quat2mat, axang2quat
    """
    import quaternion as quat
    ax_ang = vrrotmat2vec(mat, rot_type)
    q0 = np.cos(ax_ang[3, :]/2)
    q1 = ax_ang[0, :]*np.sin(ax_ang[3, :]/2)
    q2 = ax_ang[1, :]*np.sin(ax_ang[3, :]/2)
    q3 = ax_ang[2, :]*np.sin(ax_ang[3, :]/2)
    qtype = ax_ang[4, :]

    return quat.Quaternion(q0, q1, q2, q3, qtype)
# -----------------------------------------------------------------------------------------------------------


def axang2quat(ax_ang):
    """
    Create a quaternion corresponding to the rotation specified by an axis and an angle

    Parameters
    ----------
    ax_ang: numpy array or a list of (4 x 1)

    Returns
    ----------
    quaternion_rep: numpy array (5 x 1)
    """
    import quaternion as quat

    if ax_ang.ndim == 1:
        if np.size(ax_ang) == 5:
            ax_ang = np.reshape(ax_ang, (5, 1))
            msz = 1
        else:
            raise Exception('Wrong Input Type')
    elif ax_ang.ndim == 2:
        if np.shape(ax_ang)[0] == 5:
            msz = np.shape(ax_ang)[1]
        elif np.shape(ax_ang)[1] == 5:
            ax_ang = ax_ang.transpose()
            msz = np.shape(ax_ang)[1]
        else:
            raise Exception('Wrong Input Type')
    else:
        raise Exception('Wrong Input Type')

    direction = ax_ang[0:3, :]
    angle = ax_ang[3, :]

    d = np.array(direction, dtype=np.float64)
    d /= np.linalg.norm(d, axis=0)
    x = d[0, :]
    y = d[1, :]
    z = d[2, :]
    q0 = np.cos(angle/2)
    s = np.sin(angle/2)

    q1 = x*s
    q2 = y*s
    q3 = z*s

    qtype = ax_ang[4, :]
    return quat.Quaternion(q0, q1, q2, q3, qtype)
# -----------------------------------------------------------------------------------------------------------

# mats = np.zeros((20, 3, 3))
# ct1 = 0
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 1, 0, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 1, 1, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([-2, 1, 1, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, 2, 3, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([-1, 2, 3, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, -2, 3, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, 2, -3, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 2, -3, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0.34, 0.68, -0.94, -np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 1, 0, np.pi]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 1, 0, 0]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 1, 1, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([-2, 1, 1, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, 2, 3, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([-1, 2, 3, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, -2, 3, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, 2, -3, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 2, -3, np.pi/2]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([1, 0, 0, 0]))
# ct1 += 1
# mats[ct1, :, :] = vrrotvec2mat(np.array([0, 0, 1, 0]))
# ct1 += 1
#
#
# axang = vrrotmat2vec(mats)
#
# Mats = vrrotvec2mat(axang)
#
# for i in range(np.shape(Mats)[0]):
#     print vrrotmat2vec(np.dot(Mats[i, :, :],(mats[i, :, :].transpose())))
# -----------------------------------------------------------------------------------------------------------


def unique_rows_tol(data, tol=1e-12, return_index=False, return_inverse=False):
    """
    This function returns the unique rows of the input matrix within that are within the
    specified tolerance.

    Parameters
    ----------
    data: numpy array (m x n)
    tol: double
        tolerance of comparison for each rows
        Default: 1e-12
    return_index: Boolean
        flag to return the index of unique rows based on the indices of the output
    return_inverse: Boolean
        flag to return the index of unique rows based on the indices of the input

    Returns
    ----------
    unique_rows: numpy array (m' x n)
    ia: numpy array, integer (m' x 1)
        unique rows based on the indices of the output
    ic: numpy array, integer (m x 1)
        unique rows based on the indices of the input

    See Also
    --------
    unique
    """
    prec = -np.fix(np.log10(tol))
    d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
    ### fix rounds off towards zero; issues with the case of 0.9999999998 and 1.0

    ### rint solves the issue, needs extensive testing
    # prec = -np.rint(np.log10(tol))
    # d_r = np.rint(data * 10 ** prec) / 10 ** prec + 0.0

    b = np.ascontiguousarray(d_r).view(np.dtype((np.void, d_r.dtype.itemsize * d_r.shape[1])))
    _, ia = np.unique(b, return_index=True)
    _, ic = np.unique(b, return_inverse=True)

    ret_arr = data[ia, :]
    if not return_index and not return_inverse:
        return ret_arr
    else:
        if return_index and return_inverse:
            return ret_arr, ia, ic
        elif return_index:
            return ret_arr, ia
        elif return_inverse:
            return ret_arr, ic

    # if not return_index and not return_inverse:
    #     return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1])
    # else:
    #     if return_index and return_inverse:
    #         return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia, ic
    #     elif return_index:
    #         return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia
    #     elif return_inverse:
    #         return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ic
# -----------------------------------------------------------------------------------------------------------


def test_unique_rows():
    prec = 1.e-5
    mat = np.array([[-1e-6, 1, 1], [1e-7, 1, 1], [0, 1, 1], [1, 1, 1]])
    c, ia, ic = unique_rows_tol(mat, prec, True, True)
    print unique_rows_tol(mat, prec, True, True)
# -----------------------------------------------------------------------------------------------------------


def test_lll_reduction():
    Mat = np.zeros(5, dtype=[('Mat', '(3,2)float64'),('Matrix', '(3,3)float64')])
    Mat['Mat'][0] = np.array(([1, 2, 0], [-4, 5, 6])).T
    Mat['Mat'][1] = np.array(([-1, 2, 3], [-4, 0, 6])).T
    Mat['Mat'][2] = np.array(([1, 2, 0], [-4, 5, 6])).T
    Mat['Mat'][3] = np.array(([1, 2, 3], [-4, 5, 6])).T
    Mat['Mat'][4] = np.array(([1, 2, 25], [-4, 5, 6])).T
    Mat['Matrix'][0] = np.array(([1, 2, 0], [-4, 5, 6], [7, 8, 15]))
    Mat['Matrix'][1] = np.array(([-1, 2, 3], [-4, 0, 6], [3, 0, 0]))
    Mat['Matrix'][2] = np.array(([1, 2, 0], [-4, 5, 6], [3, 0, -1]))
    Mat['Matrix'][3] = np.array(([1, 2, 3], [-4, 5, 6], [2, -2, 1]))
    Mat['Matrix'][4] = np.array(([1, 2, 25], [-4, 5, 6], [3, -2, 0]))
    for j in Mat.dtype.names:
        for i in range(Mat.shape[0]):
            # a, H = lll_reduction_3by2(Mat['Matrix'][i])
            b = lll_reduction(Mat[j][i])
            # print Mat['Matrix'][i], '\n reduced: \n', H, '\n-------\n'
            print '\n______________________________________________\n'
            print Mat[j][i], '\n reduced: \n', b, '\n-------\n'
# -----------------------------------------------------------------------------------------------------------
