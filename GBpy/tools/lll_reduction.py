# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np
import integer_manipulations as int_man


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

    Args:
        delta (float): Reduction parameter. Default of 0.75 is usually
            fine.
    Returns:
        Reduced lattice.
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


def test_lll_reduction():
    import numpy as np

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

# test_lll_reduction()
