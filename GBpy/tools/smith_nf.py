# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np

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
    """
    import numpy as np
    from smith_nf import extgcd

    [c, d, g] = extgcd(a, b)
    if g:
        E = np.array([[c, d], [-b/g, a/g]])
    else:
        E = np.array([[1, 0], [0, 1]])

    return E
# -----------------------------------------------------------------------------------------------------------


def left_matrix_division(X, Y):
    import numpy as np
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



