# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
from fractions import gcd
from fractions import Fraction
# -----------------------------------------------------------------------------------------------------------

def gcd_array(input, order='all'):
    """
    The function computes the GCD of an array of numbers.

    Parameters
    ----------
    input : numpy array or list of intgers
        Input n-D array of integers (most suitable for 1D and 2D arrays)

    order : {'rows', 'columns', 'col', 'all'}, optional

    Returns
    -------
    Agcd: numpy array
        An array of greatest common divisors of the input

    Notes
    -------
    * If order = **all**, the input array is flattened and the GCD is calculated
    * If order = **rows**, GCD of elements in each row is calculated
    * If order = **columns** or **cols**, GCD of elements in each column is calculated

    See Also
    --------
    gcd: from fractions module for computing gcd of two integers

    """

    # Vectorizing the function gcd
    gcd_vec = np.vectorize(gcd)
    tmp = 0

    input = np.array(input)

    if np.ndim(input) == 1:
        input = np.reshape(input, (1, input.shape[0]))

    # Only integer values are allowed
    # if input.dtype.name != 'int64':
    if not np.issubdtype(input.dtype, np.int):
        raise Exception("Inputs must be real integers.")

    err_msg = "Not a valid input. Please choose either \"rows\" " + \
              "or \"columns\" keys for this function."

    order_options = ('rows', 'columns', 'col', 'all')
    try:
        Keys = (order_options.index(order))
    except:
        raise Exception(err_msg)

    if (Keys == 3):
        Sz = input.shape
        input = np.reshape(input, (1, Sz[0]*Sz[1]))
    if (Keys == 1) or (Keys == 2):
        input = input.T
        Agcd = gcd_array(input, 'rows')
        # handling the case of asking a column
        # vector with the 'row' key by mistake.
        tmp = 1

    Sz = input.shape
    if Sz[1] == 1:
        input.shape = (1, Sz[0])

    Agcd = gcd_vec(input[::, 0], input[::, 1])
    for i in range(Sz[1]-2):
        Agcd = gcd_vec(Agcd, input[::, i+2])

    if tmp != 1:
        Agcd.shape = (len(Agcd), 1)

    return np.absolute(Agcd)
# -----------------------------------------------------------------------------------------------------------


def lcm_vec(a, b):
    """
    The function computes the LCM of two 1D array of integers of length
    and retruns a 1D array of lcm values

    Parameters
    ----------
    a, b : numpy array
        Input 1D arrays of integers

    Returns
    -------
    lcm_vector: numpy array
        Output 1D array of integers

    See Also
    --------
    lcm_arry

    """
    gcd_vec = np.vectorize(gcd)
    lcm_vector = a * (b / gcd_vec(a, b))
    return lcm_vector
# -----------------------------------------------------------------------------------------------------------


def lcm_array(input, order='all'):
    """
    The function computes the LCM of an array of numbers.

    Parameters
    ----------
    input : numpy array or list of intgers
        Input n-D array of integers (most suitable for 1D and 2D arrays)

    order : {'rows', 'columns', 'col', 'all'}, optional

    Returns
    -------
    Alcm: numpy array
        An array of least common multiples of the input

    Notes
    -------
    * If order = **all**, the input array is flattened and the LCM is calculated
    * If order = **rows**, LCM of elements in each row is calculated
    * If order = **columns** or **cols**, LCM of elements in each column is calculated

    See Also
    --------
    gcd_array

    """
    input = np.array(input)
    tmp = 0

    if np.ndim(input) == 1:
        input = np.reshape(input, (1, input.shape[0]))

    # Only integer values are allowed
    # if input.dtype.name != 'int64':
    if not np.issubdtype(input.dtype, np.int):
        raise Exception("Inputs must be real integers.")

    err_msg = "Not a valid input. Please choose either \"rows\" " + \
              "or \"columns\" keys for this function."

    order_options = ('rows', 'columns', 'col', 'all')
    try:
        Keys = (order_options.index(order))
    except:
        raise Exception(err_msg)

    if (Keys == 3):
        Sz = input.shape
        input = np.reshape(input, (1, Sz[0]*Sz[1]))
    if (Keys == 1) or (Keys == 2):
        input = input.T
        Alcm = lcm_array(input, 'rows')
        # handling the case of asking a column vector
        # with the 'row' key by mistake.
        tmp = 1

    Sz = input.shape
    if Sz[1] == 1:
        input.shape = (1, Sz[0])

    Alcm = lcm_vec(input[::, 0], input[::, 1])
    for i in range(Sz[1]-2):
        Alcm = lcm_vec(Alcm, input[::, i+2])

    if tmp != 1:
        Alcm.shape = (len(Alcm), 1)

    return np.absolute(Alcm)
# -----------------------------------------------------------------------------------------------------------


def int_check(input, precis=6):
    """
    Checks whether the input variable (arrays) is an interger or not.
    A precision value is specified and the integer check is performed
    up to that decimal point.

    Parameters
    ----------
    input : numpy array or list
        Input n-D array of floats.

    precis : Integer
        Default = 6.
        A value that specifies the precision to which the number is an
        integer. **precis = 6** implies a precision of :math:`10^{-6}`.

    Returns
    -------
    cond: Boolean
        **True** if the element is an integer to a certain precision,
        **False** otherwise
    """

    var = np.array(input)
    tval = 10 ** -precis
    t1 = abs(var)
    cond = (abs(t1 - np.around(t1)) < tval)
    return cond
# -----------------------------------------------------------------------------------------------------------


def rat(input, tol=1e-06):
    """
    The function returns a rational (p/q) approximation of a given
    floating point array to a given precision

    Parameters
    ----------
    input : numpy array or list of real numbers

    tol : floating point tolerance value
        Default = 1e-06

    Returns
    -------
    N, D: Integer numpy arrays
        N and D contain the numerators (p) and denominators (q) of the
        rational approximations

    Notes:
    --------
    """
    input1 = np.array(input)
    if np.ndim(input1) == 1:
        input1 = np.reshape(input1, (1, input1.shape[0]))

    ## Why is this case necessary?
    if input1.ndim == 0:
        input1 = np.reshape(input1, (1, 1))

    Sz = input1.shape
    N = np.zeros((Sz[0], Sz[1]), dtype='int64')
    D = np.zeros((Sz[0], Sz[1]), dtype='int64')
    nDec = int(1/tol)
    for i in range(Sz[0]):
        for j in range(Sz[1]):
            N[i, j] = (Fraction.from_float(input1[i, j]).
                       limit_denominator(nDec).numerator)
            D[i, j] = (Fraction.from_float(input1[i, j]).
                       limit_denominator(nDec).denominator)
    return N, D
# -----------------------------------------------------------------------------------------------------------


def int_finder(input_v, tol=1e-6, order='all', tol1=1e-6):
    """
    The function computes the scaling factor required to multiply the
    given input array to obtain an integer array. The integer array is
    returned.

    Parameters
    ----------
    input1 : numpy array or list of real numbers

    tol : floating point tolerance value
        Default = 1e-06

    order : {'rows', 'columns', 'col', 'all'}
        Defualt = 'all'

    tol1:

    Returns
    -------
    output: numpy float array
    An array of integers obtained by scaling input

    See Also
    --------
    gcd_array

    Notes
    --------
    * If order = **all**, the input array is flattened and then scaled
    * If order = **rows**, elements in each row are scaled
    * If order = **columns** or **cols**, elements in each column are scaled
    """

    input1 = np.array(input_v)
    Sz = input1.shape
    if np.ndim(input1) == 1:
        input1 = np.reshape(input1, (1, input1.shape[0]))

    if int_check(input1, 15).all():
        input1 = np.around(input1)
        # Divide by LCM (rows, cols, all) <--- To Do
        tmult = gcd_array(input1.astype(dtype='int64'), order)
        if (order == 'all'):
            input1 = input1 / tmult
        elif (order == 'rows'):
            tmult = np.tile(tmult, (np.shape(input1[1])))
            input1 = input1 / tmult
        elif (order == 'col' or order == 'cols' or order == 'columns'):
            tmult = np.tile(tmult, (np.shape(input1[0])[0], 1))
            input1 = input1 / tmult
        output_v = input1
        if len(Sz) == 1:
            output_v = np.reshape(output_v, (np.size(output_v), ))
        return output_v
    else:
        #   By default it flattens the array (if nargin < 3)
        if order.lower() == 'all':
            if len(Sz) != 1:
                input1.shape = (1, Sz[0]*Sz[1])
        else:
            Switch = 0
            err_msg = "Not a valid input. For the third argument please"+ \
                      " choose either \"rows\" or \"columns\" keys for this function."
            order_options = ('rows', 'columns', 'col')
            try:
                Keys = (order_options.index(order.lower()))
            except:
                raise Exception(err_msg)

            if (Keys == 1) or (Keys == 2):
                if input1.shape[0] != 1:
                    # Handling the case of asking a row vector
                    # with the 'column' key by mistake.
                    input1 = input1.T
                    Switch = 1
            # Handling the case of asking a column
            # vector with the 'row' key by mistake.
            if (Keys == 0) and (input1.shape[1] == 1):
                input1 = input1.T
                Switch = 1

        if (abs(input1) < tol).all():
            excep1 = 'All the input components cannot' \
                     + 'be smaller than tolerance.'
            raise Exception(excep1)

        tmp = np.array((abs(input1) > tol1))
        Vec = 2 * abs(input1[::]).max() * np.ones(
            (input1.shape[0], input1.shape[1]))
        Vec[tmp] = input1[tmp]
        MIN = abs(Vec).min(axis=1)
        # Transposing a row to a column
        MIN.shape = (len(MIN), 1)
        input1 = input1 / np.tile(MIN, (1, input1.shape[1]))
        N, D = rat(input1, tol)
        N[~tmp] = 0 # <---- added
        D[~tmp] = 1 # <---- added
        lcm_rows = lcm_array(D, 'rows')
        lcm_mat = np.tile(lcm_rows, (1, input1.shape[1]))
        Rounded = (N * lcm_mat) / D
        output_v = Rounded

        # --------------------------
        if order.lower() == 'all':
            if len(Sz) != 1:
                output_v.shape = (Sz[0], Sz[1])
        else:
            if (Keys) == 1 or (Keys) == 2:
                output_v = output_v.T
            if Keys == 0 and Switch == 1:
                output_v = output_v.T

        if len(Sz) == 1:
            output_v = np.reshape(output_v, (np.size(output_v), ))

        return output_v
# -----------------------------------------------------------------------------------------------------------


def int_mult(input, tol=1e-06):
    """
    The function computes the scaling factor required to multiply the
    given input array to obtain an integer array. The integer array is
    returned.

    Parameters
    ----------
    input : numpy array or list of real numbers

    tol : floating point tolerance value
        Default = 1e-06

    Returns
    -------
    N: numpy float array
        An array of integers obtained by scaling input

    Int_Mat: numpy float array
        An array of integers obtained by scaling input

    See Also
    --------
    int_finder

    Notes
    --------
    **Change this function to accept rows and columns as input**
    """

    T = np.array(input)
    IntMat = int_finder(T, tol)
    t_ind = np.where(abs(IntMat) == abs(IntMat).max())
    t_ind_x = t_ind[0][0]
    t_ind_y = t_ind[1][0]
    N = IntMat[t_ind_x, t_ind_y] /input[t_ind_x, t_ind_y]
    return N, IntMat
# -----------------------------------------------------------------------------------------------------------