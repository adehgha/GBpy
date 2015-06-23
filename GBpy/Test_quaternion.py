# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import quaternion as quat


#############################################################
def pick_rand2():
    """
    Pick two random numbers such that x1**2 + x2**2 < 1
    
    Parameters
    -----------
    None

    Returns
    -----------
    Returns two random numbers
    x1, x2: In the range (-1, 1) and x1**2 + x2**2 < 1
    """
    x1 = 1; x2 = 1;
    while (x1**2 + x2**2 >= 1) :
        x1 = np.random.rand()
        if (np.random.rand() < 0.5):
            x1 = -x1
        x2 = np.random.rand()
        if (np.random.rand() < 0.5):
            x2 = -x2
    return x1, x2

def rand_quat(n):
    """
    Create a random quaternion

    Parameters
    -----------
    n: Integer
    number of random quaternions to be returned

    Returns
    -----------
    returns *n* random quaternions

    """
    q0 = np.zeros((n, )); q1 = np.zeros((n, ));
    q2 = np.zeros((n, )); q3 = np.zeros((n, ));
    for ct1 in range(n):
        x1, x2 = pick_rand2(); x3, x4 = pick_rand2();
        q0[ct1] = x1; q1[ct1] = x2;
        q2[ct1] = x3*(np.sqrt((1 - x1**2 - x2**2)/(x3**2 + x4**2)))
        q3[ct1] = x4*(np.sqrt((1 - x1**2 - x2**2)/(x3**2 + x4**2)))

    print q0**2 + q1**2 + q2**2 + q3**2
    return quat.Quaternion(q0, q1, q2, q3)

# test_quats = rand_quat(10)

def test_quaternion(type1, case):
    """
    Tests the quaternion class

    Parameters
    -----------
    type1: allowed values are 'Init' and 'Slice'
    * string

    case: allowed values are 0 and 5
    * integer

    Returns
    --------
    None

    Notes
    ------
    * If type1 == 'Init' and case == 0
    An empty quaternion array is created.
    * If type1 == 'Init' and case == 5
    A quaternion array of size (5 x 10) is created.
    * If type1 == 'Slice' and case == 5
    A quaternion array of size (5 x 10) is created and a new quaternion array of size (5 x 3) is created
    by slicing this array.

    See Also
    ---------
    * test_indiv_mehtods
    * rand_quat
    """
    if type1 == 'Init':
        if case == 0:
            q1 = quat.Quaternion()
            test_indiv_methods(q1)
        if case == 5:
            ctn = 10
            q1 = rand_quat(ctn)
            test_indiv_methods(q1)

    if type1 == 'Slice':
        if case == 5:
            ctn = 10
            q1 = rand_quat(ctn)
            q2 = q1[:, 0:3]

            test_indiv_methods(q1)
            test_indiv_methods(q2)

    return

def test_indiv_methods(q):
    """
    Prints the output for various methods in the Quaternion class for an input quaternion array

    Parameters
    -----------
    q: Input quaternion array
    * A quaternion array of size (5 x n)

    Returns
    --------
    None
    """
    print 'in human readable form\n', quat.display(q), '\n'
    # quat.display(q)
    
    q0 = quat.getq0(q)
    q1 = quat.getq1(q)
    q2 = quat.getq2(q)
    q3 = quat.getq3(q)

    print 'q0 component/s', q0, '\n'
    print 'q1 component/s', q1, '\n'
    print 'q2 component/s', q2, '\n'
    print 'q3 component/s', q3, '\n'

    print 'type ', quat.get_type(q), '\n'
    print 'size is ', quat.get_size(q), '\n'

    print 'antipodal ', quat.antipodal(q), '\n'
    print 'inverse', quat.inverse(q), '\n'

    st1 = quat.quat2mat(q)
    print 'quat2mat ', st1, '\n'
    print 'mat2quat, using matrices generated in test above', quat.mat2quat(st1), '\n'

    return


# test_quaternion('Init', 0)
# type1 = 'Init'; case = 5
# type1 = 'Slice'; case = 5
# test_quaternion(type1, case)

# ### Testing mtimes
# q1 = rand_quat(5); q2 = rand_quat(5);
# quat_mult = quat.mtimes(q1, q2);
# print quat.display(quat_mult);

# ### Testing mtimes
# q1 = rand_quat(5); q2 = rand_quat(1);
# quat_mult = quat.mtimes(q1, q2);
# print quat.display(quat_mult);

# ### Testing mtimes
# q1 = rand_quat(1); q2 = rand_quat(5);
# quat_mult = quat.mtimes(q1, q2);
# print quat.display(quat_mult);


### Testing eq
q1 = rand_quat(5); q2 = rand_quat(1);
print quat.eq(q1, q2);

### Testing eq
q1 = rand_quat(1); q2 = rand_quat(5);
q2[:, 2] = q1;
print quat.eq(q1, q2);


### Testing eq
q1 = rand_quat(5); q2 = rand_quat(5);
q2[:, 2] = q1[:, 2]; q2[:, 4] = q1[:, 4]; q2[4,4] = -1;
print quat.eq(q1, q2);






### Test eq




