# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import tools as tools
import geometry_tools as gmt


######################################################################
class Quaternion(np.ndarray):
    """
    The Quaternion class is defined to represent the quaternion parameterization of rotation and their operations.
    The object of this class is stored in form of a quaternion array, which inherits all the properties of
    numpy.ndarray. The dimension of this array is (5 x n). Each quaterinion in the quaterinion array is referred
    as q(i). Each quaterinion has 5 attributes saved, namely, q0, q1, q2, q3 and type. Where qn have the usual
    meaning and type implies proper or improper type of rotations. Type has allowed values +1 and -1,
    for proper and improper rotations, respectively.

    The class has the following methods 
    getq0, 
    getq1, 
    getq2, 
    getq3, 
    get_type, 
    get_size, 
    display, 
    antipodal, 
    inverse,
    mtimes, 
    eq, 
    quat2mat
    mat2quat.
    
    Please read the documentation included with each method to learn more.
    ...

    Attributes
    -----------
    Quaternion: quaternion array, has all the properties of a numpy.ndarray
    * array of dimension (5 x n)
    * data type is Quaternion.

    Methods
    --------
    getq0()
    * Returns the q0 component for each quaternion present in the input quaternion array.

    getq1()
    * Returns the q1 component for each quaternion present in the input quaternion array.

    getq2()
    * Returns the q2 component for each quaternion present in the input quaternion array.

    getq3()
    * Returns the q3 component for each quaternion present in the input quaternion array.

    get_type()
    * Returns the rotation type of each quaternion present in the input quaternion array.

    get_size()
    * Returns the size of input quaternion array.

    display()
    * Returns a string which displays the input quaternion array in human readable format.

    antipodal()
    * Returns the antipodal (or equivalent) quaternions such that q0 component is positive.

    inverse()
    * Returns the inverse quaternions for a given input quaternion array.

    mtimes()
    * Calculates the quaternion multiplication of two input quaternions.

    eq()
    * Checks whether the two input quaternions are equal or not.

    quat2mat()
    * Converts the input quaternion array to a rotation matrix array.

    mat2quat()
    * Converts the input rotation matrix array or list, to a quaternion array.

    Notes
    ------
    * An object of this class can be created using 0, 1, 4 or 5, input arguments.
    * If the rotation type is not specified explicitly then it is always set to the default value of 1, i.e.
    proper rotation type.
    * If no argument is specified, an empty quaternion array of proper rotation type is created. Note that
    this is not an actual rotation. The sum of the squares of components of a rotation quaternion must
    be unity.
    * In case at least 1 argument is specified the object may be created by one of the following allowed cases.

        Case 1:
        * narg == 1; type(arg) == Quaternion
        A copy the argument is created as an object.

        Case2:
        * narg == 1; type(arg) == int
        A (5 x arg) quaternion array is created and populated with NaN values for components.

        Case3:
        * narg ==1; type(arg) == np.ndarray; dim(arg) == 1
        A 5 x 1 quaternion array is created from input argument array only if it has 4 or 5 elements.

        Case4:
        * narg == 1; type(arg) == np.ndarray; dim(arg) == 2
        If the shape of argument array is (4 x n) or (5 x n) corresponding (5 x n) quaternion array is created.
        If the shape of argument array is (n x 4) or (n x 5) then it is transposed and corresponding 5 x n 
	quaternion array is created.
        Note that the allowed shapes for the input argument in this case are (4/5 x n) or (n x 4/5)

        Case5:
        * narg == 4 or 5; type(arg) == float or integer or tuple with numeric values or np.ndarray
        If 4 or 5 arguments are given then a quaternion array may only be created if all of the arguments are numeric.
        The total number of elements for each input argument should be the same.
        If an argument is a multi dimensional numpy array then it should be symmetric in shape and total number of
        inner elements in this array should be the same as the number of elements for rest of the arguments. In such
	cases the np.ndarray is treated as a one dimensional array, i.e. a (2 x2) array is used as a (1 x 4) array would.
        Arguments may be of different type, i.e. a tuple of numerals and a multi dimension np.ndarray may be used
        simultaneously for a quaternion array creation.
        The n (4 or 5) arguments should contain the respective qn components for the quaternions that are to be stored
        in the quaternion array.

    * The above 5 cases are the only allowed ways to create an object of this class. An error is raised for any other
     case.

     See Also
     ---------
     * numpy.ndarray
     * geometry_tools.isnumeric
    """
    def __new__(cls, *args):
        nargs = len(args)

        ##############################################################
        if nargs == 0:
            quatarray = np.empty((5, ))
            quatarray[4] = 1
            obj = quatarray.view(cls)

        ##############################################################
        elif nargs == 1:

            if type(args[0]) is Quaternion:
                quatarray = args[0]

                if np.shape(args[0])[0] == 5:
                    obj = quatarray
                else:
                    raise Exception('Wrong Quaternion Slice')

            elif type(args[0]) is int:
                quatarray = np.empty((5, args[0]))
                quatarray[:4, :] = np.NaN
                quatarray[4, :] = 1
                obj = quatarray.view(cls)

            elif type(args[0]) is np.ndarray:
                q_shape = np.shape(args[0])

                if args[0].ndim == 1:

                    if q_shape[0] == 4:
                        quat1 = args[0].reshape((4, ))
                        v_type = np.ones((1, ))
                        quatarray = np.concatenate((quat1, v_type))

                    elif q_shape[0] == 5:
                        quatarray = args[0].reshape((5, ))

                    else:
                        raise Exception('Wrong Input type')

                elif args[0].ndim == 2:

                    if q_shape[0] == 4 or q_shape[0] == 5:

                        if q_shape[0] == 4:
                            quat1 = args[0]
                            v_type = np.ones((1, q_shape[1]))
                            quatarray = np.vstack((quat1, v_type))

                        else:
                            quatarray = args[0]

                    elif q_shape[1] == 4 or q_shape[1] == 5:

                        if q_shape[1] == 4:
                            quat1 = args[0].transpose()
                            v_type = np.ones((1, q_shape[0]))
                            quatarray = np.row_stack((quat1, v_type))

                        else:
                            quatarray = args[0].transpose()

                    else:
                        raise Exception('Wrong Input type')

                else:
                    raise Exception('Wrong Input type')

                obj = quatarray.view(cls)

            else:
                err_str1 = 'Wrong Input type: Must be a'
                err_str2 = 'Quaternion or vector3d or numpy array'
                err_str = err_str1+err_str2
                raise Exception(err_str)
        ##############################################################
        elif nargs == 4:
            if gmt.isnumeric(args[0]):
                if np.size(args[0]) == 1:
                    va = np.asarray(args[0]).reshape((np.size(args[0])))
                    vb = np.asarray(args[1]).reshape((np.size(args[1])))
                    vc = np.asarray(args[2]).reshape((np.size(args[2])))
                    vd = np.asarray(args[3]).reshape((np.size(args[3])))
                    v_type = np.ones((np.size(args[3])))
                else:
                    va = np.asarray(args[0]).reshape((1, np.size(args[0])))
                    vb = np.asarray(args[1]).reshape((1, np.size(args[1])))
                    vc = np.asarray(args[2]).reshape((1, np.size(args[2])))
                    vd = np.asarray(args[3]).reshape((1, np.size(args[3])))
                    v_type = np.ones((1, np.size(args[3])))
                quatarray = np.concatenate((va, vb, vc, vd, v_type), axis=0)
                obj = quatarray.view(cls)

            else:
                raise Exception('Wrong Input Types')
        ##############################################################
        elif nargs == 5:
            if gmt.isnumeric(args[0]):
                if np.size(args[0]) == 1:
                    va = np.asarray(args[0]).reshape((np.size(args[0])))
                    vb = np.asarray(args[1]).reshape((np.size(args[1])))
                    vc = np.asarray(args[2]).reshape((np.size(args[2])))
                    vd = np.asarray(args[3]).reshape((np.size(args[3])))
                    v_type = np.asarray(args[4]).reshape((np.size(args[4])))
                else:
                    va = np.asarray(args[0]).reshape((1, np.size(args[0])))
                    vb = np.asarray(args[1]).reshape((1, np.size(args[1])))
                    vc = np.asarray(args[2]).reshape((1, np.size(args[2])))
                    vd = np.asarray(args[3]).reshape((1, np.size(args[3])))
                    v_type = np.asarray(args[4]).reshape((1, np.size(args[4])))
                if np.any(np.logical_and(v_type != 1, v_type != -1)):
                    raise Exception('Wrong Input for v_type')
                quatarray = np.concatenate((va, vb, vc, vd, v_type), axis=0)
                obj = quatarray.view(cls)
            else:
                raise Exception('Wrong Input Types')

        else:
            raise Exception('Wrong Number of Arguments')

        return obj

        # ##############################################################
        #     elif type(args[0]) is vec3d.vector3d:
        #         va = np.zeros((1, vec3d.get_size(args[0])))
        #         vec = vec3d.double(args[0])
        #         v_type = np.ones((1, vec3d.get_size(args[0])))
        #         quatarray = np.concatenate((va, vec, v_type), axis=0)
        #         obj = quatarray.view(cls)
        # ##############################################################

        # ##############################################################
        # elif nargs == 2:
        #     if isinstance(args[1], vec3d.vector3d):
        #         vec_array = vec3d.double(args[1])
        #         v_shape = np.shape(vec_array)
        #     elif type(args[1]) is np.ndarray:
        #         v_shape = np.shape(args[1])
        #         if args[1].ndim == 1:
        #             if v_shape[0] == 3:
        #                 vec_array = args[0].reshape((1, 3))
        #                 v_shape = np.shape(vec_array)
        #             elif v_shape[0] == 4:
        #                 quat1 = args[0].reshape((1, 4))
        #                 v_shape = np.shape(quat1)
        #             else:
        #                 raise Exception('Wrong Input type')
        #         elif args[1].ndim == 2:
        #             if v_shape[0] == 3 or v_shape[0] == 4:
        #                 if v_shape[0] == 3:
        #                     vec_array = args[1]
        #                     v_shape = np.shape(vec_array)
        #                 else:
        #                     quat1 = args[1]
        #                     v_shape = np.shape(quat1)
        #             elif v_shape[1] == 3 or v_shape[1] == 4:
        #                 if v_shape[1] == 3:
        #                     vec_array = args[1].transpose()
        #                     v_shape = np.shape(vec_array)
        #                 else:
        #                     quat1 = args[1].transpose()
        #                     v_shape = np.shape(quat1)
        #             else:
        #                 raise Exception('Wrong Input type')
        #         else:
        #             raise Exception('Wrong Input type')
        #     if np.size(args[0]) == 1:
        #         va = args[0]*np.ones((1, v_shape[1]))
        #     elif np.size(args[0]) == v_shape[1]:
        #         va = np.reshape(args[0], (1, v_shape[1]))
        #     else:
        #         raise Exception('Wrong Input type args[0]')
        #
        #     if isinstance(args[1], vec3d.vector3d):
        #         v_type = np.ones((1, v_shape[1]))
        #         quatarray = np.concatenate((va, vec_array, v_type), axis=0)
        #     elif type(args[1]) is np.ndarray:
        #         if v_shape[0] == 3:
        #             v_type = np.ones((1, v_shape[1]))
        #             quatarray = np.concatenate((va, vec_array, v_type),
        # axis=0)
        #         else:
        #             quatarray = np.concatenate((vec_array, va), axis=0)
        #     obj = quatarray.view(cls)
        #
        # ##############################################################
        # elif nargs == 3:
        #
        #     if isinstance(args[1], vec3d.vector3d):
        #         vec_array = vec3d.double(args[1])
        #         v_shape = np.shape(vec_array)
        #
        #     else:
        #         raise Exception('Wrong Input Type')
        #
        #     if np.size(args[0]) == 1:
        #         va = args[0]*np.ones((1, v_shape[1]))
        #
        #     elif np.size(args[0]) == v_shape[1]:
        #         va = np.reshape(args[0], (1, v_shape[1]))
        #
        #     else:
        #         raise Exception('Wrong Input type args[0]')
        #
        #     if np.size(args[2]) == 1:
        #         v_type = args[2]*np.ones((1, v_shape[1]))
        #
        #     elif np.size(args[2]) == v_shape[1]:
        #         v_type = np.reshape(args[2], (1, v_shape[1]))
        #
        #     else:
        #         raise Exception('Wrong Input type args[0]')
        #
        #     quatarray = np.concatenate((va, vec_array, v_type), axis=0)
        #     obj = quatarray.view(cls)
        #
        ##############################################################

    def __init__(self, *args, **kwargs):
        super(Quaternion, self).__init__(*args, **kwargs)

    # def __str__(q):
    #     # super(Quaternion, q).__str__()
    #     """
    #     Display Quaternion Array in human readable format
    #     """
    #     # q = self
    #     if q.ndim == 1:
    #         str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
    #         str1 += ("%f \t %f \t %f \t %f \t %d \n" %
    #                  (q[0, ], q[1, ], q[2, ], q[3, ], q[4, ]))
    #     else:
    #         s1 = np.shape(q)[1]
    #         str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
    #         for ct1 in range(s1):
    #             str1 += ("%f \t %f \t %f \t %f \t %d \n" %
    #                      (q[0, ct1], q[1, ct1],
    #                       q[2, ct1], q[3, ct1], q[4, ct1]))
    #
    #     return str1
    # def __array_finalize__(self, obj):
    #     print type(obj)
    #     if type(obj) is np.ndarray:
    #         print 'Ndarray type'
    #         return
    #     elif type(obj) is Quaternion:
    #         print 'Quaternion Type'
    #         if self.ndim == 1:
    #             print 'Ndim = 1'
    #             if np.size(self) == 5:
    #                 return
    #             else:
    #                 raise Exception('Wrong Quaternion Assignment')
    #         else:
    #             print 'Ndim = 2'
    #             if np.shape(self)[0] == 5:
    #                 return
    #             else:
    #                 raise Exception('Wrong Quaternion Slicing')

########################################################################################################################
####################################################CLASS METHODS#######################################################
########################################################################################################################


def getq0(g):
    """
    Returns the q0 component for each quaternion present in the input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    q0 components stored in a 1-D numpy array of size n.
    """
    if g.ndim == 1:
        return np.array(g[0])
    else:
        return np.array(g[0, :])
# ----------------------------------------------------------------------------------------------------------------------

def getq1(g):
    """
    Returns the q1 component for each quaternion present in the input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    q1 components stored in a 1-D numpy array of size n.
    """
    if g.ndim == 1:
        return np.array(g[1])
    else:
        return np.array(g[1, :])
# ----------------------------------------------------------------------------------------------------------------------

def getq2(g):
    """
    Returns the q2 component for each quaternion present in the input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    q2 components stored in a 1-D numpy array of size n.
    """
    if g.ndim == 1:
        return np.array(g[2])
    else:
        return np.array(g[2, :])
# ----------------------------------------------------------------------------------------------------------------------

def getq3(g):
    """
    Returns the q3 component for each quaternion present in the input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    q3 components stored in a 1-D numpy array of size n.
    """
    if g.ndim == 1:
        return np.array(g[3])
    else:
        return np.array(g[3, :])
# ----------------------------------------------------------------------------------------------------------------------

def get_type(g):
    """
    Returns the rotation type of each quaternion present in the input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    integer values (either +1 or -1) stored in a 1-D numpy array of size n.
    * +1 is returned for proper rotations
    * -1 is returned for improper rotations
    """
    if g.ndim == 1:
        return g[4]
    else:
        return g[4, :]
# ----------------------------------------------------------------------------------------------------------------------

def get_size(g):
    """
    Returns the size of input quaternion array.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    The size of the input quaternion array, n.
    """
    if g.ndim == 1:
        return 1
    else:
        return np.shape(g)[1]
# ----------------------------------------------------------------------------------------------------------------------

def display(q, p_flag=True):
    """
    Returns a string which displays the input quaternion array in human readable format.

    Parameters
    ----------
    g : input quaternion array
    * a quaternion array of size (5 x n)
    p_flag : flag to print the returned display string
    * a boolean with default value == True

    Returns
    -------
    str1: quaternion array in a human readable format

    Notes
    ------
    * The 5 components of each quaternion stored in the array are displayed under 5 columns. Each row represents a
    quaternion.
    * p_flag == True prints and returns the display string
    * p_flag == False returns the display string

    """
    if q.ndim == 1:
        str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
        str1 += ("%f \t %f \t %f \t %f \t %d \n" %
                  (q[0, ], q[1, ], q[2, ], q[3, ], q[4, ]))
        return str1

    s1 = get_size(q)
    str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
    for ct1 in range(s1):
        str1 += ("%f \t %f \t %f \t %f \t %d \n" %
                  (q[0, ct1], q[1, ct1], q[2, ct1], q[3, ct1], q[4, ct1]))
    if p_flag == True:
        print str1
    return str1
# ----------------------------------------------------------------------------------------------------------------------

def antipodal(q1, tol=1e-12):
    """
    Returns the antipodal (or equivalent) quaternions such that q0 component is positive.

    Parameters
    ----------
    q1 : input quaternion array
    * a quaternion array of size (5 x n)
    tol : tolerance to overcome floating point error
    * a float, default value is 1e-08

    Returns
    -------
    A quaternion array of size (5 x n)

    Notes
    -----
    * Antipodal quaternion for quaternions representing Pi rotations (q0==0) are obtained with a convention.
    * If q0 == 0, q3 >0; if q0 == 0 and q3 == 0, q2 >0; if q == 0, q3 == 0 and q2 == 0, q1 > 0.

    See Also
    --------
    getq0, getq1, getq2, getq3, get_type, get_size
    """
    a1 = getq0(q1)
    b1 = getq1(q1)
    c1 = getq2(q1)
    d1 = getq3(q1)
    e1 = get_type(q1)
    s1 = get_size(q1)

    if s1 == 1:
        if a1 < 0:
            a1 = -a1
            b1 = -b1
            c1 = -c1
            d1 = -d1
        if a1 < tol:
            if d1 < 0:
                b1 = -b1
                c1 = -c1
                d1 = -d1
            if d1 < tol:
                if c1 < 0:
                    b1 = -b1
                    c1 = -c1
                if c1 < tol:
                    if b1 < 0:
                        b1 = -b1

    else:
        ind1 = np.where(a1 < 0)[0]
        a1[ind1] = -a1[ind1]
        b1[ind1] = -b1[ind1]
        c1[ind1] = -c1[ind1]
        d1[ind1] = -d1[ind1]

        ind2 = np.where(a1 < tol)[0]
        b2 = b1[ind2]; c2 = c1[ind2]; d2 = d1[ind2]
        ind3 = np.where(d2 < 0)[0]
        b2[ind3] = -b2[ind3]
        c2[ind3] = -c2[ind3]
        d2[ind3] = -d2[ind3]


        ind4 = np.where(d2 < tol)[0]
        b3 = b2[ind4]; c3 = c2[ind4]
        ind5 = np.where(c3 < 0)[0]
        b3[ind5] = -b3[ind5]
        c3[ind5] = -c3[ind5]



        ind6 = np.where(c3 < tol)[0]
        b4 = b3[ind6]
        ind7 = np.where(b4 < 0)[0]
        b4[ind7] = -b4[ind7]

        b3[ind6] = b4
        b2[ind4] = b3; c2[ind4] = c3
        b1[ind2] = b2; c1[ind2] = c2; d1[ind2] = d2

    # ind8 = np.where(b1[ind6] < tol)
    # b1[ind8] = -b1[ind8]
    # c1[ind8] = -c1[ind8]
    # d1[ind8] = -d1[ind8]
    #
    # ind8 = np.where(b1[ind2] < tol)
    # ind3 = np.where(np.abs(b1[ind2]) < tol)
    # ind5 = np.where(np.abs(c1[ind3]) < tol)
    #
    #
    # if np.abs(d1) < tol:
    #     if np.abs(c1) < tol:
    #         if b1 < 0:
    #             a1[ind2] = a1[ind2]
    #             b1[ind2] = -b1[ind2]
    #             c1[ind2] = -c1[ind2]
    #             d1[ind2] = -d1[ind2]
    #     else:
    #         if c1 < 0:
    #             a1[ind1] = a1[ind1]
    #             b1[ind1] = -b1[ind1]
    #             c1[ind1] = -c1[ind1]
    #             d1[ind1] = -d1[ind1]
    # else:
    #     ind5 = np.where(d1[ind2] < 0)
    #     b1[ind5] = -b1[ind5]
    #     c1[ind5] = -c1[ind5]
    #     d1[ind5] = -d1[ind5]
    #

    return Quaternion(a1, b1, c1, d1, e1)
# ----------------------------------------------------------------------------------------------------------------------

def inverse(q1):
    """
    Returns the inverse quaternions for a given input quaternion array.

    Parameters
    ----------
    q1 : input quaternion array
    * a quaternion array of size (5 x n)

    Returns
    -------
    A quaternion array of size (5 x n)

    See Also
    --------
    getq0, getq1, getq2, getq3, get_type
    """
    a1 = getq0(q1)
    b1 = getq1(q1)
    c1 = getq2(q1)
    d1 = getq3(q1)
    e1 = get_type(q1)

    return Quaternion(a1, -b1, -c1, -d1, e1)
# ----------------------------------------------------------------------------------------------------------------------

def mtimes(q1, q2):
    """
    Calculates the quaternion multiplication of two input quaternions.

    Parameters
    ----------
    q1, q2 : Two quaternion arrays
    * Allowed options: Either size(q1) == size(q2) or size(q1) == 1 or size(q2) == 1

    Returns
    -------
    A quaternion array of size max(size(q1), size(q2))

    Notes
    -------
    * If size(q1) == size(q2), then each quaternion q1(i) is multiplied with q2(i)
    * If size(q1) > size(q2), then size(q2) must be equal to 1. Each quaternion q1(i)
    is multiplied with q2.
    * If size(q2) > size(q1), then size(q1) must be equal to 1. Each quaternion q2(i)
    is multiplied with q1.

    See Also
    --------
    antipodal, getq0, getq1, getq2, getq3, get_type, get_size
    """
    a1 = getq0(q1); b1 = getq1(q1); e1 = get_type(q1)
    c1 = getq2(q1); d1 = getq3(q1); s1 = get_size(q1)

    a2 = getq0(q2); b2 = getq1(q2); s2 = get_size(q2)
    c2 = getq2(q2); d2 = getq3(q2); e2 = get_type(q2)

    if s1 != s2:
        if s1 == 1:
            shp = np.shape(a2); a1 = a1*np.ones(shp)
            b1 = b1*np.ones(shp); c1 = c1*np.ones(shp)
            d1 = d1*np.ones(shp); e1 = e1*np.ones(shp)

        elif s2 == 1:
            shp = np.shape(a2); a2 = a2*np.ones(shp)
            b2 = b2*np.ones(shp); c2 = c2*np.ones(shp)
            d2 = d2*np.ones(shp); e2 = e2*np.ones(shp)
        else:
            str1 = 'size of q1 = ' + str(s1)
            str2 = 'size of q2 = ' + str(s2)
            raise Exception('Wrong Sizes: \n' + str1 + '\n' + str2)

    a = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c = a1*c2 + c1*a2 + d1*b2 - b1*d2
    d = a1*d2 + d1*a2 + b1*c2 - c1*b2

    quat1 = Quaternion(a, b, c, d, e1*e2)
    return antipodal(quat1)
# ----------------------------------------------------------------------------------------------------------------------

def eq(q1, q2, tol=1e-04):
    """
    Checks whether the two input quaternions are equal or not.
    Two quaternions quat1 and quat2 are equal if
    1) type(quat1) == type(quat2) and
    2) (getq0(quat1*inverse(quat2)) - 1) is less than tolerance (tol)
    
    Parameters
    ----------
    q1, q2: Two quaternion arrays
    * Allowed options: Either size(q1) == size(q2) or size(q1) == 1 or size(q2) == 1
    
    tol: The tolerance to check if two quaternions are equal
    
    Returns
    -------
    A boolean array of size max(size(q1), size(q2))

    Notes
    -------
    * If size(q1) == size(q2), then each quaternion q1(i) is checked against q2(i)
    * If size(q1) > size(q2), then size(q2) must be equal to 1. Each quaternion q1(i) is checked
    against q2
    * If size(q2) > size(q1), then size(q1) must be equal to 1. Each quaternion q2(i) is checked
    against q1

    See Also
    --------
    getq0, get_type
    """

    s1 = get_size(q1);  s2 = get_size(q2);
    e1 = get_type(q1); e2 = get_type(q2);
    if s1 != s2:
        if s1 == 1:
            e1 = e1*np.ones((s2, ))
        elif s2 == 1:
            e2 = e2*np.ones((s1, ))
        else:
            str1 = 'size of q1 = ' + str(s1)
            str2 = 'size of q2 = ' + str(s2)
            raise Exception('Wrong Sizes: \n' + str1 + '\n' + str2)


    q3 = mtimes(q1, inverse(q2))
    # print "q3 = \n"; print display(q3);
    # print "q3[0] cond = \n" + str((abs(q3[0] -1) < tol)); 
    # print "Type cond = \n" + str((e1 - e2 == 0));

    return (abs(q3[0] -1) < tol) & (e1 - e2 == 0)
# ----------------------------------------------------------------------------------------------------------------------

def quat2mat(q):
    """
    Converts the input quaternion array to a rotation matrix array.

    Parameters
    ----------
    q: input quaternion
    *  quaternion array of size (5 x n)

    Returns
    ----------
    g: rotation matrix array
    * numpy array of size (n x 3 x 3)

    See Also
    --------
    get_size, getq0, getq1, getq2, getq3, get_type
    """
    sz = get_size(q)
    q0 = getq0(q)
    q1 = getq1(q)
    q2 = getq2(q)
    q3 = getq3(q)
    qt = get_type(q)

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
# ----------------------------------------------------------------------------------------------------------------------

def mat2quat(mat, rot_type='proper'):
    """
    Converts the input rotation matrix array or list, to a quaternion array.

    Parameters
    -----------
    mat: rotation matrices, this input maybe given in the following two ways
    * numpy array of size (n x 3 x 3)
    * python list with n elements, each element is a rotation matrix with shape (3 x 3)

    rot_type: string with either of the following values, 'proper' or 'improper'
    * 'improper' if there is a possibility of having improper matrices in the input
    * 'proper' if all the rotation matrices in the input are of proper rotation type
    * default value is 'proper'

    Returns
    --------
    q: quaternion
    * quaternion array of size (5 x n)

    See Also
    ---------
    tools.vrrotmat2vec
    """
    from tools import vrrotmat2vec as vrrotmat2vec
    ax_ang = vrrotmat2vec(mat, rot_type)

    t1_vecs = ax_ang[:3].T
    new_col = np.linalg.norm(t1_vecs, axis=1)
    t1_vecs_norm = np.array([new_col, ]*3).T
    t1_vecs = np.divide(t1_vecs, t1_vecs_norm)
    ax_ang_norm = t1_vecs.T

    q0 = np.cos(ax_ang[3, :]/2)
    q1 = ax_ang_norm[0, :]*np.sin(ax_ang[3, :]/2)
    q2 = ax_ang_norm[1, :]*np.sin(ax_ang[3, :]/2)
    q3 = ax_ang_norm[2, :]*np.sin(ax_ang[3, :]/2)
    qtype = ax_ang[4, :]

    q = Quaternion(q0, q1, q2, q3, qtype)

    return q
# ----------------------------------------------------------------------------------------------------------------------
def ctranspose(g):
    """
    """
    if g.ndim == 1:
        g[1] = -g[1]
        g[2] = -g[2]
        g[3] = -g[3]
    else:
        g[1, :] = -g[1, :]
        g[2, :] = -g[2, :]
        g[3, :] = -g[3, :]
    return g


# def angle(g1, *args):
#     """
#     calcualtes the rotational angle between rotations q1 and q2
#     Syntax
#     omega = angle(q)
#     omega = angle(q1,q2)
#     Input
#     q1, q2 - @Quaternion
#
#     Output
#     omega  - double
#     """
#     nargs = len(args)
#     if nargs == 0:
#         omega = 2*np.acos(abs(dot(getq0(g1))))
#     elif nargs == 1:
#         g2 = args[0]
#         omega = 2*np.acos(abs(dot(g1, g2)))
#     else:
#         raise Exception('Wrong number of Inputs')
#     return omega
#
#
# def cross(g1, g2, g3):
#     """
#     """
#     a1 = getq0(g1)
#     b1 = getq1(g1)
#     c1 = getq2(g1)
#     d1 = getq3(g1)
#
#     a2 = getq0(g2)
#     b2 = getq1(g2)
#     c2 = getq2(g2)
#     d2 = getq3(g2)
#
#     a3 = getq0(g3)
#     b3 = getq1(g3)
#     c3 = getq2(g3)
#     d3 = getq3(g3)
#
#     # Calculate cross product
#     a = b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1
#     b = a1*c3*d2 - a1*c2*d3 + a2*c1*d3 - a2*c3*d1 - a3*c1*d2 + a3*c2*d1
#     c = a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1
#     d = a1*b3*c2 - a1*b2*c3 + a2*b1*c3 - a2*b3*c1 - a3*b1*c2 + a3*b2*c1
#
#     return Quaternion(a, b, c, d)
#
#

# def angle_outer
# def axis (g):
#     """

#     """
#     v = vec3d.vector3d(double(g[:,1:]));
#     v(np.where(getq0(g) < 0)) = -v(np.where(getq0(g) < 0));
#     return vec3d.normalize(v);

# def calcVoronoi
# def char
# def dot_angle
# def dot_outer
# def double
# def end
# def eq (g1,g2)
# def Euler
# def export
# def find
# def ge
# def horzcat
# def length
# def matrix
# def mean
# def mean_CS
# def minus
# def mldivide
# def mpower
# def mrdivide
# def mtimes
# def ndims
# def ne
# def norm
# def normalize
# def numel
# def partition
# def permute
# def plot
# def plus
# def power
# def project2FundamentalRegion
# def qmatrixL (q):
#     a = getq0(q);
#     b = getq1(q);
#     c = getq2(q);
#     d = getq3(q);

#     if get_size(q)==1:
#         Q = np.array([[a, -b, -c, -d],
#                       [b,  a, -d,  c],
#                       [c,  d,  a, -b],
#                       [d, -c,  b,  a]]);
#     else:
# def qq
# def Quaternion
# def rdivide
# def real
# def repmat
# def reshape
# def Rodrigues
# def size
# def subsasgn
# def subsref
# def sum
# def symmetrise
# def times
# def transpose
# def uminus
# def unique
# def vertcat

# def anti_inv(q1):
#     """
#     Return  for the Quaternion
#     """
#     a1 = -getq0(q1)
#     b1 = -getq1(q1)
#     c1 = -getq2(q1)
#     d1 = -getq3(q1)
#     e1 = -get_type(q1)

#     return Quaternion(a1, b1, c1, d1, e1)


# def dot(q1, q2):
#     """
#     Inner Product of Quaternions q1 and q2

#     Parameters
#     ----------
#     q1, q2: @Quaternion

#     Returns
#     -------
#     dot_product: double
#     """
#     a1 = getq0(q1); b1 = getq1(q1); e1 = get_type(q1)
#     c1 = getq2(q1); d1 = getq3(q1); s1 = get_size(q1);

#     a2 = getq0(q2); b2 = getq1(q2); s2 = get_size(q2);
#     c2 = getq2(q2); d2 = getq3(q2); e2 = get_type(q2)

#     if s1 != s2:
#         if s1 == 1:
#             shp = np.shape(a2); a1 = a1*np.ones(shp);
#             b1 = b1*np.ones(shp); c1 = c1*np.ones(shp);
#             d1 = d1*np.ones(shp); e1 = e1*np.ones(shp);

#         elif s2 == 1:
#             shp = np.shape(a2); a2 = a2*np.ones(shp);
#             b2 = b2*np.ones(shp); c2 = c2*np.ones(shp);
#             d2 = d2*np.ones(shp); e2 = e2*np.ones(shp);
#         else:
#             raise Exception('Wrong Input Types')

#     d = a1*a2 + b1*b2 + c1*c2 + d1*d2
#     if e1.ndim == 0 and e2.ndim == 0:
#         if e1 == e2:
#             return d
#         else:
#             return np.NaN
#     else:
#         d[np.where(e2 != e1*np.ones(np.shape(e2)))] = np.NaN
#         return float(d)
