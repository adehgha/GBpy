# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
# import vector3d as vec3d
import geometry_tools as gmt


######################################################################
class quaternion(np.ndarray):
    """
    The quaternion class is defined to define the quaternion parameterization
    of rotation and their operations
    ...

    Attributes
    ----------
    quaternion: numpy array
        5 x n dimensions

    """

    # 1) The first four rows are the components q0, q1, q2, q3
    # 2) The fifth component is 1 or -1:
    #     a)  1 corresponds to proper rotation
    #     b) -1 corresponds to improper rotation

    # Possible Inputs:
    # 0) nargs = 0: Output - Empty Array
    # 1) nargs = 1:
    #      a) Quaternion
    #      b) 3d vector array
    #      c) Numpy Array
    # **2 and 3 Still Need Implementation**
    # 2) nargs = 2:
    #      a) q0 array (Constant or 1xn); Vector 3d Array
    #      b) 1xn Rotation Type array; Numpy Array
    # 3) nargs = 3:
    #     a) q0 array (Constant or 1xn)
    #     b) Vector 3d Array
    #     c) 1xn Rotation Type Array
    # 4) nargs = 4: a, b, c, d
    # 5) nargs = 5: a, b, c, d, Type


    def __new__(cls, *args):
        nargs = len(args)

        ##############################################################
        if nargs == 0:
            quatarray = np.empty((5, ))
            quatarray[4] = 1
            obj = quatarray.view(cls)

        ##############################################################
        elif nargs == 1:

            if type(args[0]) is quaternion:
                quatarray = args[0]

                if np.shape(args[0])[0] == 5:
                    obj = quatarray
                else:
                    raise Exception('Wrong Quaternion Slice')

        #     elif type(args[0]) is vec3d.vector3d:
        #         va = np.zeros((1, vec3d.get_size(args[0])))
        #         vec = vec3d.double(args[0])
        #         v_type = np.ones((1, vec3d.get_size(args[0])))
        #         quatarray = np.concatenate((va, vec, v_type), axis=0)
        #         obj = quatarray.view(cls)

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
                            v_type = np.ones((1, q_shape[1]))
                            quatarray = np.vstack((quat1, v_type))

                        else:
                            quatarray = args[0].transpose()

                    else:
                        raise Exception('Wrong Input type')

                else:
                    raise Exception('Wrong Input type')

                obj = quatarray.view(cls)

            else:
                err_str1 = 'Wrong Input type: Must be a'
                err_str2 = 'quaternion or vector3d or numpy array'
                err_str = err_str1+err_str2
                raise Exception(err_str)

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

    def __init__(self, *args, **kwargs):
        pass

    # def __array_finalize__(self, obj):
    #     print type(obj)
    #     if type(obj) is np.ndarray:
    #         print 'Ndarray type'
    #         return
    #     elif type(obj) is quaternion:
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

    def __str__(self):
        """
        Display Quaternion Array in human readable format
        """
        q = self
        if q.ndim == 1:
            str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
            str1 += ("%f \t %f \t %f \t %f \t %d \n" %
                     (q[0, ], q[1, ], q[2, ], q[3, ], q[4, ]))
        else:
            s1 = np.shape(q)[1]
            str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \t \t type \n'
            for ct1 in range(s1):
                str1 += ("%f \t %f \t %f \t %f \t %d \n" %
                         (q[0, ct1], q[1, ct1],
                          q[2, ct1], q[3, ct1], q[4, ct1]))

        return str1


def double(g):
    """
    Convert to numpy array to use all the numpy array manipulation routines

    Parameters
    ----------------
    g: quaternion class
        quaternions

    Returns:
    -----------
    g viewed as a numpy array
    """
    return g.view(np.ndarray)


def getq0(g):
    """
    Return the q0 components of the quaternions
    """
    g = double(g)
    if g.ndim == 1:
        return g[0]
    else:
        return g[0, :]


def getq1(g):
    """
    Return the q1 components of the quaternions
    """
    g = double(g)
    if g.ndim == 1:
        return g[1]
    else:
        return g[1, :]


def getq2(g):
    """
    Return the q2 components of the quaternions
    """
    g = double(g)
    if g.ndim == 1:
        return g[2]
    else:
        return g[2, :]


def getq3(g):
    """
    Return the q3 components of the quaternions
    """
    g = double(g)
    if g.ndim == 1:
        return g[3]
    else:
        return g[3, :]


def get_type(g):
    """
    Return the rotation type of the quaternions
    (proper = 1, improper = -1)
    """
    g = double(g)
    if g.ndim == 1:
        return g[4]
    else:
        return g[4, :]


def get_size(g):
    """
    Return the size of the quaternion array
    """
    if g.ndim == 1:
        return 1
    else:
        return np.shape(g)[1]


def antipodal(q1):
    """
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
    else:
        ind1 = np.where(a1 < 0)
        a1[ind1] = -a1[ind1]
        b1[ind1] = -b1[ind1]
        c1[ind1] = -c1[ind1]
        d1[ind1] = -d1[ind1]

    return quaternion(a1, b1, c1, d1, e1)


def anti_inv(q1):
    """
    """
    a1 = -getq0(q1)
    b1 = -getq1(q1)
    c1 = -getq2(q1)
    d1 = -getq3(q1)
    e1 = -get_type(q1)

    return quaternion(a1, b1, c1, d1, e1)


def dot(q1, q2):
    """
    Inner Product of quaternions q1 and q2

    Parameters
    ----------
    q1, q2: @quaternion

    Returns
    -------
    dot_product: double
    """
    a1 = getq0(q1)
    b1 = getq1(q1)
    c1 = getq2(q1)
    d1 = getq3(q1)
    e1 = get_type(q1)
    s1 = get_size(q1)

    a2 = getq0(q2)
    b2 = getq1(q2)
    c2 = getq2(q2)
    d2 = getq3(q2)
    e2 = get_type(q2)
    s2 = get_size(q2)

    if s1 != s2:
        if s1 == 1:
            shp = np.shape(a2)
            a1 = a1*np.ones(shp)
            b1 = b1*np.ones(shp)
            c1 = c1*np.ones(shp)
            d1 = d1*np.ones(shp)
            e1 = e1*np.ones(shp)

        elif s2 == 1:
            shp = np.shape(a2)
            a2 = a2*np.ones(shp)
            b2 = b2*np.ones(shp)
            c2 = c2*np.ones(shp)
            d2 = d2*np.ones(shp)
            e2 = e2*np.ones(shp)
        else:
            raise Exception('Wrong Input Types')

    d = a1*a2 + b1*b2 + c1*c2 + d1*d2
    if e1.ndim == 0 and e2.ndim == 0:
        if e1 == e2:
            return d
        else:
            return np.NaN
    else:
        d[np.where(e2 != e1*np.ones(np.shape(e2)))] = np.NaN
        return d


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


def mtimes(q1, q2):
    """
    Quaternion muliplication q1*q2
    """

    a1 = getq0(q1)
    b1 = getq1(q1)
    c1 = getq2(q1)
    d1 = getq3(q1)
    e1 = get_type(q1)
    s1 = get_size(q1)

    a2 = getq0(q2)
    b2 = getq1(q2)
    c2 = getq2(q2)
    d2 = getq3(q2)
    e2 = get_type(q2)
    s2 = get_size(q2)

    if s1 != s2:
        if s1 == 1:
            shp = np.shape(a2)
            a1 = a1*np.ones(shp)
            b1 = b1*np.ones(shp)
            c1 = c1*np.ones(shp)
            d1 = d1*np.ones(shp)
            e1 = e1*np.ones(shp)

        elif s2 == 1:
            shp = np.shape(a2)
            a2 = a2*np.ones(shp)
            b2 = b2*np.ones(shp)
            c2 = c2*np.ones(shp)
            d2 = d2*np.ones(shp)
            e2 = e2*np.ones(shp)
        else:
            raise Exception('Wrong Input Types')

    a = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c = a1*c2 + c1*a2 + d1*b2 - b1*d2
    d = a1*d2 + d1*a2 + b1*c2 - c1*b2

    quat1 = quaternion(a, b, c, d, e1*e2)

    return antipodal(quat1)


def eq(q1, q2, tol):
    """
    ? q1 == q2
    """
    # s1 = get_size(q1)
    # s2 = get_size(q2)

    q3 = mtimes(q1, inverse(q2))
    q1_type = get_type(q1)
    q2_type = get_type(q2)

    if np.any(np.all(((abs(abs(getq0(q3)) - 1) < tol),
                      (q1_type == q2_type)), axis=0)):
        return True
    else:
        return False

    # if s1 == s2:
    #     if np.all(abs(abs(dot(q1, q2)-1) < tol):
    #         return True
    #     else:
    #         return False
    # else:
    #     if s1 == 1 or s2 == 1:
    #         if np.any(abs(dot(q1, q2)-1) < tol):
    #             return True
    #         else:
    #             return False
    #     else:
    #         err_str1 = 'Input dimensions incompatible: \n Either: \n'
    #         err_str2 = 'size(q1) == size(q2) or \n'
    #         err_str3 = 'size(q1) == 1 or size(q2)==5 \n'
    #         err_str4 = 'Instead'
    #         err_str5 = 'size(q1) == %d or size(q2)== %d \n' % (s1, s2)
    #         err_str = err_str1 + err_str2 + err_str3 + err_str4 + err_str5
    #         raise Exception(err_str)


def display(q):
    """
    Display Quaternion Array in human readable format
    """

    s1 = get_size(q)
    str1 = 'Quaternion: \n q0 \t \t q1 \t \t q2 \t \t q3 \n'
    for ct1 in range(s1):
        str1 += ("%f \t %f \t %f \t %f \n"
                 % (q[0, ct1], q[1, ct1], q[2, ct1], q[3, ct1]))

    return str1


def inverse(q1):
    """
    """
    a1 = getq0(q1)
    b1 = getq1(q1)
    c1 = getq2(q1)
    d1 = getq3(q1)
    e1 = get_type(q1)

    return quaternion(a1, -b1, -c1, -d1, e1)


#
# def angle(g1, *args):
#     """
#     calcualtes the rotational angle between rotations q1 and q2
#     Syntax
#     omega = angle(q)
#     omega = angle(q1,q2)
#     Input
#     q1, q2 - @quaternion
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
#     return quaternion(a, b, c, d)
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
# def quaternion
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