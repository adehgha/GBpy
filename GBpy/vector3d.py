# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import geometry_tools as gmt

######################################################################
class vector3d(np.ndarray):
    def __new__(cls, *args):
        nargs = len(args)

        ##############################################################
        if nargs == 0:
            vecarray = np.empty((3, 1))
            obj = np.asarray(vecarray).view(cls)
            # obj = np.asarray(vecarray, dtype=dt).view(cls);

        ##############################################################
        elif nargs == 1:
            if type(args[0]) is vector3d:
                vecarray = args[0]
                if np.shape(args[0])[0] == 3:
                    if args[0].ndim == 1:
                        vecarray=np.reshape(vecarray, (3, 1))
                    obj = vecarray
                else:
                    raise Exception('Wrong vector3d Type')
            elif type(args[0]) is np.ndarray:
                vshape = np.shape(args[0])
                if args[0].ndim == 1:
                    if vshape[0] == 3:
                        vecarray=args[0].reshape((3, 1))
                    else:
                        raise Exception('Wrong Input type')
                elif args[0].ndim == 2:
                    if vshape[0] == 3:
                        vecarray = args[0]
                    elif vshape[1] == 3:
                        vecarray = args[0].transpose()
                    else:
                        raise Exception('Wrong Input type')
                else:
                    raise Exception('Wrong Input type')
                obj = np.asarray(vecarray).view(cls)
                # obj = np.asarray(args[0], dtype=dt).view(cls);
            else:
                errstr1 = 'Wrong Input type: Must be a'
                errstr2 = 'vector3d or numpy array'
                errstr=errstr1+errstr2
                raise Exception(errstr)

        ##############################################################
        elif nargs == 2:
            raise Exception('For spherical coordinates use sph2vec')

        ##############################################################
        elif nargs == 3:
            if gmt.isnumeric(args[0]):
                vx = args[0].reshape([1, np.size(args[0])])
                vy = args[1].reshape([1, np.size(args[1])])
                vz = args[2].reshape([1, np.size(args[2])])
                vecarray = np.concatenate((vx, vy, vz), axis=0)
                obj = np.asarray(vecarray).view(cls)
                # obj = np.asarray(vecarray, dtype=dt).view(cls);
            elif args[0] == 'polar':
                ### sph2vec should returns numpy array
                vecarray = gmt.sph2vec(args[1], args[2])
                obj = np.asarray(vecarray).view(cls)
                # obj = np.asarray(vecarray, dtype=dt).view(cls);
        else:
            raise Exception('Wrong Number of Inputs')
        
        return obj

    def __init__(self, *args, **kwargs):
        pass

# def abs


def dot(v1, v2):
    xx = v1[0, :]*v2[0, :]
    yy = v1[1, :]*v2[1, :]
    zz = v1[2, :]*v2[2, :]

    return xx + yy + zz


def angle(v1, v2):
    a = dot(normalize(v1), normalize(v2))
    return np.acos(a)


def double(vec):
    return vec.view(np.ndarray)

# def angle_outer
# def char(v,*args):
# def check_option():
# def contourf
# def cross


def cross(v1, v2):
    """
    :param v1:
    :param v2:
    :return:
    """
    vx = v1[1, :]*v2[2, :]-v1[2, :]*v2[1, :]
    vy = v1[2, :]*v2[0, :]-v1[0, :]*v2[2, :]
    vz = v1[0, :]*v2[1, :]-v1[1, :]*v2[0, :]

    return vector3d(vx, vy, vz)
# def ctranspose
# def delete_option
# def display
# def dot:
# def dot_outer
# def end
# def eq
# def extract_option
# def find
# def get


def getx(vec):
    return vec[0, :]


def gety(vec):
    return vec[1, :]


def getz(vec):
    return vec[2, :]


def get_size(vec):
    return np.size(vec[0, :])

# def horzcat
# def kernelDensityEstimation
# def length
# def line
# def mean
# def minus
# def mtimes
# def ne
# def norm


def norm(vec):
    n = np.sqrt(np.power(vec[0,:],2) + np.power(vec[1,:],2) + np.power(vec[2,:],2))
    return n

# def normalize


def normalize(vec):
    n = np.tile(norm(vec), (3, 1))
    return vec/n


# def numel
# def orth
# def pathPatala
# def pcolor
# def plot 
# def plus
# def polar
# def project2FundamentalRegion
# def quiver
# def rdivide
# def repmat
# def reshape
# def rotate
# def scatter
# def SchmidTensor
# def set
# def set_option
# def size
# def smooth
# def splitNorthSouth
# def subsasgn
# def subsref
# def sum 
# def surf
# def symmetrise
# def text
# def times
# def transpose
# def uminus
# def unique
# def vector3d
# def vertcat


######################################################################
# # Test Cases

# # # Args = 0
# # c = vector3d();


# # Args = 1, Type 1
# P1 = np.random.rand(3,4); c = vector3d(P1); print c;

# # Args = 1, Type 0
# c1 = vector3d(c); print c1;

# # # Args = 2
# # P1 = np.random.rand(3,4);
# # c2 = vector3d(np.random.rand(3,1),np.random.rand(3,1));

# # Args = 3
# c2 = vector3d(np.random.rand(1,4), np.random.rand(1,4), np.random.rand(1,4));


# # c2 = vector3d(np.random.rand(1,12), np.random.rand(1,12), np.random.rand(1,12));
# # print c2;
