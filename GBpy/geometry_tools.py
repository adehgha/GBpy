# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2015,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.
# How to cite GBpy:
# Banadaki, A. D. & Patala, S. "An efficient algorithm for computing the primitive bases of a general lattice plane",
# Journal of Applied Crystallography 48, 585-588 (2015). doi:10.1107/S1600576715004446


import numpy as np
import quaternion as quat
import vector3d as vec3d
# -----------------------------------------------------------------------------------------------------------


def isnumeric(obj):
    try:
        obj + 0
        return True
    except TypeError:
        return False
# -----------------------------------------------------------------------------------------------------------


def sph2vec(theta, rho, *args):
    """
    Spherical to Cartesian Coordinates

    Transforms spherical into cartesian coordinates

    ## Syntax
    v = sph2vec(theta,rho)
    v = sph2vec(theta,rho,r)
    [x,y,z] = sph2vec(theta,rho,r)

    Parameters
    -----------------
    theta: spherical coordinates in radians
    rho: spherical coordinates in radians
    r: radius

    Returns
    -----------------
    v: vector3d
    x,y,z: double
    """

    nargin = len(args)

    if nargin == 0:
        r = 1
        nargout = 'Array'
    elif nargin == 1:
        if isnumeric(args[0]):
            r = args[0]
            nargout = 'Array'
        else:
            if not(agrs[1] == 'Array' or args[1] == 'vector3d'):
                raise TypeError('Additional input should be a string: ndarray or Vector3d')
            else:
                nargout = args[1]
    elif nargin == 2:
        r = args[0]
        if not(agrs[1] == 'Array' or args[1] == 'vector3d'):
            raise TypeError('Additional input should be a string: ndarray or Vector3d')
        else:
            nargout = args[1]

    if isinstance(theta, int) or isinstance(theta, float):
        theta = np.array([theta])
    if isinstance(theta, int) or isinstance(theta, float):
        rho = np.array([rho])

    x = r*np.sin(theta)*np.cos(rho)
    y = r*np.sin(theta)*np.sin(rho)
    z = r*np.cos(theta)

    if nargout == 'vector3d':
        return vec3d.vector3d(x, y, z)
    elif nargout == 'Array':
        return np.column_stack((x, y, z))
# -----------------------------------------------------------------------------------------------------------


def idquaternion():
    return quat.Quaternion(1, 0, 0, 0, 1)
# -----------------------------------------------------------------------------------------------------------
