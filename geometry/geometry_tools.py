import numpy as np

def isnumeric(obj):
    try:
        obj + 0
        return True
    except TypeError:
        return False


import vector3d as vec3d
def sph2vec(theta, rho, *args):
    """
    Spherical to Cartesian Coordinates

    Transforms spherical into cartesian coordinates

    ## Syntax
    v = sph2vec(theta,rho)
    v = sph2vec(theta,rho,r)
    [x,y,z] = sph2vec(theta,rho,r)
    #

    ## Input
    theta, rho - spherical coordinates in radians
    r          - radius
    nargout    - Type of output arguments
                   1 - vector3d
                   2 - x,y,z
    #

    ## Output
    v          - vector3d
    x,y,z      - double
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


import quaternion as quat
def idquaternion():
    return quat.quaternion(1, 0, 0, 0)
