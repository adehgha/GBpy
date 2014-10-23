# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

import numpy as np
# import quaternion_v1 as quat

#############################################################
# Test Cases
# Case 1: No Arguments


def test_quaternion(type1, case):
    if type1 == 'Init':
        if case == 0:
            q1 = quat.quaternion()
            print q1
        if case == 5:
            ctn = 10
            a1 = np.random.rand(ctn)
            b1 = np.random.rand(ctn)
            c1 = np.random.rand(ctn)
            d1 = np.random.rand(ctn)

            e1 = np.ones(ctn)

            q1 = quat.quaternion(a1, b1, c1, d1, e1)

            print np.shape(q1)
            print q1
    if type1 == 'Slice':
        if case == 5:
            ctn = 10
            a1 = np.random.rand(ctn)
            b1 = np.random.rand(ctn)
            c1 = np.random.rand(ctn)
            d1 = np.random.rand(ctn)

            e1 = np.ones(ctn)

            q1 = quat.quaternion(a1, b1, c1, d1, e1)
            q2 = q1[:, 0:3]

            print q1

            q1

            print np.shape(q2)
            print q2

# type1 = 'Init'; case = 5
##type1 = 'Slice'; case = 5
##test_quaternion(type1, case)

# # Case 2: One Argument 2
# P1 = np.random.rand(4,3); v1 = vec3d.vector3d(P1); q1 = quaternion(v1);
# # Case 2: One Argument quaternion
# q2 = quaternion(q1[:,0]);
# # Case 2: One argument
# q3 = quaternion(np.random.rand(4,))
# q3 = quaternion(np.random.rand(4,3));

# # Case 3: Two Arguments
# # Vector3d and scalar
# q3 = quaternion(0.5, v1);
# P1 = np.random.rand(10,1);
# V2 = vec3d.vector3d(np.random.rand(10,3));
# q3 = quaternion(P1, V2);

# P1 = np.random.rand(10,);
# V2 = vec3d.vector3d(np.random.rand(10,3));
# q3 = quaternion(P1, V2);


# P1 = np.random.rand(10,);
# V2 = np.random.rand(10,3);
# q3 = quaternion(P1, V2);

# # Case 3: 4 arguments
# V1 = np.random.rand(10,);
# V2 = np.random.rand(10,);
# V3 = np.random.rand(10,);
# V4 = np.random.rand(10,);
# q3 = quaternion(V1, V2, V3, V4);

# print q3;
##############################################################
