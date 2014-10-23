import numpy as np


def eq(m1, m2, tol):
    """
    Check if the two rotation matrices are the same
    """
    if m1.ndim == 2 and m2.ndim == 2:
        m = abs(m1 - m2)

        if np.amax(m) < tol:
            return True
        else:
            return False
    elif m1.ndim == 2:
        msz = np.shape(m2)[0]
        tmat1 = m1.reshape((1, 9))
        tmat2 = np.tile(tmat1, (msz, 1))
        tmat3 = tmat2.reshape(msz, 3, 3)

        m = abs(tmat3 - m2)
        max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
        if np.any(max1):
            return True
        else:
            return False

    elif m2.ndim == 2:
        msz = np.shape(m1)[0]
        tmat1 = m2.reshape(msz, (1, 9))
        tmat2 = np.tile(tmat1, (msz, 1))
        tmat3 = tmat2.reshape(msz, 3, 3)

        m = abs(m1 - tmat3)
        max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
        if np.any(max1):
            return True
        else:
            return False
    else:
        if np.shape(m1)[0] == np.shape(m2)[0]:
            m = abs(m1 - m2)
            max1 = np.amax(np.amax(m, axis=1), axis=1) < tol
            return np.where(max1)
        else:
            raise Exception('Wrong Input Types')
