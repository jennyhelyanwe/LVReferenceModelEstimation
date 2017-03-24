__author__ = 'zwan145'

"""
functions for transforming points
"""

import scipy


def transformRigid3D(x, t):
    """ applies a rigid transform to list of points x.
    T = (tx,ty,tz,rx,ry,rz)
    """

    X = scipy.vstack((x.T, scipy.ones(x.shape[0]) ))

    T = scipy.array([[1.0, 0.0, 0.0, t[0]], \
                     [0.0, 1.0, 0.0, t[1]], \
                     [0.0, 0.0, 1.0, t[2]], \
                     [0.0, 0.0, 0.0, 1.0]])

    Rx = scipy.array([[1.0, 0.0, 0.0], \
                      [0.0, scipy.cos(t[3]), -scipy.sin(t[3])], \
                      [0.0, scipy.sin(t[3]), scipy.cos(t[3])]])

    Ry = scipy.array([[scipy.cos(t[4]), 0.0, scipy.sin(t[4])], \
                      [0.0, 1.0, 0.0], \
                      [-scipy.sin(t[4]), 0.0, scipy.cos(t[4])]])

    Rz = scipy.array([[scipy.cos(t[5]), -scipy.sin(t[5]), 0.0], \
                      [scipy.sin(t[5]), scipy.cos(t[5]), 0.0], \
                      [0.0, 0.0, 1.0]])

    T[:3, :3] = scipy.dot(scipy.dot(Rx, Ry), Rz)
    return scipy.dot(T, X)[:3, :].T


def transformScale3D(x, S):
    """ applies scaling to a list of points x. S = (sx,sy,sz)
    """
    return scipy.multiply(x, S)


def transformRigidSize3D(x, t):
    return transformScale3D(transformRigid3D(x, t[:6]), t[6])


def transformAffine(x, t):
    """ applies affine transform t (shape = (3,4)) to list of points x
    """
    return scipy.dot(t, scipy.vstack((x.T, scipy.ones(x.shape[0]))))[:3, :].T