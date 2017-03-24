__author__ = 'zwan145'

"""
functions and classes for fitting/registering data points

each function returns the fitted transform parameters and the fitted
data cloud.

Ju Zhang
2011-11-09
"""

import scipy
from scipy.spatial import cKDTree
from scipy.optimize import leastsq, fmin
from util_transform import transformRigid3D, transformScale3D, transformRigidSize3D, transformAffine

# ======================================================================#
# correspondent data fitting										   #
#======================================================================#


def fitTranslation(data, target, xtol=1e-5, maxfev=0, verbose=0):
    """ fits for tx,ty for transforms points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """

    def obj(x):
        dataT = data + x
        d = ((dataT - target) ** 2.0).sum(1)
        return d

    x0 = scipy.array([0.0, 0.0, 0.0])

    if verbose:
        rms0 = scipy.sqrt(obj(x0).mean())
        print 'initial RMS:', rms0

    xOpt = leastsq(obj, x0, xtol=xtol, maxfev=maxfev)[0]

    if verbose:
        rmsOpt = scipy.sqrt(obj(xOpt).mean())
        print 'final RMS:', rmsOpt

    dataFitted = data + xOpt
    return xOpt, dataFitted


def fitRigid(data, target, x0=None, xtol=1e-3, maxfev=0, verbose=0):
    """ fits for tx,ty,tz,rx,ry,rz to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """
    if x0 == None:
        x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def obj(x):
        dataT = transformRigid3D(data, x)
        d = ( ( dataT - target ) ** 2.0 ).sum(1)
        return d

    x0 = scipy.array(x0)
    if verbose:
        rms0 = scipy.sqrt(obj(x0).mean())
        print 'initial RMS:', rms0

    xOpt = leastsq(obj, x0, xtol=xtol, maxfev=maxfev)[0]

    if verbose:
        rmsOpt = scipy.sqrt(obj(xOpt).mean())
        print 'final RMS:', rmsOpt

    dataFitted = transformRigid3D(data, xOpt)
    return xOpt, dataFitted


def fitRigidFMin(data, target, x0=None, xtol=1e-3, maxfev=0, verbose=0):
    """ fits for tx,ty,tz,rx,ry,rz to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order. uses fmin instead of leastsq
    """
    if x0 == None:
        x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def obj(x):
        dataT = transformRigid3D(data, x)
        d = ( ( dataT - target ) ** 2.0 ).sum(1)
        rmsD = scipy.sqrt(d.mean())
        return rmsD

    x0 = scipy.array(x0)
    if verbose:
        rms0 = scipy.sqrt(obj(x0).mean())
        print 'initial RMS:', rms0

    xOpt = fmin(obj, x0, xtol=xtol, maxiter=maxfev)

    if verbose:
        rmsOpt = scipy.sqrt(obj(xOpt).mean())
        print 'final RMS:', rmsOpt

    dataFitted = transformRigid3D(data, xOpt)
    return xOpt, dataFitted


def fitRigidSize(data, target, x0=None, xtol=1e-3, maxfev=0, verbose=0):
    """ fits for tx,ty,tz,rx,ry,rz,s to transform points in data to points
    in target. Points in data and target are assumed to correspond by
    order
    """
    if x0 == None:
        x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

    def obj(x):
        dataT = transformRigidSize3D(data, x)
        d = ( ( dataT - target ) ** 2.0 ).sum(1)
        return d

    x0 = scipy.array(x0)
    if verbose:
        rms0 = scipy.sqrt(obj(x0).mean())
        print 'initial RMS:', rms0

    xOpt = leastsq(obj, x0, xtol=xtol, maxfev=maxfev)[0]

    if verbose:
        rmsOpt = scipy.sqrt(obj(xOpt).mean())
        print 'final RMS:', rmsOpt

    dataFitted = transformRigidSize3D(data, xOpt)
    return xOpt, dataFitted


#======================================================================#
# Non correspondent data fitting									   #
#======================================================================#
def fitDataRigidNoCorr(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    dataTree = cKDTree(data)
    X = scipy.array(X)

    def obj(t):
        x = transformRigid3D(X, t)
        d = dataTree.query(list(x))[0]
        #~ print d.mean()
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3D(X, tOpt)
    finalRMSE = scipy.sqrt(obj(tOpt).mean())
    print 'fitDataRigidEPDP finalRMSE:', finalRMSE

    return tOpt, XOpt


def fitDataTranslateNoCorr(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0])

    dataTree = cKDTree(data)
    X = scipy.array(X)

    def obj(t):
        x = transformRigid3D(X, scipy.hstack((t, [0.0, 0.0, 0.0])))
        d = dataTree.query(list(x))[0]
        #~ print d.mean()
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3D(X, scipy.hstack((tOpt, [0.0, 0.0, 0.0])))

    return tOpt, XOpt


def fitDataRigidScaleNoCorr(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])

    dataTree = cKDTree(data)
    X = scipy.array(X)

    def obj(t):
        xR = transformRigid3D(X, t[:6])
        xRS = transformScale3D(xR, scipy.ones(3) * t[6])
        d = dataTree.query(list(xRS))[0]
        #~ print d.mean()
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3D(X, tOpt[:6])
    XOpt = transformScale3D(XOpt, tOpt[6:])

    return tOpt, XOpt
