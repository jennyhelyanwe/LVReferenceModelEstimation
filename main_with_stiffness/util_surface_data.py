__author__ = 'zwan145'
"""
This function contains subroutines designed to calculate transformation between two sets of data with
no correspondance
--Author: Vicky Wang
--Date: 10/11/2011
"""

import os
import scipy
import numpy
import math
import re
import string
from scipy import array, concatenate
from scipy import spatial

from scipy.spatial import cKDTree
from scipy.optimize import leastsq, fmin
from util_transform import *
import util_fitting_tools

# ==================================================================#
def readSurfaceData(fileName):
    try:
        file = open(fileName, 'r')
    except IOError:
        print 'ERROR: readSurfaceData: unable to open', fileName
        return
    ## Read a single line at a time
    no_line = 0
    data = file.readline()
    no_line = no_line + 1
    cont_read = 1

    ## Initialise an array to store data
    surface_data = []

    ## Start to read the file
    while cont_read == 1:
        check_match = re.match(data.strip(), "<Points>")
        if (check_match != None):
            ## Start to read in the actual data
            ## First line is junk
            junk = file.readline()
            no_line = no_line + 1
            while cont_read == 1:
                data = file.readline()
                no_line = no_line + 1
                check_match = re.match(data.strip(), "</DataArray>")
                if (check_match != None):
                    cont_read = 0
                    break
                if (len(data) > 50):
                    surface_data.append([float(string) for string in data.split()[0:3]])
                    surface_data.append([float(string) for string in data.split()[3:6]])
                else:
                    surface_data.append([float(string) for string in data.split()[0:3]])
        else:
            data = file.readline()
            no_line = no_line + 1

    surface_data = array(surface_data)

    return surface_data


#==================================================================#

#======================================================================#
def readIpdata(fileName):
    """ reads ipdata file and returns the x y z coords on data points
    and the header if there is one
    """

    try:
        file = open(fileName, 'r')
    except IOError:
        print 'ERROR: readIpdata: unable to open', fileName
        return

    lines = file.readline()
    lines = file.readlines()

    coords = []
    for l in lines:
        coords.append([float(string) for string in l.split()[1:4]])

    coords = array(coords)

    return coords


#======================================================================#

#==========================================================================#
def constructTransformationMatrix(t):
    """ construct full transformation matrix T from transformation vector t.

    """
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

    return T


#==========================================================================#

#==========================================================================#
def inverseTransformRigid3D(x, t):
    """ applies a inverse rigid transform to list of points x.
    T = (tx,ty,tz,rx,ry,rz)
    """

    X = scipy.vstack((x.T, scipy.ones(x.shape[0]) ))

    T = constructTransformationMatrix(t)
    #print "T Matrix"
    #print T

    Tinverse = scipy.linalg.inv(T)
    #print "T Matrix Inverse "
    #print Tinverse
    return scipy.dot(Tinverse, X)[:3, :].T


#==========================================================================#

#==========================================================================#
def transformScale3D(x, S):
    """ applies scaling to a list of points x. S = (sx,sy,sz)
    """
    #print 'Scaling Matrix'
    #print S
    return scipy.multiply(x, S)


#==========================================================================#

#==========================================================================#
def transformRigid3DFinal(x, t):
    """ applies a rigid transform to list of points x.
    T = (tx,ty,tz,rx,ry,rz)
    """

    X = scipy.vstack((x.T, scipy.ones(x.shape[0]) ))
    T = constructTransformationMatrix(t)

    print 'Final Transformation Matrix from Canine to UPF with Anisotropic Scaling'
    print T

    return scipy.dot(T, X)[:3, :].T


#==========================================================================#

#======================================================================#
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


#==========================================================================#

#======================================================================#
def fitDataRigidScaleNoCorr_TwoSurfaces(X1, X2, data1, data2, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])

    data1Tree = cKDTree(data1)
    data2Tree = cKDTree(data2)
    X1 = scipy.array(X1)
    X2 = scipy.array(X2)

    def obj(t):
        x1R = transformRigid3D(X1, t[:6])
        x1RS = transformScale3D(x1R, scipy.ones(3) * t[6])
        d1 = data1Tree.query(list(x1RS))[0]
        x2R = transformRigid3D(X2, t[:6])
        x2RS = transformScale3D(x2R, scipy.ones(3) * t[6])
        d2 = data2Tree.query(list(x2RS))[0]
        d = concatenate((d1, d2), 0)
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    X1Opt = transformRigid3D(X1, tOpt[:6])
    X2Opt = transformRigid3D(X2, tOpt[:6])
    X1Opt = transformScale3D(X1Opt, tOpt[6:])
    X2Opt = transformScale3D(X2Opt, tOpt[6:])

    return tOpt, X1Opt, X2Opt


#==========================================================================#

#==========================================================================#
def fitDataRigidAnisotropicScaleNoCorr(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    dataTree = cKDTree(data)
    X = scipy.array(X)

    def obj(t):
        xR = transformRigid3D(X, t[:6])
        xRS = transformScale3D(xR, t[6:])
        d = dataTree.query(list(xRS))[0]
        #~ print d.mean()
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3D(X, tOpt[:6])
    XOpt = transformScale3D(XOpt, tOpt[6:])

    return tOpt, XOpt


#==========================================================================#

#==========================================================================#
def fitDataRigidScaleNoCorr_ModelTree(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])

    data = scipy.array(data)

    def obj(t):
        xR = transformRigid3D(X, t[:6])
        xRS = transformScale3D(xR, scipy.ones(3) * t[6])
        xRSTree = cKDTree(xRS)
        d = xRSTree.query(list(data))[0]
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3D(X, tOpt[:6])
    XOpt = transformScale3D(XOpt, tOpt[6:])

    return tOpt, XOpt


#==========================================================================#

#==========================================================================#
def fitDataRigidAnisotropicScaleNoCorr_ModelTree(X, data, xtol=1e-5, maxfev=0, t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    data = scipy.array(data)

    def obj(t):
        xR = transformRigid3D(X, t[:6])
        xRS = transformScale3D(xR, t[6:])
        xRSTree = cKDTree(xRS)
        d = xRSTree.query(list(data))[0]
        return d * d

    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    XOpt = transformRigid3DFinal(X, tOpt[:6])
    XOpt = transformScale3D(XOpt, tOpt[6:])

    return tOpt, XOpt


#==========================================================================#

#==========================================================================#
def fitDataRigidAnisotropicScaleNoCorr_ModelTree_TwoSurfaces(X1, X2, data1, data2, filename, xtol=1e-5, maxfev=0,
                                                             t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    data1 = scipy.array(data1)
    data2 = scipy.array(data2)

    def obj(t):
        x1R = transformRigid3D(X1, t[:6])
        x1SR = transformScale3D(x1R, t[6:])
        x1SRTree = cKDTree(x1SR)
        d1 = x1SRTree.query(list(data1))[0]
        x2R = transformRigid3D(X2, t[:6])
        x2SR = transformScale3D(x2R, t[6:])
        x2SRTree = cKDTree(x2SR)
        d2 = x2SRTree.query(list(data2))[0]
        d = concatenate((d1, d2), 0)
        return d * d


    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    X1Opt = transformRigid3DFinal(X1, tOpt[:6])
    X2Opt = transformRigid3DFinal(X2, tOpt[:6])
    X1Opt = transformScale3D(X1Opt, tOpt[6:])
    X2Opt = transformScale3D(X2Opt, tOpt[6:])

    writeTransformation_Scaling(tOpt[6:], filename)

    print 'Final SSQ'
    Final_SSQ = scipy.sum(obj(tOpt))
    print Final_SSQ

    return tOpt, X1Opt, X2Opt


#==========================================================================#

#==========================================================================#
def fitDataAnisotropicScaleRigidNoCorr_ModelTree_TwoSurfaces(X1, X2, data1, data2, filename, xtol=1e-5, maxfev=0,
                                                             t0=None):
    """ fit list of points X to list of points data by minimising
    least squares distance between each point in X and closest neighbour
    in data
    """
    if t0 == None:
        t0 = scipy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0])

    data1 = scipy.array(data1)
    data2 = scipy.array(data2)

    def obj(t):
        x1S = transformScale3D(X1, t[6:])
        x1RS = transformRigid3D(x1S, t[:6])
        x1RSTree = cKDTree(x1RS)
        d1 = x1RSTree.query(list(data1))[0]
        x2S = transformScale3D(X2, t[6:])
        x2RS = transformRigid3D(x2S, t[:6])
        x2RSTree = cKDTree(x2RS)
        d2 = x2RSTree.query(list(data2))[0]
        d = concatenate((d1, d2), 0)
        return d * d


    tOpt = leastsq(obj, t0, xtol=xtol, maxfev=maxfev)[0]
    X1Opt = transformScale3D(X1, tOpt[6:])
    X2Opt = transformScale3D(X2, tOpt[6:])
    X1Opt = transformRigid3DFinal(X1Opt, tOpt[:6])
    X2Opt = transformRigid3DFinal(X2Opt, tOpt[:6])

    writeTransformation_Scaling(tOpt[6:], filename)

    print 'Final SSQ'
    Final_SSQ = scipy.sum(obj(tOpt))
    print Final_SSQ

    return tOpt, X1Opt, X2Opt


#==========================================================================#

#==========================================================================#
def writeIpdata(d, filename, header=None):
    """ write the coordinates of points in 3D space to ipdata file. Each
    row in d is [x,y,z] of a datapoint. filename is a string, header is
    a string. if ex!=False, uses cmConvert to conver to ex formates. ex
    can be 'data', 'node' or 'both'
    """

    outputFile = open(filename, 'w')
    if header:
        outputFile.write(header + '\n')

    n = 1
    for i in d:
        outputFile.write(
            " " + str(n) + "\t   %(x)3.10f\t%(y)3.10f\t%(z)3.10f\t1.0 1.0 1.0\n" % {'x': i[0], 'y': i[1], 'z': i[2]})
        n += 1
    outputFile.close()


#==========================================================================#

#==========================================================================#
def writeTransformation(t, filename):
    """ Write the transformaion to a file
    """
    outputFile = open(filename, 'w')
    string = ""
    print t
    for i in t:
        for j in i:
            string = string + "%s," % ( j)
    outputFile.write(string)
    outputFile.close()


#==========================================================================#
#==========================================================================#
def writeTransformation_Scaling(t, filename):
    """ Write the transformaion to a file
    """
    outputFile = open(filename, 'w')

    for i in t:
        print i
        outputFile.write("%(x)3.10f," % {'x': i})

    outputFile.close()


#==========================================================================#
#==========================================================================#
def writeTransformation_Rotation(t, filename):
    """ Write the transformaion to a file
    """
    outputFile = open(filename, 'a+')
    outputFile.write(str(t))
    outputFile.close()

#==========================================================================#
