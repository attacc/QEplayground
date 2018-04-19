# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
import numpy as np

atol = 1e-6

def vec_in_list(veca,vec_list):
    """ check if a vector exists in a list of vectors
    """
    return np.array([ np.allclose(veca,vecb,rtol=atol,atol=atol) for vecb in vec_list ]).any()

def red2car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    new_pos=[]
    for row in red:
        new_pos.append(np.matmul(row,lat))
    return np.array(new_pos)

def car2red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    new_pos=[]
    for row in car:
        new_pos.append(np.linalg.solve(np.array(lat).T,row))
    return np.array(new_pos)

def rec_lat(lat):
    """
    Calculate the reciprocal lattice vectors
    """
    a1,a2,a3 = np.array(lat)
    v = np.dot(a1,np.cross(a2,a3))
    b1 = 2.0*np.pi*np.cross(a2,a3)/v
    b2 = 2.0*np.pi*np.cross(a3,a1)/v
    b3 = 2.0*np.pi*np.cross(a1,a2)/v
    return np.array([b1,b2,b3])
