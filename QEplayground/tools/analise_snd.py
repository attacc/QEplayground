#!/usr/bin/python3
#
# Copyright (c) 2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
import numpy as np

eps_name="o.eps_q1_ip"
noph_folder="NOPH"
prefix="bn"
sigma_fname="sigma.dat"

sigma_file=open(sigma_fname,'r')

modes=numpy.genfromtxt(sigma_fname,delimiter=",")
print(modes[:,0])




