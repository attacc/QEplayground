#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Calculate frozen phonon starting from pwscf input and dynmat output
#

from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.pwout  import *
from QEplayground.units  import autime2s,amu2au,thz2cm1,au2kelvin
from QEplayground.utilities import bose

import math

def second_deriv(qe_input, qe_dyn, delta, r_order=2, modes=None,T=0.0):
    #
    # Check orthogonality 
    #
    if(not qe_dyn.check_orthogonality()):
        print(" ERROR ERROR ERROR!! ")
        print(" Use the dynamical matrix eigenvectors as input!! ")
        print(" Not the one normalized with the masses!! ")
        exit(1)

    masses=qe_input.get_masses()
    qe_dyn.normalize_with_masses(masses)

    string="\n\n* * * 2nd derivaves calculations * * *\n\n"

    if modes == None:
        modes = range(3, qe_dyn.nmodes) #skyp acustic modes at q=0

    sigma_file=open("sigma.dat","w")

    for im in modes:
        print(" Calculating mode %d .... " % (im+1))
        w_au = qe_dyn.eig[0,im]*(2.0*math.pi)/thz2cm1*autime2s*1e12
        q_0  = 1.0/math.sqrt(2.0*w_au)
        q_T  = q_0*math.sqrt(1.0+2.0*bose(w_au,T/au2kelvin))

        sigma=q_T**2

        sigma_file.write("Sigma(M"+str(im)+")= "+str(sigma)+"\n")

        if r_order == 1:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            qe_right.write(qe_input.filename+"_M"+str(im).zfill(3)+"_R1")


        elif r_order == 2 or r_order == 3:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            qe_left =qe_dyn.generate_displacement(0, im, -delta)
            #
            qe_left.write(qe_input.filename+"_M"+str(im).zfill(3)+"_L1")
            #
            qe_right.write(qe_input.filename+"_M"+str(im).zfill(3)+"_R1")

            if r_order == 3:
                der2_large=der2
                qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
                qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
                #
                qe_left.write(qe_input.filename+"_M"+str(im).zfill(3)+"_L2")
                #
                qe_right.write(qe_input.filename+"_M"+str(im).zfill(3)+"_R2")

    sigma_file.close()


