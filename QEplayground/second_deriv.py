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
from QEplayground.units  import autime2s,amu2au,thz2cm1
import math

def second_deriva(qe_input, qe_dyn, delta, r_order=2, modes=None, nscf=None, nscf_input=None):
    #
    # Check orthogonality 
    #
    if(not qe_dyn.check_orthogonality()):
        print(" ERROR ERROR ERROR!! ")
        print(" Use the dynamical matrix eigenvectors as input!! ")
        print(" Not the one normalized with the masses!! ")
        exit(1)

    qe_dyn.normalize_with_masses(masses)

    string="\n\n* * * 2nd derivaves calculations * * *\n\n"

    # DFTP results
    eig = np.array(qe_dyn.eig)

    if modes == None:
        modes = range(3, qe_dyn.nmodes) #skyp acustic modes at q=0

    for im in modes:
        print(" Calculating mode %d .... " % (im+1))
        if r_order == 1:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            folder="RIGHT_"+str(im)
            qe_right.write(scf_filename,folder)

        elif r_order == 2 or r_order == 3:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            qe_left =qe_dyn.generate_displacement(0, im, -delta)
            #
            folder="LEFT_"+str(im)
            qe_left.write(scf_filename,folder)
            #
            folder="RIGHT_"+str(im)
            qe_right.write(scf_filename,folder)

            if r_order == 3:
                der2_large=der2
                qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
                qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
                #
                folder="LEFT_bis_"+str(im)
                qe_left.write(scf_filename,folder)
                #
                folder="RIGHT_bis_"+str(im)
                qe_right.write(scf_filename,folder)
