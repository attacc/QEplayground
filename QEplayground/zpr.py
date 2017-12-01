#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Calculate zpr starting from pwscf input and dynmat output
#

from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.pwout  import *
from QEplayground.units  import ha2ev
import math


def zpr(qe_input, qe_dyn, delta, r_order=2, modes=None):
    
    string="\n\n* * * Zero Point Motion calculations * * *\n\n"
    string+="Reference  mass : %12.8f \n" % ref_mass
    string+="Reduced    mass : %12.8f [ref_mass units]  %12.8f  [amu]\n" % (M/ref_mass,M )

    ofile=open("zpr.log","w")
    ofile.write(string)

    # output reader
    qe_output=Pwout(qe_input)

    scf_filename="scf.in"

    #Equilibrium calculation   
    folder="EQUIL"
    qe_input.run(scf_filename,folder)
    qe_output.read_output(scf_filename+".log", path=folder)
    dir_gap_eq, ind_gap_eq=qe_output.find_the_gap()

    print("\n\n * * * ZPR contribution for the different phonons * * *\n")
    print("Displacement    : %10.8f " % delta)
    print("Richardson order: %d " % r_order)

    print("Equilibrium indirect gap: %12.8f " % (ind_gap_eq*ha2ev))
    print("Equilibrium direct gap  : %12.8f " % (dir_gap_eq*ha2ev))

    # DFTP results
    eig = np.array(qe_dyn.eig)

    if modes == None:
        modes = range(3, qe_dyn.nmodes) #skyp acustic modes at q=0

    for im in modes:
        print(" Calculating mode %d .... " % (im+1))
        if r_order == 1:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            folder="RIGHT_"+str(im)
            qe_right.run(scf_filename,folder)
            qe_output.read_output(scf_filename+".log", path=folder)
            dir_gap_right, ind_gap_right=qe_output.find_the_gap()
            der2=2.0*(dir_gap_right-dir_gap_eq)/delta**2

        elif r_order == 2 or r_order == 3:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            qe_left =qe_dyn.generate_displacement(0, im, -delta)
            #
            folder="LEFT_"+str(im)
            qe_left.run(scf_filename,folder)
            qe_output.read_output(scf_filename+".log", path=folder)
            dir_gap_left, ind_gap_left=qe_output.find_the_gap()
            #
            folder="RIGHT_"+str(im)
            qe_right.run(scf_filename,folder)
            qe_output.read_output(scf_filename+".log", path=folder)
            dir_gap_right, ind_gap_right=qe_output.find_the_gap()

            der2=(dir_gap_right+dir_gap_left-2.0*dir_gap_eq)/delta**2

            if r_order == 3:
                der2_large=der2
                qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
                qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
                #
                folder="LEFT_bis_"+str(im)
                qe_left.run(scf_filename,folder)
                qe_output.read_output(scf_filename+".log", path=folder)
                dir_gap_left, ind_gap_left=qe_output.find_the_gap()
                #
                folder="RIGHT_bis_"+str(im)
                qe_right.run(scf_filename,folder)
                qe_output.read_output(scf_filename+".log", path=folder)
                dir_gap_right, ind_gap_right=qe_output.find_the_gap()

                der2_small=(dir_gap_right+dir_gap_left-2.0*dir_gap_eq)/(0.5*delta)**2

                der2=(4.0*der2_small-der2_large)/3.0

        string = "Mode %d   der2 gap   %12.8f "  % (im, der2*ha2ev)
        print(string)

        ofile.write(string)

    ofile.close()
