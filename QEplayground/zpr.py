#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Calculate frozen phonon starting from pwscf input and dynmat output
#

from pwscf  import *
from matdyn import *
from pwout  import *
from units  import ha2ev
import math


r_order=3    # Richardson extrapolation order 
delta=0.01  # Displacement in a.u.

scf_filename    ="diamond.scf.in"
dynmat_filename ="dynmat.out"

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"
#pw="/home/elena/sources/qe-6.2/bin/pw.x"

qe_input=Pwscf(scf_filename)
qe_output=Pwout(qe_input)

# Setup the band gap calculation

qe_input.electrons['diago_full_acc']='.true.'
qe_input.system['nbnd']             =7

qe_input.set_run_options(pw=pw)
folder="EQUIL"
qe_input.run(scf_filename,folder)

qe_output.read_output(scf_filename+".log", path=folder)
dir_gap_eq, ind_gap_eq=qe_output.find_the_gap()

qe_dyn=Matdyn(qe_input,dynmat_filename)
masses=qe_input.get_masses()
qe_dyn.normalize()
qe_dyn.normalize_with_masses(masses)

print("\n\n * * * ZPR contribution for the different phonons * * *\n")
print("Displacement    : %10.8f " % delta)
print("Richardson order: %d " % r_order)

print("Equilibrium indirect gap: %12.8f " % (ind_gap_eq*ha2ev))
print("Equilibrium direct gap  : %12.8f " % (dir_gap_eq*ha2ev))

for im in range(3,6):

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

    print("Mode %d   der2 gap   %12.8f " % (im, der2*ha2ev))



