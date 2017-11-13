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


r_order=2    # Richardson extrapolation order 
delta=0.01  # Displacement in a.u.

scf_filename    ="diamond.scf.in"
dynmat_filename ="dynmat.out"

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"

qe_input=Pwscf(scf_filename)
qe_output=Pwout(qe_input)

qe_input.set_run_options(pw=pw)

qe_dyn=Matdyn(qe_input,dynmat_filename)

for im in range(0,qe_dyn.nmodes):
    if r_order == 2:
        qe_right=qe_dyn.generate_displacement(0, im,  delta)
        qe_left =qe_dyn.generate_displacement(0, im, -delta)
        #
        folder="LEFT_"+str(im)
        qe_left.run(scf_filename,folder)
        qe_output.read_output(folder+"/"+scf_filename+".log")
        en_left=qe_output.tot_energy
        #
        folder="RIGHT_"+str(im)
        qe_right.run(scf_filename,folder)
        qe_output.read_output(folder+"/"+scf_filename+".log")
        en_right=qe_output.tot_energy

        der2=(en_right-en_left)*2.0/delta**2

        print("Mode %d   frequency %f12.8    \n" % (im,der2))
        print("Eright %12.8f        Eleft %12.8f   \n" % (en_right,en_left))

