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
import math


r_order=2    # Richardson extrapolation order 
delta=0.01  # Displacement in a.u.

scf_filename    ="diamond.scf.in"
dynmat_filename ="dynmat.out"

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"
#pw="/home/elena/sources/qe-6.2/bin/pw.x"

qe_input=Pwscf(scf_filename)
qe_output=Pwout(qe_input)

qe_input.set_run_options(pw=pw)
folder="EQUIL"
qe_input.run(scf_filename,folder)
qe_output.read_output(folder+"/"+scf_filename+".log")
en_equil=qe_output.tot_energy

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

        der2=(en_right+en_left-2.0*en_equil)/delta**2

        M=1./2.*12.0107*1.660539040e-27/9.109382e-31
        autosi=1.0/2.4188843265e-17
        omega=sqrt(der2/M)*autosi

        print("Mode %d   fr.(THz)   %12.8f " % (im,omega/(2.0*math.pi)/1e12))
        print("Mode %d   fr.(cm^-1) %12.8f \n" % (im,omega/(2.0*math.pi)/1e12*33.3564))



