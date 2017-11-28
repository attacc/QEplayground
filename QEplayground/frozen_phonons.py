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


r_order=3    # Richardson extrapolation order 
delta=0.02  # Displacement in a.u.

scf_filename    ="hBN.scf.in"
dynmat_filename ="dynmat.out"



qe_input=Pwscf(scf_filename)
qe_output=Pwout(qe_input)

# Setup pw.x executable
pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"
#pw="/home/elena/sources/qe-6.2/bin/pw.x"

#Serial job
#qe_input.set_run_options(pw=pw)

#Parallel job
qe_input.set_run_options(pw=pw, nprocs=2, npool=2)


# Pseudo-potential directory
qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"


print("\n\n * * * Frozen phonon calculation * * *\n")
print("Displacement    : %10.8f " % delta)
print("Richardson order: %d " % r_order)

# calculate the reduced mass
masses=qe_input.get_masses()
ref_mass=max(masses)
M       =1.0/(np.sum(np.reciprocal(masses)))
print("Reference  mass : %12.8f " % ref_mass)
print("Reduced    mass : %12.8f [ref_mass units]  %12.8f  [amu]" % (M/ref_mass,M ))

# convert Mass to a.u.
M       =M*float(qe_input.system['nat'])*2.0
M       =M*amu2au

folder="EQUIL"
qe_input.run(scf_filename,folder)
qe_output.read_output(scf_filename+".log", path=folder)
en_equil=qe_output.tot_energy
qe_dyn=Matdyn(qe_input,dynmat_filename)
eig = np.array(qe_dyn.eig)


for im in range(3,qe_dyn.nmodes):   #skyp acustic modes at q=0
    if r_order == 1:
        qe_right=qe_dyn.generate_displacement(0, im,  delta)
        folder="RIGHT_"+str(im)
        qe_right.run(scf_filename,folder)
        qe_output.read_output(scf_filename+".log", folder)
        en_right=qe_output.tot_energy

        der2=2.0*(en_right-en_equil)/delta**2

    elif r_order == 2 or r_order == 3:
        qe_right=qe_dyn.generate_displacement(0, im,  delta)
        qe_left =qe_dyn.generate_displacement(0, im, -delta)
        #
        folder="LEFT_"+str(im)
        qe_left.run(scf_filename,folder)
        qe_output.read_output(scf_filename+".log", folder)
        en_left=qe_output.tot_energy
        #
        folder="RIGHT_"+str(im)
        qe_right.run(scf_filename,folder)
        qe_output.read_output(scf_filename+".log", folder)
        en_right=qe_output.tot_energy

        der2=(en_right+en_left-2.0*en_equil)/delta**2

        if r_order == 3:
            der2_large=der2
            qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
            qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
            #
            folder="LEFT_bis_"+str(im)
            qe_left.run(scf_filename,folder)
            qe_output.read_output(scf_filename+".log", folder)
            en_left=qe_output.tot_energy
            #
            folder="RIGHT_bis_"+str(im)
            qe_right.run(scf_filename,folder)
            qe_output.read_output(scf_filename+".log", folder)
            en_right=qe_output.tot_energy

            der2_small=(en_right+en_left-2.0*en_equil)/(0.5*delta)**2

            der2=(4.0*der2_small-der2_large)/3.0


    autosi=1.0/autime2s
    omega=sqrt(der2/M)*autosi # We should understand this factor

    print("Mode %d   fr.(THz)   %12.8f   [DFTP %12.8f]" % (im+1,omega/(2.0*math.pi)/1e12,eig[0,im]/thz2cm1))
    print("Mode %d   fr.(cm^-1) %12.8f  [DFTP %12.8f] \n" % (im+1,omega/(2.0*math.pi)/1e12*thz2cm1,eig[0,im]))



