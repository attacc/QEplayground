#!/usr/bin/python3
#
# Copyright (c) 2017-2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Calculate thermal lines on single phonon modes or 
# on linear combination of the different phonons
#

from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.pwout  import *
from QEplayground.units  import autime2s,amu2au,thz2cm1
from itertools import *

import math

def single_mode_thermal_line(qe_input, qe_dyn, modes):

    masses=qe_input.get_masses()
    ref_mass=max(masses)
    M       =1.0/(np.sum(np.reciprocal(masses)))
    M       =M*float(qe_input.system['nat'])

    string="\n\n* * * Thermal line * * *\n\n"
    string+="Reference  mass : %12.8f \n" % ref_mass
    string+="Reduced    mass : %12.8f [ref_mass units]  %12.8f  [amu]\n" % (M/ref_mass,M )

    print(string)

    # convert Mass to a.u.
    M       =M*amu2au

    # output reader

    scf_filename="scf.in"

    # DFTP results
    eig = np.array(qe_dyn.eig)

    qe_dyn.normalize_with_masses(masses)
    qe_dyn.normalize()
    qe_dyn.write_modes("mass_norm.modes")

    print("\n")

    for im in modes:

        w_atomic_units=eig[0,im]*(2.0*math.pi)/thz2cm1*autime2s*1e12

        delta=1.0/np.sqrt(2.0*w_atomic_units)

        print("Displacement mode %d = %12.8f a.u. \n" % (im+1,delta))
        qe_dyn.print_atomic_sigma_amplitude(0,im, delta)

        qe_dyn.print_atoms_sigma(0,im, delta)

#        qe_right=qe_dyn.generate_displacement(0, im,  delta)
#        qe_left =qe_dyn.generate_displacement(0, im, -delta)
        #
#        folder="LEFT_"+str(im)
#        qe_left.write(scf_filename,folder)
        #
#        folder="RIGHT_"+str(im)
#        qe_right.write(scf_filename,folder)



def print_thermal_lines(qe_input, qe_dyn):

    masses=qe_input.get_masses()
    ref_mass=max(masses)
    M       =1.0/(np.sum(np.reciprocal(masses)))
    M       =M*float(qe_input.system['nat'])

    string="\n\n* * * Thermal line * * *\n\n"
    string+="Reference  mass : %12.8f \n" % ref_mass
    string+="Reduced    mass : %12.8f [ref_mass units]  %12.8f  [amu]\n" % (M/ref_mass,M )

    print(string)

    # convert Mass to a.u.
    M       =M*amu2au

    scf_filename="scf.in"

    # DFTP results
    eig = np.array(qe_dyn.eig)

    qe_dyn.normalize_with_masses(masses)
    qe_dyn.normalize()
    qe_dyn.write_modes("mass_norm.modes")

    tl_list=qe_dyn.generate_thermal_lines(mode_range=[3])

    for tl in tl_list:
        print(tl)
        print("\n")
    #sign generation
#    single_mode_sign = [-1.0, 1.0]
#
#    nmodes=qe_dyn.nmodes-3 # remove acoustic modes
#    for mode_signs in product(single_mode_sign,repeat=nmodes):
#        qe_new=qe_dyn.generate_thermal_lines(mode_signs,mode_range=[3,4])
#        new_atoms=qe_new.get_atoms()
#        print("Termal line: "+str(mode_signs))
#        print(new_atoms)
##        print("\n")
