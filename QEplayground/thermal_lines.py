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
from QEplayground.units  import autime2s,amu2au,thz2cm1,au2kelvin
from QEplayground.utilities import bose
from itertools import *
import math
import random


def generate_thermal_lines(qe_dyn, T=0.0, folder="TL", n_tlines=None, tl2_lines=True, mode_range=None, debug=None):

    atoms      = qe_dyn.qe_input.get_atoms("bohr")
    new_atoms  = np.empty((qe_dyn.natoms,3),dtype=float)
    masses     = qe_dyn.qe_input.get_masses()
    #
    # Initialize random number
    #
    random.seed(a=100)
    #
    #
    # Check ortogonaly of the phonon eigenvectors
    # 
    if not qe_dyn.check_orthogonality():
        print("Error phonon eigenvectors not orthogonal!!! ")
        exit(0)
    #
    # Folder name with temperature
    #
    folder=folder+str(T)+"K"
    #
    if mode_range == None:
        mode_range=range(3, qe_dyn.nmodes) # Exclude the first 3 acustic modes

    tl_list=[]  # Thermal lines list

    if(len(mode_range)>10 and n_tlines == None):
        print(" WARNING! Too many modes! Number of thermal line set to 100! ")
        n_tlines=100

    if n_tlines == None:
        #
        # Generate all possible thermal lines
        # this part may explode in large/medium systems
        # lines_sign is a list of all possible sign combination 2^nmodes
        #
        single_mode_sign = [-1.0, 1.0]
        #
        all_sign=product(single_mode_sign,repeat=len(mode_range))

        for l_sign in all_sign:
            tl_list.append(l_sign)
    else:
        for i in range(n_tlines):
            tmp_list=[random.randint(0,1) for p in range(len(mode_range))]
            l_sign = [-1.0 if x == 0 else float(x) for x in tmp_list]
            tl_list.append(l_sign)
            if tl2_lines:
                minus_line=[-float(x) for x in l_sign]
                tl_list.append(minus_line)
                # Add the thermal line with opposite sign

    qe_new=copy.deepcopy(qe_dyn.qe_input)

    if debug:
        print("\n\n Sign list \n")
        for line_sign in tl_list:
            print(line_sign)
        print("\n\n Themal lines \n")

    ic=0
    for tl_line in tl_list:
        new_atoms  = atoms.copy()
        for im,im_sign in zip(mode_range,tl_line):
            #
            # Gaussian width
            #
            # see Eq. 12 of arXiv:1512.06377v1
            #
            w_au = qe_dyn.eig[0,im]*(2.0*math.pi)/thz2cm1*autime2s*1e12
            q_0  = 1.0/math.sqrt(2.0*w_au)
            q_T  = q_0*math.sqrt(1.0+2.0*bose(w_au,T/au2kelvin))
            #
            if debug:
                print("W and T atomic units : %14.10f, %14.10f " % (w_au,T/au2kelvin))
                print("Amplitude at T=0     : %14.10f " % q_0)
                print("Amplitude at finite T: %14.10f " % q_T)
            #
            delta =q_T*im_sign
            for a in range(qe_dyn.natoms):
                e = qe_dyn.eiv[0,im,a*3:(a+1)*3]
                new_atoms[a][:]=new_atoms[a][:]+e.real*delta/np.sqrt(masses[a]*amu2au)
        ic=ic+1

        qe_new.control['prefix']=qe_dyn.qe_input.control['prefix'].strip("'")+"_TL"+str(ic)
        qe_new.control['prefix']="'"+qe_new.control['prefix']+"'"

        if not debug:
            qe_new.set_atoms(new_atoms,units="bohr")
            qe_new.write(qe_dyn.qe_input.filename+"_TL"+str(ic),folder)
        else:
            print("Thermal line: %d " %(ic))
            print(new_atoms)



