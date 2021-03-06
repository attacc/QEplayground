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
    """ qe_dyn is an instance of Matdyn() ????
"""
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
            #Pierre# doable with random.choice() in only one lign
            tmp_list=[random.randint(0,1) for p in range(len(mode_range))] 
            l_sign = [-1.0 if x == 0 else float(x) for x in tmp_list]
            #eg : l_signP=[random.choice(single_mode_sign) for p in range(len(mode_range))]
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
        print("\n\n Thermal lines \n")

    ic=0
    #
    # Normalize with masses
    #
    #Pierre# masses is already defined
    masses=qe_dyn.qe_input.get_masses()
    qe_dyn.normalize_with_masses(masses)
    #
    #
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
                new_atoms[a][:] += e.real*delta/math.sqrt(amu2au) #/np.sqrt(masses[a]*amu2au) 
        ic+=1  

        qe_new.control['prefix']=qe_dyn.qe_input.control['prefix'].strip("'")+"_TL"+str(ic)
        qe_new.control['prefix']="'"+qe_new.control['prefix']+"'"

        if not debug:
            qe_new.set_atoms(new_atoms,units="bohr")
            qe_new.write(qe_dyn.qe_input.filename+"_TL"+str(ic),folder)
        else:
            print("Thermal line: %d " %(ic))
            print(new_atoms)


def generate_ZG_conf(qe_dyn, T=0.0, folder="ZG", new_filename=None, excluded_freq = None, mode_range=None, debug=None):
    """ #Pierre# added a extra argument new_filename that specifies the name of the file to be written. In the case when the 
    input file is given with its path, you can't write the new file name as the input + suffix
    
        #Pierre# other argument : excluded_freq = number of frequencies that should be zero, depends on the ASR used in dynmat
    """
    atoms      = qe_dyn.qe_input.get_atoms("bohr")
    new_atoms  = np.empty((qe_dyn.natoms,3),dtype=float)
    masses     = qe_dyn.qe_input.get_masses()
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
        mode_range=range(excluded_freq, qe_dyn.nmodes) # Exclude the first frequencies that should be zero with 
                                                       # zero with ASR (3 translations + up to 3 rotations)
 
    qe_new=copy.deepcopy(qe_dyn.qe_input)

    ic=0
    #
    # Normalize with masses
    #
    masses=qe_dyn.qe_input.get_masses()
    qe_dyn.normalize_with_masses(masses)
    #
    #
    new_atoms  = atoms.copy()
    for im in mode_range:
       #
       # Gaussian width
       #
       # see Eq. 12 of arXiv:1512.06377v1
       #
       w_au = qe_dyn.get_phonon_freq(0,im+1,unit='Ha')
       q_0  = 1.0/math.sqrt(2.0*w_au)
       q_T  = q_0*math.sqrt(1.0+2.0*bose(w_au,T/au2kelvin))
       #
       if debug:
          print("W and T atomic units : %14.10f, %14.10f " % (w_au,T/au2kelvin))
          print("W in cm-1 %14.10f " % qe_dyn.get_phonon_freq(0,im+1,unit='cm-1'))
          print("Amplitude at T=0     : %14.10f " % q_0)
          print("Amplitude at finite T: %14.10f " % q_T)
       #
       if (im % 2 ) ==0:
           delta =q_T
       else:
           delta =-q_T

       for a in range(qe_dyn.natoms):
           e = qe_dyn.eiv[0,im,a*3:(a+1)*3]
           new_atoms[a][:]=new_atoms[a][:]+e.real*delta/math.sqrt(amu2au) #/np.sqrt(masses[a]*amu2au)
        
       qe_new.control['prefix']=qe_dyn.qe_input.control['prefix'].strip("\"")+"_ZG"
       qe_new.control['prefix']="'"+qe_new.control['prefix']+"'"

       if not debug:
          qe_new.set_atoms(new_atoms,units="bohr")
          qe_new.write(str(new_filename)+"_ZG",folder)
       else:
          print("ZG line: ")
          print(new_atoms)
