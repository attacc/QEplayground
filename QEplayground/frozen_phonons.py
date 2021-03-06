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
from QEplayground.units  import autime2s,amu2au,thz2cm1,default_freq_thr
import math
import pathlib

def frozen_phonons(qe_input, qe_dyn, delta, r_order=2, modes=None, freq_thr=default_freq_thr):
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

    string="\n\n* * * Frozen phonon calculations * * *\n\n"
    string+="Min and max  mass : %12.8f  -  %12.8f  [amu]\n" % (min(masses),max(masses))

    print(string)
    ofile=open("frozen_phonons.log","w")
    ofile.write(string)

    # fictitious mass
    M       =amu2au #1.0

    # Rescale delta because masses are in amu
    delta=delta

    # output reader
    qe_output=Pwout(qe_input)

    scf_filename="scf.in"

    #Equilibrium calculation
    
    folder="EQUIL"
    if pathlib.Path(folder).exists():
        print("WARNING! EQUIL folder exists, skipping equilibrium calculation ")
    else:
        qe_input.run(scf_filename,folder)

    qe_output.read_output(scf_filename+".log", path=folder)
    en_equil=qe_output.tot_energy

    # DFTP results
    eig = np.array(qe_dyn.eig)

    if modes == None:
        modes = range(qe_dyn.nmodes)  # PIERRE : added a threshold on the frequency instead of skippping the first 3

    print("\nNumber of modes "+str(len(modes))+"\n")

    for im in modes:
        w_au = qe_dyn.get_phonon_freq(0,im+1,unit='Ha')
        if w_au < freq_thr:
            print(f" Mode {im+1} skipped too low frequency")
            continue
         
        print(" Calculating mode %d .... " % (im+1))
        w_au = qe_dyn.eig[0,im]*(2.0*math.pi)/thz2cm1*autime2s*1e12
        if w_au > 1e-6:
            q_0  = 1.0/math.sqrt(2.0*w_au)
        else:
            continue
        print("Mode Amplitude at T=0     : %14.10f  a.u." % q_0)
        if r_order == 1:
            qe_right=qe_dyn.generate_displacement(0, im,  delta)
            folder="RIGHT_"+str(im)
            #
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
        omega=math.sqrt(der2/M)*autosi/(2.0*math.pi)/1e12

        diff=abs(omega/(eig[0,im]/thz2cm1)*100.0-100.0)

        string ="Mode %d   fr.(THz)   %12.8f  [DFTP %12.8f]\n" % (im+1,omega,eig[0,im]/thz2cm1)
        string+="Mode %d   fr.(cm^-1) %12.8f  [DFTP %12.8f]  DIFF  %7.3f%% \n" % (im+1,omega*thz2cm1,eig[0,im],diff)

        print(string)
        ofile.write(string)

    ofile.close()

