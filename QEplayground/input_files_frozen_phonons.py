#!/usr/bin/python3
#
# 9/03/2020
#
# generates the Pwscf input files to compute the total energies of
# slightly moved configurations, later used to compute the phonon frequency

from QEplayground.pwscf import *
from QEplayground.matdyn import *
from QEplayground.units import autime2s, amu2au, thz2cm1
from QEplayground.pwout import *
import sys
from math import pi, sqrt

def input_files_frozen_phonons(qe_input, qe_dyn, delta, r_order=2, queue=6, nb_cores=12):
    #
    # Check orthogonality 
    #
    if(not qe_dyn.check_orthogonality()):
        print(" ERROR ERROR ERROR!! ")
        print(" Use the dynamical matrix eigenvectors as input!! ")
        print(" Not the one normalized with the masses!! ")
        sys.exit(1)

    masses=qe_input.get_masses()
    qe_dyn.normalize_with_masses(masses)

    scf_filename = "scf.in"

    modes = range(qe_dyn.nmodes)

    freqs = [qe_dyn.get_phonon_freq(0,im,unit="THz") for im in modes]
    freqs.sort()
    
    with open('qe.fp.job','w') as job :
        job.write("#/bin/csh -f\n\n")
        job.write("# Variables d'environnement SGE\n#$ -S /bin/csh\n#$ -cwd\n#$ -j y\n")
        job.write(f"#$ -pe qlogic {nb_cores}\n")
        job.write(f"#$ -q ib{queue}.q\n\n")
        #
        for im in modes :
            w_au = freqs[im]
            if w_au > 1e-6 :
                print(f"\n Calculating mode {im} .... ")
                print(f"with associated frequency     :    {w_au:.5f}  THz")
            else :
                continue
            if r_order == 1 :
                #
                qe_right = qe_dyn.generate_displacement(0,im,delta)
                folder = f"RIGHT_{im}"
                qe_right.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_{im}\n")
                print(f"File written in {folder}")
                
            elif r_order == 2 or r_order == 3:
                #
                qe_left =qe_dyn.generate_displacement(0, im, -delta)
                folder=f"LEFT_{im}"
                qe_left.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_l_{im}\n")
                print(f"File written in {folder}")
                #
                qe_right=qe_dyn.generate_displacement(0, im,  delta)
                folder=f"RIGHT_{im}"
                qe_right.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_r_{im}\n")
                print(f"File written in {folder}")
                      
                if r_order == 3:
                    #
                    qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
                    folder=f"LEFT_bis_{im}"
                    qe_left.write(scf_filename,folder)
                    job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_lb_{im}\n")
                    print(f"File written in {folder}")
                    #
                    qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
                    folder=f"RIGHT_bis_{im}"
                    qe_right.write(scf_filename,folder)
                    job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_lb_{im}\n")
                    print(f"File written in {folder}")
        #
        print(f"\n Job input file written as qe.fp.job \n")



def compute_derivative(r_order=1):
