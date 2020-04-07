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
        #sys.exit(1)

    masses=qe_input.get_masses()
    qe_dyn.normalize_with_masses(masses)
    
    # manual override of the non-orthogonality
    _natoms = int(qe_input.system['nat'])
    _nmodes = 3*_natoms
    for nq in range(qe_dyn.nqpoints):
                for n in range(_nmodes):
                    for a in range(_natoms):
                        qe_dyn.eiv[nq,n,a*3:(a+1)*3] *= 1.0/math.sqrt(masses[a])

    scf_filename = "scf.in"

    modes = range(qe_dyn.nmodes)
    list_modes=[]
    
    freqs = [qe_dyn.get_phonon_freq(0,im,unit="THz") for im in modes]
    freqs.sort()

    job_name = "qe.fp.job"
    with open(job_name,'w') as job :
        job.write("#/bin/csh -f\n\n")
        job.write("# Variables d'environnement SGE\n#$ -S /bin/csh\n#$ -cwd\n#$ -j y\n")
        job.write(f"#$ -pe qlogic {nb_cores}\n")
        job.write(f"#$ -q ib{queue}.q\n\n")
        #
        for im in modes :
            w_au = freqs[im]
            if w_au > 1e-6 :
                print(f"\n Calculating mode {im+1} .... ")
                print(f"with associated frequency     :    {w_au:.5f}  THz")
                list_modes.append(im)
            else :
                continue
            if r_order == 1 :
                #
                #
                qe_right = qe_dyn.generate_displacement(0,im,delta)
                folder = f"RIGHT_{im+1}"
                qe_right.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x -inp {Folder}/{scf_filename} > output_r_{im+1}\n")
                print(f"File written in {folder}")
                #
                #
            elif r_order == 2 or r_order == 3:
                #
                #
                qe_left =qe_dyn.generate_displacement(0, im, -delta)
                folder=f"LEFT_{im+1}"
                qe_left.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x -inp {folder}/{scf_filename} > output_l_{im+1}\n")
                print(f"File written in {folder}")
                #
                #
                qe_right=qe_dyn.generate_displacement(0, im,  delta)
                folder=f"RIGHT_{im+1}"
                qe_right.write(scf_filename, folder = folder)
                job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_r_{im+1}\n")
                print(f"File written in {folder}")
                #
                #   
                if r_order == 3:
                    #
                    #
                    qe_left =qe_dyn.generate_displacement(0, im, -delta/2.0)
                    folder=f"LEFT_bis_{im+1}"
                    qe_left.write(scf_filename,folder)
                    job.write(f"mpirun -np $NSLOTS pw.x -inp {folder}/{scf_filename} > output_lb_{im+1}\n")
                    print(f"File written in {folder}")
                    #
                    #
                    qe_right=qe_dyn.generate_displacement(0, im,  delta/2.0)
                    folder=f"RIGHT_bis_{im+1}"
                    qe_right.write(scf_filename,folder)
                    job.write(f"mpirun -np $NSLOTS pw.x < {folder}/{scf_filename} > output_rb_{im+1}\n")
                    print(f"File written in {folder}")
        #
        #
        job.write("touch finished_job")
        print(f"\n Job input file written as {job_name} \n")
    #
    return(list_modes, job_name)


def compute_frequency(qe_input_equil, output_equil_filename, qe_dyn, list_modes, delta, r_order=1):
    #
    scf_filename = "scf.in"
    qe_output_equil = Pwout(qe_input_equil)
    qe_output_equil.read_output(output_equil_filename,path='.')
    en_equil = qe_output_equil.tot_energy
    #
    dfpt_freqs = [qe_dyn.get_phonon_freq(0,im,unit="THz") for im in range(qe_dyn.nmodes)]
    dfpt_freqs.sort()
    new_freqs = []
    for im in list_modes :
        print(f"Computing frequency of mode {im+1} .... ")
        #
        if r_order == 1 :
            #
            #
            folder = f"RIGHT_{im+1}"
            output_filename = f"output_r_{im+1}"
            qe_input = Pwscf(folder+"/"+scf_filename)
            qe_output= Pwout(qe_input)
            qe_output.read_output(output_filename,".")
            en_right = qe_output.tot_energy
            #
            der2=2.0*(en_right-en_equil)/delta**2
            #
            #
        elif r_order == 2 or r_order == 3:
            #
            #
            l_folder = f"LEFT_{im+1}"
            l_output = f"output_l_{im+1}"
            qe_input = Pwscf(l_folder+"/"+scf_filename)
            qe_output= Pwout(qe_input)
            qe_output.read_output(l_output,".")
            en_left = qe_output.tot_energy
            #
            r_folder = f"RIGHT_{im+1}"
            r_output = f"output_r_{im+1}"
            qe_input = Pwscf(r_folder+"/"+scf_filename)
            qe_output= Pwout(qe_input)
            qe_output.read_output(r_output,".")
            en_right = qe_output.tot_energy
            #
            der2=(en_right+en_left-2.0*en_equil)/delta**2
            #
            #
            if r_order == 3:
                #
                #
                der2_large = der2
                #
                lb_folder = f"LEFT_bis_{im+1}"
                lb_output = f"output_lb_{im+1}"
                qe_input = Pwscf(lb_folder+"/"+scf_filename)
                qe_output= Pwout(qe_input)
                qe_output.read_output(lb_output,".")
                en_left = qe_output.tot_energy
                #
                rb_folder = f"RIGHT_bis_{im+1}"
                rb_output = f"output_rb_{im+1}"
                qe_input = Pwscf(rb_folder+"/"+scf_filename)
                qe_output= Pwout(qe_input)
                qe_output.read_output(rb_output,".")
                en_right = qe_output.tot_energy
                #
                der2_small = (en_right+en_left-2.0*en_equil)/(0.5*delta)**2
                der2 = (4.0*der2_small-der2_large)/3.0
        #
        #
        #
        
        omega = sqrt(der2/amu2au)/autime2s/(2*pi)/1e12   # in THz
        new_freqs.append(omega)
        diff = abs(omega/sorted(dfpt_freqs)[im]*100.0-100.0)
        print(f"Mode {im+1} :    frequency {omega:.8f} (THz)    [DFPT {dfpt_freqs[im]:.8f} (THz)]   DIFF : {diff} %\n")
        #
    return(new_freqs)

    
