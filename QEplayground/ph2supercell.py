#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
from itertools import product
from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.supercell import *

class map_phonons():
    """
    Map phonon in a supercell
    """
    def __init__(self,qe_input, qe_dyn, q_grid):
        self.qe_input=qe_input
        self.qe_dyn  =qe_dyn
        self.q_grid  =q_grid
        print(" \n\n\n * * * Map phonons in a supercell * * *\n")
        print(" WARNING! This code works only without symmetries!!! nosym=.true. and noinv=.true. \n")

    def build_mapping(self):
        n_qpoints=self.q_grid[0]*self.q_grid[1]*self.q_grid[2]
        print(" Q grid "+str(self.q_grid))
        print(" Number of q-points "+str(n_qpoints))

        if(self.qe_dyn.nqpoints != n_qpoints):
            print("\n ERROR! number of q-points in matdyn.modes different from the input q-grid!! \n")
            exit(0)

        #
        # Build the supercell
        #
        superc=supercell(self.qe_input,R=self.q_grid,mode='diagonal')
        self.qe_s=superc.write()
        self.qe_s.write(self.qe_input.filename+"_supercell")
        print(" Supercell scf file: "+self.qe_input.filename+"_supercell")
        self.qe_dyn_s=Matdyn(self.qe_s)
        #
        # Mapping the phonons
        #
        nmodes_old=self.qe_dyn.nmodes
        nmodes_new=self.qe_dyn_s.nmodes
        #
        self.qe_dyn_s.nqpoints=1
        self.qe_dyn_s.qpoints =np.zeros([1,3])
        self.qe_dyn_s.eig     =np.zeros([1,nmodes_new])
        self.qe_dyn_s.eiv     =np.zeros([1,nmodes_new,nmodes_new],dtype=complex)
        #
        # Copy old eivenvalues and eigenvectors (no phase yet)
        #
        for iq in range(n_qpoints):
            for im in range(nmodes_old):
                im_q=im+iq*nmodes_old
                self.qe_dyn_s.eig[0,im_q]=self.qe_dyn.eig[iq,im]
                for iq2 in range(n_qpoints):
                    self.qe_dyn_s.eiv[0,im_q,iq2*nmodes_old:(iq2+1)*nmodes_old]=self.qe_dyn.eiv[iq,im,:]

        # 
        # Add phases to the eigenvectors
        #
        new_atoms =self.qe_s.get_atoms(units="alat")
        new_natoms=int(self.qe_s.system["nat"])
        #
        for iq in range(n_qpoints):
            for im in range(nmodes_old):
                im_q=im+iq*nmodes_old
                for a in range(new_natoms):
                    sprod=np.dot(self.qe_dyn.qpoints[iq][:],new_atoms[a][:])
                    phase=np.exp(complex(0.0,1.0)*sprod*2.0*math.pi)
#                    self.qe_dyn_s.eiv[0,im_q,a*3:(a+1)*3]=self.qe_dyn_s.eiv[0,im_q,a*3:(a+1)*3]*phase

        self.qe_dyn_s.write_modes(filename="dynmat_supercell.out")




