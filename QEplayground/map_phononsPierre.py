#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
from itertools import product, groupby
from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.supercell import *
import sys
import numpy as np
import math

class map_phonons():
    """
    Map phonon in a supercell
    no_invar_ph = remove phonon modes invariant under inversion symmetry
                  modulo a reciprocal lattice vector
    """
    def __init__(self,qe_input, qe_dyn, ff, new_supercell_name, new_dynmat_name):           # ff stands for folding factor
        
        print(" \n\n\n * * * Map phonons in a supercell * * *\n")
        print(" This code works only without symmetries!!! \n")
        
        self.qe_input=qe_input

        if not qe_input.is_it_true(qe_input.system['noinv']):
            self.qe_input.system['noinv']='.true.'
            print(" WARNING! noinv flag set to .true. !!! \n")

        if not qe_input.is_it_true(qe_input.system['nosym']):
            self.qe_input.system['nosym']='.true.'
            print(" WARNING! nosym flag set to .true. !!! \n")

        self.qe_dyn  =qe_dyn
        self.ff      = ff
        self.new_supercell_name = new_supercell_name
        self.new_dynmat_name = new_dynmat_name
        superc=supercell(qe_input,R=ff, mode='keep_kpoints')
        self.qe_s=superc.write()
        self.qe_s.write(new_supercell_name)


    def get_translation_vectors(self):
        #
        # Works only for diagonal supercell a priori
        #
        translation_vectors=np.zeros((self.ff[0]*self.ff[1]*self.ff[2],3))
        latvec             =np.array(self.qe_input.cell_parameters)
        alat               =np.linalg.norm(self.qe_input.cell_parameters[0])
        #
        # Notice that I generate translation in the same order of the supercell.py
        # code
        #
        for nz,ny,nx in product(range(int(self.ff[2])),range(int(self.ff[1])),range(int(self.ff[0]))):
            cell=nx+ny*int(self.ff[0])+nz*int(self.ff[0]*self.ff[1])
            translation_vectors[cell,:]=(nx*latvec[0]+ny*latvec[1]+nz*latvec[2])/alat
            
        return(translation_vectors)

    def build_mapping(self,sort_ph=True,print_eig=False,norm_eig=True):
        #
        # Select all the q points to be folded
        #
        qpoints_all = self.qe_dyn.qpoints
        tr = self.get_translation_vectors()
        #
        print(" Translation vectors in alat ")
        for it in range(self.ff[2]*self.ff[1]*self.ff[0]):
            print(str(tr[it,:]))

        print("\n Q-points in 2 pi/alat ")
        for qpoint in qpoints_all:
            print(str(qpoint))
        
        #
        # Check orthogonality 
        #              
        if(not self.qe_dyn.check_orthogonality()):
            print(" ERROR ERROR ERROR!! ")
            print(" Use the dynamical matrix eigenvectors as input!! ")
            print(" Not the one normalized with the masses!! ")
            sys.exit(1)
        
        #
        # Build the supercell
        #
        superc=supercell(self.qe_input,R=self.ff,mode='keep_kpoints')
        self.qe_s=superc.write()
        self.qe_s.write(self.new_supercell_name)
        print("\nSupercell scf file: "+self.new_supercell_name+"\n")
        #
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
        # Copy old eivenvalues and eigenvectors (no phase yet) (phase is exp(1j*q.T) where q is 0)
        #
        #
        #
        n_qpoints=self.qe_dyn.nqpoints
        #
        # Print eigenvalues and eigenvectors
        # 
        if(print_eig):
            self.qe_dyn.write_modes()
        
        N = self.ff[0]*self.ff[1]*self.ff[2]
        
        for iq in range(N):
            for im in range(nmodes_old):
                im_q=im+iq*nmodes_old
                self.qe_dyn_s.eig[0,im_q]=self.qe_dyn.eig[iq,im]
                for iq2 in range(N):
                    self.qe_dyn_s.eiv[0,im_q,iq2*nmodes_old:(iq2+1)*nmodes_old]=self.qe_dyn.eiv[iq,im,:]
        # 
        #
        if(print_eig):
            self.qe_dyn_s.write_modes()
        #
        # Add phases to the eigenvectors
        #
        tpiba=2.0*math.pi/np.linalg.norm(self.qe_input.cell_parameters[0])
        #
        # I assume there is an unknown number of q-points
        #
        #
        # I order the q points in the same way than the translation vectors for simplicity reasons
        #
        qpoints_all_sorted = sorted(qpoints_all, key = lambda q: q[2])
        heights = [list(el) for _,el in groupby(qpoints_all, key = lambda q: q[2])]
        depths = []
        for i in range(len(heights)):
            depths.append([list(el) for _,el in groupby(heights[i], key = lambda q: q[1])])
        q_ordered = []
        for i in range(self.ff[2]):
            for j in range(self.ff[1]):
                for k in range(self.ff[0]):
                    q_ordered.append(depths[-i][-j][-k])
        #
        for iq in range(N):
            for im in range(nmodes_old):
                im_q=im+iq*nmodes_old
                for cell in range(N):
                    # q in units of 2pi/alat, Tr in units of alat
                    sprod=np.dot(tr[cell][:],q_ordered[cell][:]*2.0*math.pi) 
                    phase=np.exp(1j*sprod)
                    # Add phase
                    self.qe_dyn_s.eiv[0,im_q,cell*nmodes_old:(cell+1)*nmodes_old] *= phase
                    # Make it real
                    self.qe_dyn_s.eiv[0,im_q,cell*nmodes_old:(cell+1)*nmodes_old] = np.real(self.qe_dyn_s.eiv[0,im_q,cell*nmodes_old:(cell+1)*nmodes_old])

                    
        #
        if(sort_ph):
            #
            # Sort phonons
            #
            sort_idx=np.argsort(self.qe_dyn_s.eig,axis=1)
            self.qe_dyn_s.eig=np.sort(self.qe_dyn_s.eig,axis=1)
            new_eig=np.empty_like(self.qe_dyn_s.eiv)
            for im in range(nmodes_new):
                new_eig[0,im,:]=self.qe_dyn_s.eiv[0,sort_idx[0][im],:]
            self.qe_dyn_s.eiv=new_eig
        #
        # Normalize eigenvectors again
        #
        if(norm_eig):
            print("Normalize the new eigenvectors ")
            self.qe_dyn_s.normalize()
            print("Check normalization...")
            self.qe_dyn_s.check_orthogonality()
        else:
            for n in range(nmodes_new):
                print("New norm "+str(np.linalg.norm(self.qe_dyn_s.eiv[0,n])))
        #
        # Write output
        #
        self.qe_dyn_s.write_modes(filename=self.new_dynmat_name)
