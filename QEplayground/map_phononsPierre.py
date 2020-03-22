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
        self.qe_input=qe_input
        self.qe_dyn  =qe_dyn
        self.ff      = ff
        self.new_supercell_name = new_supercell_name
        self.new_dynmat_name = new_dynmat_name
        superc=supercell(qe_input,R=ff, mode='keep_kpoints')
        self.qe_s=superc.write()
        self.qe_s.write(new_supercell_name)
        
        print(" \n\n\n * * * Map phonons in a supercell * * *\n")
        print(" WARNING! This code works only without symmetries!!! nosym=.true. and noinv=.true. \n")

    def get_translation_vectors(self):
        #
        # Works only for diagonal supercell a priori
        #
        qe_in = self.qe_input
        qe_in_s = self.qe_s
        #old_atoms = qe_in.get_atoms()
        old_natoms = int(qe_in.system["nat"])
        new_atoms = qe_in_s.get_atoms(units="alat")
        new_natoms = int(qe_in_s.system["nat"])
        translation_vectors = np.zeros((new_natoms), dtype = list)

        for a in range(new_natoms) :
            translation_vectors[a] = np.ones(3) + np.array(new_atoms[a]) - np.array(new_atoms[a%old_natoms])
        
        return(translation_vectors)

    def build_mapping(self):
        #
        # Select all the q points to be folded
        #
        qpoints_all = self.qe_dyn.qpoints
        folded_qpoints = [np.array([0.0, 0.0, 0.0])]
        index_folded_qpoints = [0]
        ffinv = np.reciprocal(self.ff, dtype = float)
        ffinv = np.around(ffinv,decimals=4)
        for ind,q in enumerate(qpoints_all) : 
            if np.any(np.greater_equal(q,ffinv)):
                folded_qpoints.append(q)
                index_folded_qpoints.append(ind)
        n_qpoints = len(folded_qpoints)
        #
        # Check orthogonality 
        #
        if(not self.qe_dyn.check_orthogonality()):
            print(" ERROR ERROR ERROR!! ")
            print(" Use the dynamical matrix eigenvectors as input!! ")
            print(" Not the one normalized with the masses!! ")
            sys.exit(1)
        #
        #
        # Build the supercell
        #
        #superc=supercell(self.qe_input,R=self.ff, mode='keep_kpoints')
        #self.qe_s=superc.write()
        #self.qe_s.write(self.new_supercell_name)
        print(" Supercell scf file: "+self.new_supercell_name)
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
        for iq,ifolded in enumerate(index_folded_qpoints):
            for im in range(nmodes_old):
                im_q=im+iq*nmodes_old
                self.qe_dyn_s.eig[0,im_q]=self.qe_dyn.eig[ifolded,im]
                for iq2 in range(n_qpoints):
                    self.qe_dyn_s.eiv[0,im_q,iq2*nmodes_old:(iq2+1)*nmodes_old]=self.qe_dyn.eiv[ifolded,im,:]

        # 
        # Add phases to the eigenvectors
        #
        new_atoms =self.qe_s.get_atoms(units="alat")
        new_natoms=int(self.qe_s.system["nat"])
        #
        tpiba=2.0*math.pi/float(self.qe_input.system['celldm(1)'])
        #
        phases=np.zeros([n_qpoints,new_natoms],dtype=float)

        eps=1e-5
        tr = self.get_translation_vectors()
        print("translation vectors :",tr)
        for iq,ifolded in enumerate(index_folded_qpoints): 
            for a in range(new_natoms):
                sprod=np.dot(self.qe_dyn.qpoints[ifolded][:],tr[a][:]*2.0*math.pi) # q in units of 2pi/alat, Tr in units of alat
                print("scalar product =",sprod)
                print("exponential value :",np.exp(1j*sprod))
                phases[iq,a]=np.real(np.exp(1j*sprod))
                print(f" Phase [q= {iq}, a= {a} ] = {phases[iq,a]}\n ")
                if iq !=0 and abs(phases[iq,a])<=eps:
                    print("Zero phase for atom %d at q= %iq ! Please check the code! ")
                    sys.exit(1)

        for im in range(nmodes_old):
            for iq in range(n_qpoints):
                im_q=im+iq*nmodes_old
                for a in range(new_natoms):
                    self.qe_dyn_s.eiv[0,im_q,a*3:(a+1)*3] *= phases[iq,a]

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
        self.qe_dyn_s.normalize()
        #
        # Write output
        #
        self.qe_dyn_s.write_modes(filename=self.new_dynmat_name)
