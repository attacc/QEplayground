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
        # Check orthogonality 
        #
        if(not self.qe_dyn.check_orthogonality()):
            print(" ERROR ERROR ERROR!! ")
            print(" Use the dynamical matrix eigenvectors as input!! ")
            print(" Not the one normalized with the masses!! ")
            exit(1)
        #
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
        new_atoms =self.qe_s.get_atoms(units="bohr")
        new_natoms=int(self.qe_s.system["nat"])
        #
        tpiba=2.0*math.pi/float(self.qe_input.system['celldm(1)'])
        #
        phases=np.zeros([n_qpoints,new_natoms],dtype=float)

        eps=1e-5

        for iq in range(n_qpoints):
            for a in range(new_natoms):
                sprod=np.dot(self.qe_dyn.qpoints[iq][:],new_atoms[a][:]*tpiba)
                phases[iq,a]=np.real(np.exp(1j*sprod))
                print(" Phase [q= %d, a= %d ] = %lf " % (iq,a,phases[iq,a]))
                if iq !=0 and abs(phases[iq,a])<=eps:
                    print("Zero phase for atom %d at q= %iq ! Please check the code! ")
                    exit(0)

        for im in range(nmodes_old):
            for iq in range(n_qpoints):
                im_q=im+iq*nmodes_old
                for a in range(new_natoms):
                    self.qe_dyn_s.eiv[0,im_q,a*3:(a+1)*3]=self.qe_dyn_s.eiv[0,im_q,a*3:(a+1)*3]*phases[iq,a]

        #
        # Sorf phonons
        #
        sort_idx=np.argsort(self.qe_dyn_s.eig,axis=1)
        self.qe_dyn_s.eig=np.sort(self.qe_dyn_s.eig,axis=1)
        new_eig=np.empty_like(self.qe_dyn_s.eiv)
        for im in range(nmodes_new):
            new_eig[0,im,:]=self.qe_dyn_s.eiv[0,sort_idx[0][im],:]
        self.qe_dyn_s.eiv=new_eig
        #
        # Normalize eigevectors again
        #
        self.qe_dyn_s.normalize()
        #
        # Write output
        #
        self.qe_dyn_s.write_modes(filename="dynmat_supercell.out")
