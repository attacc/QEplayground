#
# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy projec
#
# Modified by C. Attaccalite and E. Cannuccia
#
#
import os
import re
import math
import numpy as np
from pwscf import Pwscf

class supercell():
    """
    A class to generate custom supercells from a quantum espresso input file
    """
    def __init__(self,qe_input,R,mode='diagonal'):
        """
        qe_input: a Pwscf() instance of an input file
        R:
        if default:
            R is a list of the repetitions of the unit-cell in each direction
        else:
            R contains the fractional coordinates of the q-point to be folded at Gamma in a nondiagonal supercell like [[m1,m2,m3],[n1,n2,n3]]
        units: atomic positions in fractional(/angstrom/bohr)
        """
        self.qe_input = qe_input
        self.latvec   = np.array(qe_input.cell_parameters)
        self.basis    = int(qe_input.system['nat'])
        self.atoms    = qe_input.atoms
        self.b2a      = 0.529177
        #Case of nondiagonal supercell
        if mode!='diagonal':
            self.Q = np.array(R)
            print('Nondiagonal supercell')
            if (qe_input.kpoints % self.Q[1] != 0).any():
                print('ERROR: You must set a unit cell k-point mesh where%s\
       Nx,Ny,Nz are multiples of %d,%d,%d, respectively.'%('\n',self.Q[1,0],self.Q[1,1],self.Q[1,2])) 
                exit()
            self.R, self.new_latvec = self.find_nondiagonal()
        #Case of diagonal supercell    
        else: 
            self.R = R
            self.sup_size = R[0]*R[1]*R[2]
            self.new_latvec = np.array([self.latvec[i]*R[i] for i in range(3)])

        self.sup_size = self.R[0]*self.R[1]*self.R[2]
        new_atoms = self.build_supercell()        
        if write:
            #PwIn() object that can be printed, written to file, etc.
            self.qe = self.write(new_atoms,mode)

