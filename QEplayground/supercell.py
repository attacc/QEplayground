#!/usr/bin/python3
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
from math import *
import numpy as np
import copy
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
            R contains the fractional coordinates of the q-point to 
            be folded at Gamma in a nondiagonal supercell like [[m1,m2,m3],[n1,n2,n3]]
        """
        self.qe_input    = qe_input
        self.latvec      = np.array(qe_input.cell_parameters)
        self.old_nat     = int(qe_input.system['nat'])
        self.atoms       = qe_input.atoms
        self.R           = R
        self.new_latvec  = np.array([self.latvec[i]*R[i] for i in range(3)])
        self.sup_size    = self.R[0]*self.R[1]*self.R[2]
        self.new_atoms   = self.build_supercell()

    def build_supercell(self):
        #
        # Return new atomic position in bohr
        #
        from lattice import red2car, car2red
        from itertools import product
        latvec     = self.latvec
        R          = self.R
        atoms      = self.qe_input.get_atoms("bohr")
        new_atoms      = np.empty((self.sup_size*self.old_nat,3),dtype=float)
        for nz,ny,nx in product(range(self.R[2]),range(self.R[1]),range(self.R[0])):
            cell=nx+ny*self.R[0]+nz*self.R[0]*self.R[1]
            for b in range(self.old_nat):
                new_atoms[cell*self.old_nat+b]=atoms[b] +nx*latvec[0] +ny*latvec[1] +nz*latvec[2]
        return new_atoms

    def lattice_constants(self,vec):
        return [np.linalg.norm(vec[0]),np.linalg.norm(vec[1]),np.linalg.norm(vec[2])]

    def atoms_input(self, new_atoms):
        """ Put the atomic element labels in the right order
        """
        positions_input = new_atoms.tolist()
        elements_input  = [[self.qe_input.atoms[i][0] for i in range(self.old_nat)] for j in range(self.sup_size)]
        elements_input  = [ item for sublist in elements_input for item in sublist ]
        atoms_input     = [[elements_input[i], positions_input[i]] for i in range(self.sup_size*self.old_nat)]
        return atoms_input


    def write(self):
        new_latvec = self.new_latvec
        alat = self.lattice_constants(new_latvec)
        qe = self.qe_input
        #A suggestion for a consistent new kpoint mesh 
        new_kpoints = [ceil(qe.kpoints[0]/self.R[0]), ceil(qe.kpoints[1]/self.R[1]), ceil(qe.kpoints[2]/self.R[2])]
        qe_s = copy.deepcopy(qe)
        qe_s.system['nat'] = self.old_nat*self.sup_size
        qe_s.atoms = self.atoms_input(self.new_atoms)
        qe_s.control['prefix'] = qe.control['prefix'][:-1]+"_s'"
        qe_s.system['ibrav']=0
        qe_s.atomic_pos_type = 'bohr'
        qe_s.cell_units      = 'bohr'
        qe_s.convert_atoms(qe.atomic_pos_type)
        qe_s.cell_parameters = new_latvec
        #Just a suggestion for the new bands
        if qe.system['nbnd'] != None: qe_s.system['nbnd'] = self.sup_size*int(qe.system['nbnd'])
        qe_s.kpoints = new_kpoints
        return qe_s
