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
from QEplayground.pwscf import Pwscf
import fractions as frc
from QEplayground.lattice import rec_lat


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
        self.mode        = mode

        #Case of nondiagonal supercell
        if mode =='nondiagonal':
            self.Q = np.array(R)
            print('Nondiagonal supercell')
            if (qe_input.kpoints % self.Q[1] != 0).any():
                print(self.Q)
                print('ERROR: You must set a unit cell k-point mesh where Nx,Ny,Nz ')
                print('       are multiples of %d, %d, %d, respectively.' % (self.Q[1,0],self.Q[1,1],self.Q[1,2]))
                exit(0)
            self.R, self.new_latvec = self.find_nondiagonal()
        #Case of diagonal supercell
        else:
            self.new_latvec  = np.empty_like(self.latvec)
            self.R           = R
            for i in range(3):
                self.new_latvec[i]  = self.latvec[i]*R[i]

        self.sup_size    = int(self.R[0]*self.R[1]*self.R[2])
        self.new_atoms   = self.build_supercell()

    def find_integers(self,nums,g23,g12,g31,g123):
        """Compute integers for off-diagonal supercell matrix elements 
           Called by find_nondiagonal()
        """
        #Compute p (it's a modulo equation)
        if g23 == 1: p = 0
        else:
            for i in range(1,g23):
                if (nums[1]+i*nums[2]) % g23 == 0:
                    p=i
                    break
         #Compute q
        g12_r = int(g12/g123)
        g23_r = int(g23/g123)
        g31_r = int(g31/g123)
        if g12_r == 1: q = 0
        else:
            for i in range(1,g12_r):
                if (g23_r*nums[0]+i*g31_r*nums[1]) % g12_r == 0:
                    q=i
                    break
        #Compute r
        gg_r = int(g31*g23/g123)
        z = g23*nums[0]/g12+g31*q*nums[1]/g12
        if gg_r == 1: r = 0
        else:
            for i in range(1,gg_r):
                if (z+i*nums[2]) % gg_r == 0:
                    r=i
                    break
        return p,q,r

    def find_nondiagonal(self):
        """Nondiagonal supercell, based on [Phys. Rev. B 92, 184301]
        """
        Q = self.Q
        #Take care of components already at Gamma
        Q[1,np.where(Q[0]==0)]=1
        #Shift the q-point into the positive quadrant of the reciprocal unit cell
        Q[0,np.where(Q[0]<0)]+=Q[1,np.where(Q[0]<0)]
        #GCDs of Q[1] (in the logical order of the derivation)
        g23  = frc.gcd(Q[1,1],Q[1,2])
        g12  = frc.gcd(Q[1,0],Q[1,1])
        g31  = frc.gcd(Q[1,2],Q[1,0])
        g123 = frc.gcd(Q[1,0],frc.gcd(Q[1,1],Q[1,2]))
        #Integers needed to solve the supercell matrix equation    
        p,q,r = self.find_integers(Q[0],g23,g12,g31,g123)
        #Matrix elements (in order of derivation) and supercell matrix
        S_33 =        Q[1,2]
        S_22 =        Q[1,1]/g23
        S_23 =      p*Q[1,2]/g23
        S_11 =   g123*Q[1,0]/(g12*g31)
        S_12 = q*g123*Q[1,1]/(g12*g23)
        S_13 = r*g123*Q[1,2]/(g31*g23)
        self.S = np.array([[S_11,S_12,S_13],[0,S_22,S_23],[0,0,S_33]])
        #New lattice vectors and actual supercell size
        new_latvec = np.einsum('ij,jx->ix',self.S,self.latvec)
        R          = [self.S[0,0],self.S[1,1],self.S[2,2]]
        print(self.S)
        return R, new_latvec

    def build_supercell(self):
        #
        # Return new atomic position in bohr
        #
        from QEplayground.lattice import red2car, car2red
        from itertools import product
        latvec     = self.latvec
        R          = self.R
        atoms      = self.qe_input.get_atoms("bohr")
        new_atoms      = np.empty((self.sup_size*self.old_nat,3),dtype=float)
        for nz,ny,nx in product(range(int(self.R[2])),range(int(self.R[1])),range(int(self.R[0]))):
            cell=nx+ny*int(self.R[0])+nz*int(self.R[0]*self.R[1])
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

    def reciprocal(self,mode):
        """Function to compute reciprocal lattice
        """
        #Unit cell
        repvec = rec_lat(self.latvec)
        alat0=self.lattice_constants(self.latvec)[0]
        self.repvec = repvec/(2.0*np.pi/alat0)
        print('Old lattice vectors (alat):\n'+str(self.latvec/alat0))
        print('New lattice vectors (alat):\n'+str(self.new_latvec/alat0))
        print('Old reciprocal lattice vectors (2 pi/alat):\n'+str(self.repvec))
        #Supercell
        if mode=='diagonal' or mode =='keep_kpoints': self.new_repvec = np.array([self.repvec[i]/float(R[i]) for i in range(3)])
        else:
            self.S_inv_T = np.linalg.inv(self.S).T
            self.new_repvec = np.einsum('ij,jx->ix',self.S_inv_T,self.repvec)
        print('New reciprocal lattice vectors (2 pi/alat):\n'+str(self.new_repvec))

    def posint(self,value):
        return abs(int(round(value)))

    def write(self):
        R = self.R
        new_latvec = self.new_latvec
        alat = self.lattice_constants(new_latvec)
        qe = self.qe_input
        if self.mode=='diagonal':
            #A suggestion for a consistent new kpoint mesh 
            new_kpoints = [ceil(qe.kpoints[0]/self.R[0]), ceil(qe.kpoints[1]/self.R[1]), ceil(qe.kpoints[2]/self.R[2])]
        elif self.mode == 'keep_kpoints':
            new_kpoints = qe.kpoints
        else:
            #The compulsory new kpoint mesh - (sub)multiples of it are also fine but not consistent
            self.reciprocal('nondiagonal')
            new_kpoints = [1, 2, 3]
            new_kpoints[0] = self.posint(qe.kpoints[0])*np.linalg.norm(self.new_repvec[0])/np.linalg.norm(self.repvec[0])
            new_kpoints[1] = self.posint(qe.kpoints[1])*np.linalg.norm(self.new_repvec[1])/np.linalg.norm(self.repvec[1])
            new_kpoints[2] = self.posint(qe.kpoints[2])*np.linalg.norm(self.new_repvec[2])/np.linalg.norm(self.repvec[2])
            new_kpoints = [self.posint(new_kpoints[0]),self.posint(new_kpoints[1]),self.posint(new_kpoints[2])]


        qe_s = copy.deepcopy(qe)
        qe_s.system['nat'] = self.old_nat*self.sup_size
        qe_s.atoms = self.atoms_input(self.new_atoms)
        qe_s.control['prefix'] = qe.control['prefix'][:-1]+f"_{R[0]}x{R[1]}x{R[2]}'"
        qe_s.system['ibrav']=0
        if qe.system['celldm(1)'] != None:
            #
            # I keep the same alat if present
            #
            qe_s.system['celldm(1)']    = float(qe.system['celldm(1)'])
            qe_s.cell_units_output = 'alat'
        else:
            qe_s.cell_units_output = 'bohr'

        qe_s.system['celldm(2)']    = None
        qe_s.system['celldm(3)']    = None

        qe_s.atomic_pos_type   = 'bohr'
        qe_s.cell_parameters = new_latvec
        qe_s.convert_atoms(qe.atomic_pos_type)
        #Just a suggestion for the new bands
        if qe.system['nbnd'] != None: qe_s.system['nbnd'] = self.sup_size*int(qe.system['nbnd'])
        qe_s.kpoints = new_kpoints
        return(qe_s)
