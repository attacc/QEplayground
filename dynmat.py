# Copyright (C) 2017 Elena Cannuccia, Claudio Attaccalite
# All rights reserved.
#
import numpy as np
import re, sys
import math


class Dynmat:

    def __init__(self):
        #self._pol_vector       = np.zeros((int(nat),3,nmodes,nqpoints),dtype=complex)
        #self._pol_vector       = np.zeros((2,3),dtype=complex)
        #self._norma_polv = 0.0
        self._qvec             = np.zeros(3)
        #self._PHfreq       = 0.0
        self._qpattern    = 'q ='
        self._freq_pattern    = 'freq ('
        self._omega      = 0.1

    def __del__(self):
        print("Destroy class Dynmat")

    def read_dyn_output(self,ifile,nat):
        #
        nmodes=3*int(nat)
        nqpoints=0;  qpoints = []
        #
        PHfreq=[]
        lines=ifile.readlines()
        for line in lines:
            if self._qpattern in line:
                nqpoints += 1
            if self._freq_pattern in line:
               match = re.search(r'\[THz\] =\s*(\d*\.\d*)',line)
               freq = match.group(1) 
               PHfreq.append(float(freq)) 
        if len(PHfreq)==nmodes*nqpoints:
           print('*** #freq found is EQUIVALENT to the #modes ***')
        else:
           print('Error: #freq is DIFFERENT from #Modes')
             
        for iq in range(nqpoints):
            q_line = 2 + iq*(4 + nmodes*(int(nat)+1))
            qpoints.append(lines[q_line])
            self._pol_vector = np.zeros((nqpoints,nmodes,int(nat),3),dtype=complex)
            self._norma_polv = np.zeros((nqpoints,nmodes),dtype=float)
            mode_line = 5  
            for im in range(nmodes):
                for ia in range(int(nat)):
                    polv = np.zeros(int(nat)*3,dtype=float)
                    mode_line = mode_line + ia
                    polv = re.findall(r'-?0.\d*', lines[mode_line])
                    #
                    #  The phonon eigenvectors are in Cartesian axis, i.e. along x, y, z
                    #
                    self._pol_vector[iq][im][ia].real  = polv[0::2]
                    self._norma_polv[iq][im] = self._norma_polv[iq][im] + np.sum(np.square(self._pol_vector[iq][im][ia].real))
                    #
                    self._pol_vector[iq][im][ia].imag  = polv[1::2]
                    #print(ia, self._pol_vector[iq][im][ia].real)
                self._norma_polv[iq][im]=np.around(self._norma_polv[iq][im])
                #print(self._norma_polv[iq][im])
                mode_line = mode_line + 2
        return nqpoints, nmodes, PHfreq

    def set_displacement(self, iq, im, nat, it):
        self._displacement = np.zeros((int(nat),3),dtype=float)
        delta = 0.007   # Let's assume here delta in **Bohr** units because I will always set cell_param in **Bohr**
        for ia in range(int(nat)):
            self._displacement[ia,:] = np.around(self._pol_vector[iq][im][ia].real/self._norma_polv[iq][im]*delta*it, decimals=3)

    def get_new_atomic_positions(self, iq, im, nat, start_pos_bohr, it, units="bohr"):
        new_pos_bohr = [[],np.zeros((int(nat),3),dtype=float)]
        for ia in range(int(nat)):
            new_pos_bohr[ia] = [start_pos_bohr[ia][0],(start_pos_bohr[ia][1]+ self._displacement[ia]).tolist()]
        return new_pos_bohr

    def set_linear_combination(self, iq, modes, p_o_m, nat, start_pos_bohr, it, units="bohr"):
        new_pos_bohr = [[],np.zeros((int(nat),3),dtype=float)]
        ThermalLine = np.zeros((len(modes),int(nat),3),dtype=float)
        delta = 0.007   # Let's assume here delta in **Bohr** units because I will always set cell_param in **Bohr**
        for ia in range(int(nat)):
            for im in range(len(modes)):
                ThermalLine[ia,:] = ThermalLine[ia,:] + p_o_m[im]*np.around(self._pol_vector[iq][modes[im]][ia].real*delta*it, decimals=5)
                new_pos_bohr[ia] = [start_pos_bohr[ia][0],(start_pos_bohr[ia][1] +ThermalLine[ia,:]).tolist()]
        return new_pos_bohr
