#!/usr/bin/python3
# Copyright (C) 2015 Henrique Pereira Coutada Miranda, Alejandro Molina Sanchez
# All rights reserved.
#
# This file is part of yambopy
# (Slightly modified by C. Attaccalite)
#
import os
import re
import math
import numpy as np
import random
from QEplayground.auxiliary import *
from QEplayground.units import *
from itertools import product
import copy

#matdyn class

class MatdynIn():
    """
    Generate and manipulate quantum espresso input files for matdyn.x
    """
    def __init__(self):
        self.variable = dict()
        self.qpoints     = [] 

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        s = '&input\n'
        s += self.stringify_group('',self.variable) #print variable
        if len(self.qpoints) > 0:
          s+=("%d\n"%len(self.qpoints))
          for q in self.qpoints:
            s+=("%12.8lf %12.8lf %12.8lf")%tuple(q[:3])+"\n"
        return s

    def __setitem__(self,key,value):
        self.variable[key] = value

    def __getitem__(self,key):
        return self.variable[key]

    def stringify_group(self, keyword, group):
        if group != {}:
            string='\n'
            for keyword in group:
                string += "%20s = %s\n" % (keyword, group[keyword])
            string += "/\n"
            return string
        else:
            return ''


class Matdyn():
    """
    Class to read and plot the data from matdyn.modes files 
    """

    def __init__(self,qe_input, filename=None, cutoff_ph=100.0):
        """
        natoms is to be removed, but for now is left for legacy purposes
        """
        self.qe_input = qe_input
        self.filename = filename
        self.natoms     = int(qe_input.system['nat'])
        self.nmodes     = 3*int(self.natoms)
        self.cutoff_ph  = cutoff_ph # cutoff_ph on phonon energy, default 100 cm-1


        if(filename != None):
            self.read_modes(filename)

    def read_modes(self,filename):
        """
        read the modes
        """
        f = open(self.filename,'r')
        data_phon = f.readlines()
        f.close()

        #detect dimensions of the file
        #qpoints
        qpattern="q ="
        lines = iter(data_phon)
        nqpoints=sum(qpattern in line for line in lines)
        self.nqpoints= nqpoints
        #empty stuff
        eig = []
        eiv = []
        qpoints = []

        #read qpoints, modes and energies
        for j in range(nqpoints):
            frec, v_frec = [], []
            k=2 + j*(self.nmodes*(self.natoms+1)+5)
            qpoints.append( float_from_string(data_phon[k]) )
            for i in range(self.nmodes):
                k=4 + j*(self.nmodes*(self.natoms+1)+5) + i*(self.natoms+1)
                y = float_from_string(data_phon[k])
                v_mode = []
                for ii in range(1,self.natoms+1):
                    z      = float_from_string(data_phon[k+ii])
                    v_atom = [complex(z[0],z[1]),complex(z[2],z[3]),complex(z[4],z[5])]
                    v_mode.append(v_atom)
                v_frec.append(v_mode)
                frec.append(y[1])
            eig.append(frec)
            eiv.append(v_frec)

        #store info
        self.qpoints = np.array(qpoints)
        self.eig     = np.array(eig)
        self.eiv     = np.array(eiv).reshape(nqpoints,self.nmodes,self.nmodes)

    def write_modes(self,filename=None, wmode='w'):
        """
        save the phonon modes in a file
        """
        s = " matrix written with qepy\n\n"
        for nq in range(self.nqpoints):
            s += ("q =  "+"%12.6lf "*3+"\n")%tuple(self.qpoints[nq])
            s += "*"*81+"\n"
            for n,mode in enumerate(self.eiv[nq]):
                phfreqmev = self.get_phonon_freq(nq,n+1,unit='THz')
                phfreqcm1 = self.get_phonon_freq(nq,n+1,unit='cm-1')
                s += '    freq ( %4d) = %12.6lf [ThZ] = %12.6lf [cm-1]\n'%(n+1,phfreqmev,phfreqcm1)
                for a in range(self.natoms):
                    xr,yr,zr = mode[a*3:(a+1)*3].real
                    xi,yi,zi = mode[a*3:(a+1)*3].imag
                    s += ("( "+"%12.6lf "*6+')\n')%(xr,xi,yr,yi,zr,zi)
            s += "*"*81+"\n"

        if filename:
            f = open(filename, wmode)
            f.write(s)
            f.close()
        else:
            print(s)

    def rotate_phonons(self,eps=1e-5,debug=False):
        """
        Rotate the degenerate states to align them with the x,y axis
        Arguments:
        eps -> threshold for finding degeneracies
        """
        eig = self.eig
        eiv = self.eiv

        #the eigenvalues are probably sorted but just in case...
        eig, eiv = zip(*sorted(zip(eig,eiv), key=lambda x: x[0]))
        eig = np.array(eig)
        eiv = np.array(eiv)

        #iterate over qpoints
        for nq,(eigq,eivq) in enumerate(zip(eig,eiv)):
            #detect the degeneracies
            keys = [group.mean() for group in np.split(eigq, np.where(np.diff(eigq) > eps)[0]+1)]
            unique_frequencies = dict([(key,0) for key in keys])
            #count degeneracies
            for eig in eigq:
                for key in keys:
                    if np.isclose(key,eig,atol=eps):
                        unique_frequencies[key] += 1
            for val,deg in sorted(unique_frequencies.items()):
                if debug: print("frequency: %12.8lf [Thz] %12.8lf [cm-1] degeneracy: %2d" % (val,val,deg))
                if deg > 1:
                    #get the indexes of the modes with the same frequency
                    indexes = [n for n,f in enumerate(eigq) if np.isclose(f,val,atol=eps)]
                    r = np.array(eivq[indexes])
                    if debug:
                        print("input basis:")
                        for n in range(deg):
                            print("mode: %3d"%n)
                            for i in range(self.natoms):
                               print("atom %3d"%i+("%12.8lf"*3)%tuple(r[n,i*3:(i+1)*3].real))
                    #we make sure the first column vector the matrix r in non zero
                    rows,cols = r.shape
                    n = 0
                    while np.isclose(r[:,n],np.zeros([len(r)])).all():
                        n += 1
                    q, a = np.linalg.qr(r[:,n:])
                    r[:,n:] = a
                    if debug:
                        print("canonical basis:")
                        for n in range(deg):
                            print("mode: %3d"%n)
                            for i in range(self.natoms):
                                print("atom %3d"%i+("%12.8lf"*3)%tuple(r[n,i*3:(i+1)*3].real))
                    eivq[indexes] = r
                self.eiv[nq] = eivq

    def plot_eigen(self,path=[]):
        """ plot the phonon frequencies using matplotlib
        """
        import matplotlib.pyplot as plt

        if path:
            if isinstance(path,Path):
                path = path.get_indexes()
            plt.xticks( *zip(*path) )
        plt.ylabel('\\omega (cm$^{-1}$)')

        #plot vertical line
        for point in path:
            x, label = point
            plt.axvline(x)

        #plot bands
        eig = np.array(self.eig)
        for ib in range(self.nmodes):
           plt.plot(range(self.nqpoints),eig[:,ib], 'r-', lw=2,marker="o")
        plt.show()

    def get_phonon_freq(self,nq,n,unit="eV"):
        from QEplayground.units import cm1toeV,eV2ha,cm12Thz
        """
        Get the value of the phonon frequency
        nq -> q-point from where to get the frequency from
        n  -> mode of the phonon
        """
        if   unit == "eV":
            factor = cm1toeV
        elif unit == "Ha":
            factor = cm1toeV*eV2ha
        elif unit == "THz":
            factor = cm12Thz
        elif unit == "cm-1":
            factor = 1
        else:
            raise ValueError('Unit %s not known'%unit)

        return self.eig[nq][n-1]*factor

    def normalize(self):
        """
        Normalize the displacements u^n_{ai} according to:
        sum_ai ( u^n_{ai} )**2 = 1
        """

        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                self.eiv[nq,n] /= np.linalg.norm(self.eiv[nq,n])

    def normalize_with_masses(self,masses):
        """
        Normalize the displacements u^n_{ai} according to:
        sum_{ai} M_a u^n_{ai} u^m_{ai} = delta_{nm}

        u -> displacement
        n -> phonon mode
        a -> atom index
        i -> direction
        M -> mass
        """

        masses = np.array(masses)

        #divide by masses
        if self.check_orthogonality():
            for nq in range(self.nqpoints):
                for n in range(self.nmodes):
                    for a in range(self.natoms):
                        self.eiv[nq,n,a*3:(a+1)*3] *= 1.0/math.sqrt(masses[a])
        else:
            print("These eigenvectors are non-orthogonal, probably they are already scaled by the masses so I won't do it")

        #enforce delta_nm
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                s = 0
                for a in range(self.natoms):
                    e = self.eiv[nq,n,a*3:(a+1)*3]
                    #get normalization constant
                    s += masses[a]*np.vdot(e,e).real
                self.eiv[nq,n] *= 1.0/math.sqrt(s)

    def check_orthogonality(self,atol=1e-5):
        """
        Check if the eigenvectors are orthogonal
        """

        orth_check=True
        orth = np.zeros([self.nmodes,self.nmodes])
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                e1 = self.eiv[nq,n]
                for m in range(self.nmodes):
                    e2 = self.eiv[nq,m]
                    orth[n,m] = np.vdot(e1,e2).real
            orth_check=orth_check and np.allclose(orth,np.eye(self.nmodes),atol=atol)
        return orth_check

    def check_normalization(self,masses,atol=1e-5):
        """
        Check if the displacements are normalized according to:
        sum_{ai} M_a u^n_{ai} u^m_{ai} = delta_{nm}
        """

        masses = np.array(masses)
        ref_mass = max(masses)
        masses = masses/ref_mass

        #check normalization
        norm_check=True
        norm = np.zeros([self.nmodes])
        for nq in range(self.nqpoints):
            for n in range(self.nmodes):
                s = 0
                for a in range(self.natoms):
                    e = self.eiv[nq,n,a*3:(a+1)*3]
                    #get normalization constant
                    s += masses[a]*np.vdot(e,e).real
                norm[n] = s
            norm_check = norm_check and np.allclose(norm,np.ones(self.nmodes),atol=atol)
        return norm_check

    def __str__(self):
        s = ""
        for nq in range(self.nqpoints):
            for n,mode in enumerate(self.eiv[nq]):
                phfreqmev = self.get_phonon_freq(nq,n+1,unit='eV')*1000
                phfreqcm1 = self.get_phonon_freq(nq,n+1,unit='cm-1')
                s+= 'mode: %d freq: %8.2lf meV %8.2lf cm-1\n'%(n+1,phfreqmev,phfreqcm1)
                for a in range(self.natoms):
                    s += ("%12.8lf "*3)%tuple(mode[a*3:(a+1)*3].real)
                    s += ("  + i ("+"%12.8lf "*3+")\n")%tuple(mode[a*3:(a+1)*3].imag)
        return s




    def generate_displacement(self, iq, imode, delta):
        #
        # Displace atoms along a phonon mode [iq,imode]
        # the displacement 'delta' is in atomic units
        # The new atomic positions are written in bohr
        #
        qe_new = copy.deepcopy(self.qe_input)

        atoms      = self.qe_input.get_atoms("bohr")
        new_atoms  = np.empty((self.natoms,3),dtype=float)

        for a in range(self.natoms):
            e = self.eiv[iq,imode,a*3:(a+1)*3]
            new_atoms[a][:]=atoms[a][:]+e.real*delta

        qe_new.set_atoms(new_atoms,units='bohr')
        return qe_new

    def print_atoms_sigma(self, iq, imode, delta):
        #
        # Print atoms sigma for a given phonon mode
        # in a.u. 
        #
        masses = self.qe_input.get_masses()
        for a in range(self.natoms):
            e = self.eiv[iq,imode,a*3:(a+1)*3]
            norm  = np.sqrt(np.vdot(e,e).real)
            sigma = float(norm*delta/np.sqrt(masses[a]*amu2au))
            print("Atom %d  mass %12.8f sigma %12.8f" % (a,masses[a], sigma))
