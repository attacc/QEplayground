#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
# Part of the code has been taken from the yambopy projec
# Copyright (c) 2015, Henrique Miranda
#

import sys
import numpy as np
import re
import subprocess


class Pwscf:

    def __init__(self, filename=None):
        self.control    = {}
        self.system     = {}
        self.electrons  = {}
        self.atypes      = {}
        self.atoms       = []
        self.cell_parameters = []

        self.ktype = "automatic"
        self.kpoints = [1,1,1]
        self.shiftk = [0,0,0]
        self.klist = []

        self.set_pwscf_default()

        self.atomicp_pattern=r'\s*ATOMIC_POSITIONS\s*\{?\s*(\w*)\s*\}?'
        self.cellp_pattern  =r'\s*CELL_PARAMETERS\s*\{?\s*(\w*)\s*\}?'

        if filename:
           self.read(filename)

        self.set_run_options()

    def set_run_options(self, nprocs=1, nthreads=1, npool=1, mpi='mpirun', pw='pw.x'):
        self._pw      =pw
        self._mpi     =mpi
        self._nprocs  =nprocs
        self._nthreads=nthreads
        self._npool   =npool

    def run(self,filename,folder='.'):
        """ this function is used to run this job locally
        """
        subprocess.call('mkdir -p %s'%folder,shell=True, cwd='./')
        self.write("%s/%s"%(folder,filename))
        pwjob="%s -npool %d  -inp %s > %s.log" % (self._pw,self._npool,filename,filename)
        if self._nprocs == 1:
            subprocess.call('OMP_NUM_THREADS=%d %s' % (self._nthreads,pwjob),shell=True, cwd=folder)
        else:
            subprocess.call('OMP_NUM_THREADS=%d mpirun -np %d %s' % (self._nthreads,self._nprocs,pwjob), shell=True, cwd=folder)

#    def __del__(self):
#        print("Destroy class PwscfIn")

    def set_pwscf_default(self):
        self.control['calculation']  = 'scf'
        self.control['restart_mode'] = None
        self.control['prefix']       = 'test'
        self.control['pseudo_dir']   = './'
        self.control['outdir']       = '\'./\''
        self.control['wf_collect']   = None
        self.control['tprnfor']      = None

        self.system['ibrav']        = None
        self.system['celldm(1)']    = None
        self.system['celldm(2)']    = None
        self.system['celldm(3)']    = None
        self.system['nat']          = 0
        self.system['ntyp']         = 0
        self.system['nbnd']         = None
        self.system['ecutwfc']      = float(30.0)
        self.system['force_symmorphic'] = None
        self.system['occupations']  = None

        self.electrons['mixing_mode']     = None
        self.electrons['mixing_beta']     = None
        self.electrons['conv_thr']        = float(1E-7)
        self.electrons['diago_full_acc']  = None

        self._kpoints           = [0,0,0,0,0,0]

        self.cell_units_output  = 'bohr'
        self.atomic_pos_type = 'bohr'

    def write(self,filename):
        f = open(filename,'w')
        f.write(str(self))
        f.close()


    def stringify_group(self, keyword, group):
        string='&%s\n' % keyword
        for keyword in group:
            if group[keyword] != None:
                string += "%20s = %s\n" % (keyword, group[keyword])
        string += "/&end\n"
        return string


    def __str__(self):
        string=''
        string += self.stringify_group("control",self.control)
        string += self.stringify_group("system",self.system)
        string += self.stringify_group("electrons",self.electrons)

        if int(self.system['ibrav']) == 0:
            string += self.write_cell_parameters()

        string += self.write_atomicspecies()
        string += self.write_atoms()
        string += self.write_kpoints()
        return string

    def read_atomicspecies(self):
        lines = iter(self.file_lines)
        for line in lines:
            if "ATOMIC_SPECIES" in line:
                for i in range(int(self.system["ntyp"])):
                    atype, mass, psp = next(lines).split()
                    self.atypes[atype] = [mass,psp]

    def write_atomicspecies(self):
        string = "ATOMIC_SPECIES\n"
        for atype in self.atypes:
            string += " %3s %8s %20s\n" % (atype, self.atypes[atype][0], self.atypes[atype][1])
        return string

    def read_atoms(self):
        lines = iter(self.file_lines)
        regex=re.compile(self.atomicp_pattern)
        for line in lines:
            match = regex.match(line)
            if match is not None:
                self.atomic_pos_type = match.group(1).lower()
                for i in range(int(self.system["nat"])):
                    atype, x,y,z = next(lines).split()
                    self.atoms.append([atype,[float(i) for i in [x,y,z]]])

    def write_atoms(self):
        string = "ATOMIC_POSITIONS { %s }\n"%self.atomic_pos_type
        for atom in self.atoms:
            string += "%3s %14.10lf %14.10lf %14.10lf\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2])
        return string

    def read_cell_parameters(self):
        from units import ang2au
        #
        # Internal cell-paramters always in atomic units
        # self.cell_units_output used only in output
        #
        self.cell_parameters = [[1,0,0],[0,1,0],[0,0,1]]
        ibrav = int(self.system['ibrav'])
        regex=re.compile(self.cellp_pattern)
        if ibrav == 0:
            lines = iter(self.file_lines)
            for line in lines:
                match = re.match(regex, line)
                if match is not None:
                    cell_units = match.group(1)
                    if cell_units == "alat":
                        scale=float(self.system['celldm(1)'])
                    elif cell_units == "bohr":
                        scale=1.0
                    elif cell_units == "angstrom":
                        scale=ang2au
                    self.cell_units_output=cell_units

                    for i in range(3):
                        self.cell_parameters[i] = [ float(x)*scale for x in next(lines).split() ]
        elif ibrav == 4:
            a = float(self.system['celldm(1)'])
            c = float(self.system['celldm(3)'])
            self.cell_parameters = [[   a,          0,  0],
                                    [-a/2,np.sqrt(3)/2*a,  0],
                                    [   0,          0,c*a]]
        elif ibrav == 2:
            a = float(self.system['celldm(1)'])
            self.cell_parameters = [[ -a/2,   0, a/2],
                                    [    0, a/2, a/2],
                                    [ -a/2, a/2,   0]]
        else:
            print('ibrav = %d not implemented'%ibrav)
            exit(1)

    def write_cell_parameters(self):
        string = "CELL_PARAMETERS { %s }\n"%self.cell_units_output
        if self.cell_units_output == "alat":
            scale=1.0/float(self.system['celldm(1)'])
        elif self.cell_units_output == "bohr":
            scale=1.0
        elif self.cell_units_output == "angstrom":
            scale=1.0/ang2au
        else:
            exit(1)
        for i in range(3):
            vec=np.asarray(self.cell_parameters[i])*scale
            string += ("%14.10lf "*3+"\n")%(vec[0],vec[1],vec[2])
        return string

    def read_namelist(self,group):
        for line in self.file_lines:
            if "="  not in line:
                continue
            key, value=line.split("=")
            for keyword in group:
                if key.strip().lower() == keyword:
                    group[keyword] = value.strip().strip(',')

    def read(self, filename):
        ifile = open(filename,'r')
        self.file_lines=ifile.readlines()

        self.read_namelist(self.control)
        self.read_namelist(self.system)
        self.read_namelist(self.electrons)

        self.read_cell_parameters()
        self.read_atomicspecies()
        self.read_atoms()
        self.read_kpoints()

        ifile.close()


    def read_kpoints(self):
        lines = iter(self.file_lines)
        #find K_POINTS keyword in file and read next line
        for line in lines:
            if "K_POINTS" in line:
                #chack if the type is automatic
                if "automatic" in line.lower():
                    self.ktype = "automatic"
                    vals = list(map(float, next(lines).split()))
                    self.kpoints, self.shiftk = vals[0:3], vals[3:6]
                #otherwise read a list
                elif "gamma" in line.lower():
                    self.ktype = "gamma"
                else:
                    #read number of kpoints
                    nkpoints = int(next(lines).split()[0])
                    self.klist = []
                    self.ktype = ""
                    try:
                        lines_list = list(lines)
                        for n in range(nkpoints):
                            vals = lines_list[n].split()[:4]
                            self.klist.append( map(float,vals) )
                    except IndexError:
                        print("wrong k-points list format")
                        exit(1)

    def write_kpoints(self):
        string = "K_POINTS { %s }\n"%self.ktype
        if self.ktype == "automatic":
            string += ("%3d"*6+"\n")%tuple(self.kpoints + self.shiftk)
        elif self.ktype == "crystal" or self.ktype == "tpiba" :
            string += "%d\n" % len(self.klist)
            for i in self.klist:
              string += ('%12.8lf '*4+'\n') % tuple(i)
        return string

    def set_atoms(self, new_atoms, units=None):
        if units != None:
            self.atomic_pos_type=units

        positions_input = new_atoms.tolist()
        elements_input  = [self.atoms[i][0] for i in range(int(self.system['nat']))]
        self.atoms = [[elements_input[i], positions_input[i]] for i in range(int(self.system['nat']))]

    def get_atoms(self, units=None):
        from units     import ang2au,au2ang
        from lattice   import red2car,car2red

        atoms= np.array([atom[1] for atom in self.atoms])   #atom[0] is the atomic symbol; atom[1] is the 3 coord list
        units = units if units is not None else self.atomic_pos_type   #self.atomic_pos_type = crystal in my case

        if units == self.atomic_pos_type:
            return atoms

        scale_in =1.0
        if self.atomic_pos_type == "angstrom":
            scale_in = ang2au
        elif self.atomic_pos_type == "alat":
            scale_in = float(self.system['celldm(1)'])
        elif self.atomic_pos_type == "crystal":
            atoms = red2car(atoms, np.array(self.cell_parameters))
 
        for i in range(int(self.system["nat"])):
            atoms[i]=atoms[i]*scale_in # transform in bohr
 
        scale_out=1.0
        if units == "alat":
            scale_out=1.0/float(self.system['celldm(1)'])
        elif units == "angstrom":
            scale_out=au2ang
        elif units == "crystal":
            atoms = car2red(atoms, np.array(self.cell_parameters))

        for i in range(int(self.system["nat"])):
            atoms[i]=atoms[i]*scale_out

        return atoms

    def get_masses(self):
        """ Get an array with the masses of all the atoms
        """
        masses = []
        for atom in self.atoms:
            atype = self.atypes[atom[0]]
            mass = float(atype[0])
            masses.append(mass)
        return np.asarray(masses)

    def convert_atoms(self, units):
        #
        # Convert units of the internal atoms
        # 
        new_atoms=self.get_atoms(units).tolist()
        for atom, self_atom in zip(new_atoms,self.atoms):
            self_atom[1]=atom
        self.atomic_pos_type = units

    def bring_atoms_in_the_cell(self):
        import math
        old_units=self.atomic_pos_type
        self.convert_atoms('crystal')
        for atom in self.atoms:
            atom[1] = map( lambda coord: coord - math.floor(coord), atom[1])
        self.convert_atoms(old_units)
