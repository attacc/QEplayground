#!/usr/bin/python3

import sys
import numpy as np
import re
from utilities import   decimal,real_sn,r_or_d,AngToAu


class Pwscf:

    def __init__(self, filename=None):
        self.control    = {}
        self.system     = {}
        self.electrons  = {}
        self._kpoints    = []
        self._atoms      = {}
        self._forces     = {}
        self._atoms_positions = []
        self._cell_parameters = []
        self.set_pwscf_default()
        if filename != None:
            self.read(filename)


    def __del__(self):
        print("Destroy class PwscfIn")

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
        self._kpoints           = [0,0,0,0,0,0]

    def _read_word(self, value):
        value=value.strip().strip(',')
        return value

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
## Control section
        string += self.stringify_group("control",self.control)
        string += self.stringify_group("system",self.system)
        string += self.stringify_group("electrons",self.electrons)

## CELL_PARAMETERS section
        if len(self._cell_parameters)!=0:
           string += "CELL_PARAMETERS {bohr}\n"
#           for vec in self._cell_parameters:
#               self._write_word(string, None, vec[0:3])

## ATOMIC SPECIES section
        string += "ATOMIC_SPECIES\n"
#        for word in self._atoms:
#            self._write_word(string, word, self._atoms[word])

## ATOMIC_POSITIONS section
        string += "ATOMIC_POSITIONS {bohr}\n"
        print("Warning: ATOMIC_POSITIONS are going to be written in bohr")
        for atom in self._atoms_positions:
            string += str(atom[0])+" "+str(atom[1][0])+" "+str(atom[1][1])+" "+str(atom[1][2])+"\n"

## KPOINTS section
        string += "K_POINTS automatic\n"
        for ik in self._kpoints:
            string += str(ik)+""
##
        return string

    def _write_line(self,string, values):
        for value in values:
            string += "5%"+str(values)+"%5s"
        string += "\n"



    def read(self, filename):
        ifile = open(filename,'r')
        cellp_pattern  =r'\s*CELL_PARAMETERS\s*\{?\s*(\w*)\s*\}?'
        lines=ifile.readlines()
        for i, line in enumerate(lines):
            if re.search(cellp_pattern, line):
                match = re.search(cellp_pattern, line)
                cell_units = match.group(1)

                if cell_units == "bohr":
                    scale=1.0
                elif cell_units == "angstrom":
                    scale=AngToAu
                elif cell_units == "alat":
                    scale=float(self.system['celldm(1)'])
                
                for line in lines[i+1:i+4]:
                    cell_vectors = line.split()
                    vec=np.asarray(cell_vectors[0:3], dtype=np.float)*scale
                    self._cell_parameters.append(vec.tolist())

            if "ATOMIC_SPECIES" in line:
                for line in lines[i+1:i+1+int(self.system['ntyp'])]:
                    info_atoms = line.split()
                    self._atoms[info_atoms[0]] = info_atoms[1:]

            atomicp_pattern=r'\s*ATOMIC_POSITIONS\s*\{?\s*(\w*)\s*\}?'
            if re.search(atomicp_pattern, line):
                match = re.search(atomicp_pattern, line)
                pos_units = match.group(1)
                self._atoms_positions = [[],np.zeros((int(self.system['nat']),3),dtype=float)]
                
                for ia in range(int(self.system['nat'])):
                    info_atoms_pos = lines[i+1+ia].split()
                    self._atoms_positions[ia] = [info_atoms_pos[0],info_atoms_pos[1:4]]

                self._convert_to_bohr(pos_units)

            if "K_POINTS" in line:
                self._kpoints = lines[i+1]
            if "="  not in line:
                continue
            flag, value=line.split("=")
            for word in self.control:
                if flag.strip().lower() == word:
                    self.control[word] = self._read_word(value)
            for word in self.electrons:
                if flag.strip().lower() == word:
                    self.electrons[word] = self._read_word(value)
            for word in self.system:
                if flag.strip().lower() == word:
                    self.system[word] = self._read_word(value)
        ifile.close()
#
    def read_qe_output(self,ofile):
        lines=ofile.read()
        pattern_nelec=r'number'+'\s'+'of'+'\s'+'electrons'+'\s+'+'='+'\s+(\d\.\d*)'
        try:
            match = re.search(pattern_nelec, lines)
            nelec = match.group(1)
        except:
            print("Error in matching number of electrons")
        nbnd = int(0.5*float(nelec))
        ofile.seek(0)
        lines=ofile.readlines()
        pattern_egap='k = 0.0000 0.0000 0.0000'
        for i, line in enumerate(lines):
            if pattern_egap in line:
              info_gap = lines[i+2].split()
              top_val = float(info_gap[nbnd-1])
              bot_cond = float(info_gap[nbnd])
              egap = (bot_cond-top_val)/27.2114
        return nelec, nbnd, egap
        
        #for line in lines:
        #    if "Forces" in line:
        #        for line in lines[i+2:i+2+int(self.system['ntyp'])]:
        #            info_force_atoms = line.split()
        #            print(info_force_atoms[1], info_force_atoms[6:])
        #            self._forces[info_force_atoms[1]] = info_force_atoms[6:]
        #    if "Total force =" in line:
        #            info_total_force = line.split()
        #            Total_Force = info_total_force[3]
        #            print("Total Force = "+str(Total_Force))
        #    if "!    total" in line:
        #            info_total_energy = line.split()
        #            Total_Energy = info_total_energy[4]
        #            print("Total Energy = "+str(Total_Energy))

    def _convert_to_bohr(self, units):
        if units == "bohr":
            return
        elif units == "angstrom":
            self._rescale_pos(AngToAu)
        elif units == "alat":
            self._rescale_pos(float(self.system['celldm(1)']))
        elif units == "crystal":
            vec = np.zeros((3,3),dtype=float);    i=0
            atom_pos=         np.zeros((int(self.system['nat']),3), dtype=float)
            atom_pos_bohr=[[],np.zeros((int(self.system['nat']),3), dtype=float)]
            for cp in self._cell_parameters:
                vec[i,:]= np.asarray(cp[:], dtype=np.float) #in **Bohr**
                i=i+1
            for ia in range(int(self.system['nat'])):
                atom_pos[ia,:] = np.asarray(self._atoms_positions[ia][1:4], dtype=np.float)
                atom_pos_bohr[ia] = [self._atoms_positions[ia][0],np.dot(vec,atom_pos[ia]).tolist()]
                self._atoms_positions[ia]=atom_pos_bohr[ia]

    def get_cell_parameters(self):
        cell = np.zeros((3,3),dtype=float)
        for cp in self._cell_parameters:
                cell[i,:]= np.asarray(cp[:], dtype=np.float)

    def _rescale_pos(self, scale):
        for ia in range(int(self.system['nat'])):
            self._atoms_positions[ia][1:4]=scale*np.asarray(self._atoms_positions[ia][1:4], dtype=np.float)

    def get_starting_positions(self):
        start_atom_pos= [[],np.zeros((int(self.system['nat']),3), dtype=float)]
        for ia in range(int(self.system['nat'])):
            start_atom_pos[ia] = self._atoms_positions[ia]
        return start_atom_pos

    def set_positions(self,positions):
        if isinstance(positions, np.ndarray):
            for ia in range(int(self.system['nat'])):
                self._atoms_positions[ia][1:4]=positions[ia]
        else:
            print('passo qui')
            self._atoms_positions = positions


# Recently assigned and default values can be shown thanks to the 
# following function
#
    def show_values(self):
        for word in self.control:
            print (word+"="+str(self.control[word]))
        for word in self.system:
            print (word+"="+str(self.system[word]))
        for word in self.electrons:
            print (word+"="+str(self.electrons[word]))


