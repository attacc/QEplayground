#!/usr/bin/python3

import sys
import numpy as np
import re

class PwscfIn:

#
#   Useful regular expression for QE
#
    decimal=r'-?\d+\.\d*|-?\.\d+'
    real_sn=r'[+\-]?\d*\.\d*[E][+\-]\d\d?'
    r_or_d =r'('+real_sn + '|' + decimal + ')'

#
#   Units
# 
    AngToAu=1.8897261


    def __init__(self):
        self._control    = {}
        self._system     = {}
        self._electrons  = {}
        self._kpoints    = []
        self._atoms      = {}
        self._forces     = {}
        self._atoms_positions = []
        self._cell_parameters = []
        self.set_pwscf_default()

    def __del__(self):
        print("Destroy class PwscfIn")

    def set_pwscf_default(self):
        self._control['calculation']  = 'scf'
        #self._control['restart_mode'] = 'from_scratch'
        self._control['prefix']       = 'test'
        self._control['pseudo_dir']   = './'
        #self._control['outdir']       = './'
        self._control['wf_collect']   = True
        #self._control['tprnfor']      = True

        self._system['ibrav']        = None
        self._system['celldm(1)']    = None
        self._system['celldm(2)']    = None
        self._system['celldm(3)']    = None
        self._system['nat']          = 0
        self._system['ntyp']         = 0
        self._system['nbnd']         = None
        self._system['ecutwfc']      = float(30.0)
        #self._system['force_symmorphic'] = False
        self._system['occupations']  = None

        self._electrons['mixing_mode']     = None
        self._electrons['mixing_beta']     = None
        self._electrons['conv_thr']        = float(1E-7)
        self._kpoints           = [0,0,0,0,0,0]

    def _write_word(self, fout, word, value):
        if word == None:
            if isinstance(value, list) and len(value)==1:
                fout.write("    "+str(value[0])+"\n")
            if isinstance(value, list) and len(value)==2:
                fout.write("    "+str(value[0])+"  "+str(value[1])+"\n")
            elif isinstance(value, list) and len(value)==3:
                fout.write("    "+str(value[0])+"  "+str(value[1])+"  "+str(value[2])+"\n")
            return

        if value != None:
            if value == False:
                fout.write("    "+word+" "+"="+" "+".false."+",\n")
            elif value == True:
                fout.write("    "+word+" "+"="+" "+".true."+",\n")
            elif isinstance(value, str):
                 if re.search(r'\d', value):
                     fout.write("    "+word+" "+"="+" "+str(value)+",\n")
                 else:
                     fout.write("    "+word+" = '"+str(value)+"',\n")
            elif isinstance(value, list) and len(value)==2:
                fout.write("    "+word+"  "+str(value[0])+"  "+str(value[1])+"\n")
            elif isinstance(value, list) and len(value)==3:
                fout.write("    "+word+"  "+str(value[0])+"  "+str(value[1])+"  "+str(value[2])+"\n")
            else:
                fout.write("    "+word+" = "+str(value)+",\n")

    def _read_word(self, value):
        value=value.strip().strip(',').strip("'")
        if value == ".false.":
            return False
        elif value == ".true.":
            return True
        else:
            return value

    def write_pwscf(self, filename):

        # This is because the file() builtin function is removed in Python 3
        # 
        #if isinstance(filename,file):
        #    print(" passo qui")
        #    fout=filename
        #else:
        try:
            fout=open(filename,"w")
        except:
            print("Error opening %s " % filename)
            sys.exit(1)

## Control section
        fout.write("&control\n")
        for word in sorted(self._control):
            self._write_word(fout, word, self._control[word])
        fout.write("/ \n")
##
## System section
        fout.write("&system\n")
        for word in sorted(self._system):
            self._write_word(fout, word, self._system[word])
        fout.write("/ \n")
##
## Electrons section
        fout.write("&electrons\n")
        for word in self._electrons:
            self._write_word(fout, word, self._electrons[word])
        fout.write("/ \n")
##
## CELL_PARAMETERS section
        if len(self._cell_parameters)!=0:
           fout.write("CELL_PARAMETERS {bohr}\n")
           for vec in self._cell_parameters:
               self._write_word(fout, None, vec[0:3])
##
## ATOMIC SPECIES section
        fout.write("ATOMIC_SPECIES\n")
        for word in self._atoms:
            self._write_word(fout, word, self._atoms[word])
##
## ATOMIC_POSITIONS section
        fout.write("ATOMIC_POSITIONS {bohr}\n")
        print("Warning: ATOMIC_POSITIONS are going to be written in bohr")
        for atom in self._atoms_positions:
            fout.write(str(atom[0])+" "+str(atom[1][0])+" "+str(atom[1][1])+" "+str(atom[1][2])+"\n")
##
## KPOINTS section
        fout.write("K_POINTS automatic\n")
        for ik in self._kpoints:
            fout.write(str(ik)+"")
##
        fout.close()
#
# Values coming from QE input file are assigned to the corresponding
# word in the dictionary. The words in the dictionary are then 
# eventually re-defined.
#

    def read_qe_input(self, ifile):
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
                    scale=float(self._system['celldm(1)'])
                
                for line in lines[i+1:i+4]:
                    cell_vectors = line.split()
                    vec=np.asarray(cell_vectors[0:3], dtype=np.float)*scale
                    self._cell_parameters.append(vec.tolist())

            if "ATOMIC_SPECIES" in line:
                for line in lines[i+1:i+1+int(self._system['ntyp'])]:
                    info_atoms = line.split()
                    self._atoms[info_atoms[0]] = info_atoms[1:]

            atomicp_pattern=r'\s*ATOMIC_POSITIONS\s*\{?\s*(\w*)\s*\}?'
            if re.search(atomicp_pattern, line):
                match = re.search(atomicp_pattern, line)
                pos_units = match.group(1)
                self._atoms_positions = [[],np.zeros((int(self._system['nat']),3),dtype=float)]
                
                for ia in range(int(self._system['nat'])):
                    info_atoms_pos = lines[i+1+ia].split()
                    self._atoms_positions[ia] = [info_atoms_pos[0],info_atoms_pos[1:4]]

                self._convert_to_bohr(pos_units)

            if "K_POINTS" in line:
                self._kpoints = lines[i+1]
            if "="  not in line:
                continue
            flag, value=line.split("=")
            for word in self._control:
                if flag.strip().lower() == word:
                    self._control[word] = self._read_word(value)
            for word in self._electrons:
                if flag.strip().lower() == word:
                    self._electrons[word] = self._read_word(value)
            for word in self._system:
                if flag.strip().lower() == word:
                    self._system[word] = self._read_word(value)
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
        #        for line in lines[i+2:i+2+int(self._system['ntyp'])]:
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
            self._rescale_pos(float(self._system['celldm(1)']))
        elif units == "crystal":
            vec = np.zeros((3,3),dtype=float);    i=0
            atom_pos=         np.zeros((int(self._system['nat']),3), dtype=float)
            atom_pos_bohr=[[],np.zeros((int(self._system['nat']),3), dtype=float)]
            for cp in self._cell_parameters:
                vec[i,:]= np.asarray(cp[:], dtype=np.float) #in **Bohr**
                i=i+1
            for ia in range(int(self._system['nat'])):
                atom_pos[ia,:] = np.asarray(self._atoms_positions[ia][1:4], dtype=np.float)
                atom_pos_bohr[ia] = [self._atoms_positions[ia][0],np.dot(vec,atom_pos[ia]).tolist()]
                self._atoms_positions[ia]=atom_pos_bohr[ia]

    def get_cell_parameters(self):
        cell = np.zeros((3,3),dtype=float)
        for cp in self._cell_parameters:
                cell[i,:]= np.asarray(cp[:], dtype=np.float)

    def _rescale_pos(self, scale):
        for ia in range(int(self._system['nat'])):
            self._atoms_positions[ia][1:4]=scale*np.asarray(self._atoms_positions[ia][1:4], dtype=np.float)

    def get_starting_positions(self):
        start_atom_pos= [[],np.zeros((int(self._system['nat']),3), dtype=float)]
        for ia in range(int(self._system['nat'])):
            start_atom_pos[ia] = self._atoms_positions[ia]
        return start_atom_pos

    def set_positions(self,positions):
        if isinstance(positions, np.ndarray):
            for ia in range(int(self._system['nat'])):
                self._atoms_positions[ia][1:4]=positions[ia]
        else:
            print('passo qui')
            self._atoms_positions = positions


# Recently assigned and default values can be shown thanks to the 
# following function
#
    def show_values(self):
        for word in self._control:
            print (word+"="+str(self._control[word]))
        for word in self._system:
            print (word+"="+str(self._system[word]))
        for word in self._electrons:
            print (word+"="+str(self._electrons[word]))


