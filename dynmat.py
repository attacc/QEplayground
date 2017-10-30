# Copyright (C) 2017 Elena Cannuccia, Claudio Attaccalite
# All rights reserved.
#
import numpy as np
import re, sys
import math
import pwscf
from utilities import r_or_d


class dynmat:
    """
    A class to read dynmat.x output (only at Gamma point for the moment)
    """
    def __init__(self, qe_input, filename=None):
        self.qe_input=qe_input
        self.freq_pattern=r'\s*freq \(\s*(\d*)\) =\s*('+r_or_d+')\s\[THz\]'
        self.phonons    = []
        self.nat        = int(self.qe_input.system['nat'])
        self.n_modes    = 3*int(self.nat)
        self.nqpoints   = 1  #Only at gamma for the momemt
        if(filename):
            self.read_dyn_output(filename)
        

    def __del__(self):
        print("Destroy class Dynmat")

    def read_dyn_output(self, filename):
        ifile = open(filename,'r')
        self.file_lines=ifile.readlines()
        lines = iter(self.file_lines)
        ifile.close()
        regex=re.compile(self.freq_pattern)
        self.pol_vec=np.zeros((self.nqpoints,self.n_modes,self.nat,3),dtype=complex)

        for line in lines:
            match=regex.match(line)
            if match is not None:
                imode=int(match.group(1))-1  # Mode number
                for ia in range(self.nat):
                    tmp_vec     =re.findall(r_or_d, next(lines))
                    self.pol_vec[0][imode][ia].real=tmp_vec[0::2]
                    self.pol_vec[0][imode][ia].imag=tmp_vec[1::2]

