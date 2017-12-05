#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#

import sys
import numpy as np
import re
from QEplayground.pwxml import *
from QEplayground.utilities import r_or_d

class Pwout:

    def __init__(self, qe_input):
        self.tot_energy   = 0.0
        self.nkpoints     = 0
        self.nel          = 0
        self.nbnd         = 0
        self.qe_input     =qe_input

    def read_output(self, filename, path="./"):
        """
        read QE output
        """
        f = open(path+"/"+filename,'r')
        pw_output = f.readlines()
        f.close()

        lines=iter(pw_output)

        en_pattern=r'!\s*total energy\s*=\s*('+r_or_d+') Ry'
        en_regexp =re.compile(en_pattern)

        nk_pattern=r'\s*number of k points=\s*([0-9]+)'
        nk_regexp =re.compile(nk_pattern)

        nel_pattern=r'\s*number of electrons\s*=\s*([0-9]+)'
        nel_regexp =re.compile(nel_pattern)

        for line in lines:
            match = en_regexp.match(line)
            if match is not None:
                self.tot_energy=float(match.group(1))/2.0 # from Rydberg to Hartree
            match = nk_regexp.match(line)
            if match is not None: self.nkpoints  =int(match.group(1))
            match = nel_regexp.match(line)
            if match is not None: self.nel       =int(match.group(1))
        #
        # Now reads the XML
        #
        prefix=self.qe_input.control['prefix'].strip('\'')
        self.qe_xml=PwXML(prefix=prefix,path=path)
  
    def get_dft_energy(self, kp, band):
        eigen=np.array(self.qe_xml.eigen)
        dft_energy = eigen[kp, band]
        return dft_energy

    def find_the_gap(self):
        ival =int(self.nel/2)-1
        icond=int(self.nel/2)
        if self.qe_xml.nbands < icond:
            print("Too few bands for the gap!! ")
            return
        eigen=np.array(self.qe_xml.eigen)
        top_val    =np.amax(eigen[:,ival])
        bottom_cond=np.amin(eigen[:,icond])

        ind_gap=bottom_cond-top_val
        dir_gap=np.amin(eigen[:,icond]-eigen[:,ival])

        return ival, icond, dir_gap,ind_gap



    def __str__(self):
        string=''
        string+='Total energy       :  %12.8f Ry\n' % self.tot_energy
        string+='Number of k-points :  %d \n' % self.nkpoints
        string+='Number of electrons:  %d \n' % self.nel

        return string

