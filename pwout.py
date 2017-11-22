#!/usr/bin/python3
#
# Copyright (c) 2017, E. Cannuccia and C. Attaccalite
# All rights reserved.
#

import sys
import numpy as np
import re
from utilities import r_or_d

class Pwout:

    def __init__(self, qe_input):
        self.tot_energy   = 0.0
        self.nkpoints     = 0
        self.nel          = 0
        self.nbnd         = 0
        self.qe_input     =qe_input

    def read_output(self, filename):
        """
        read QE output
        """
        f = open(filename,'r')
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
            if match is not None: self.tot_energy=float(match.group(1))
            match = nk_regexp.match(line)
            if match is not None: self.nkpoints  =int(match.group(1))
            match = nel_regexp.match(line)
            if match is not None: self.nel       =int(match.group(1))

        if self.qe_input.system['nbnd'] != None:
            self.nbnd=int(self.qe_input.system['nbnd'])
        else:
        # Only valid for insulators without spin-orbit coupling
            self.nbnd=self.nel/2

    def __str__(self):
        string=''
        string+='Total energy       :  %12.8f Ry\n' % self.tot_energy
        string+='Number of k-points :  %d \n' % self.nkpoints
        string+='Number of electrons:  %d \n' % self.nel

        return string

