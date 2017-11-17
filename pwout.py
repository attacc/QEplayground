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

        nk_pattern=r'number of k points=\s*([0-9]+)'
        nk_regexp =re.compile(nk_pattern)

        for line in lines:
            match = en_regexp.match(line)
            if match is not None: self.tot_energy=float(match.group(1))
            match = nk_regexp.match(line)
            if match is not None:
                print("PASSO QUI ")
                self.nkpoints  =int(match.group(1))

        print(self)

        self.read_bands(pw_output)

    def read_bands(self, pw_output):
        lines=iter(pw_output)

        for line in lines:
            if "End of self-consistent calculation" in line: break

        for ik in range(0,self.nkpoints):
            print("Ciao ")
        exit(0)

    def __str__(self):
        string=''
        string+='Total energy       :  %12.8f Ry\n' % self.tot_energy
        string+='Number of k-points :  %d \n' % self.nkpoints

        return string

