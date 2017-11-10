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

    def read_output(self, filename):
        """
        read QE output
        """
        f = open(self.filename,'r')
        pw_output = f.readlines()
        f.close()

        energy_pattern=r'!\s*total energy\s*=\s*('+r_or_d+')\s*Ry'
        match=re.search(pw_output, lines, re.MULTILINE)
        self.tot_energy=float(match.group(1))

    def __str__(self):
        string=''
        string+='Total energy :  %12.8f \n' % self.tot_energy

