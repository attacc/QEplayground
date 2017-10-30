# Copyright (C) 2017 Elena Cannuccia, Claudio Attaccalite
# All rights reserved.
#
import numpy as np
import re, sys
import math
import pwscf


class dynmat:
    """
    A class to read dynmat.x output 
    """
    def __init__(self, qe_input, filename=None):
        self.qe_input=qe_input
        if(filename):
            self.read_dyn_output(filename)

    def __del__(self):
        print("Destroy class Dynmat")

    def read_dyn_output(self, filename):
        ifile = open(filename,'r')
        self.file_lines=ifile.readlines()
        ifile.close()

