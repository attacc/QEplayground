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
        self.nat        = int(self.qe_input.system['nat'])
        self.nmodes     = 3*int(self.nat)
        self.qvec       = []
        #
        if(filename):
            self.read_dyn_output(filename)

    def __del__(self):
        print("Destroy class Dynmat")

    def __str__(self):
        string=""
        for iq in range(0,self.nqpoints):
            string+="q = \n"
            for im in range(0,self.nmodes):
                string+="\tmode (   %d) \n" % (im+1)
                for ia in range(0,self.nat):
                    string+="\t\t "+str(self.pol_vec[iq][im][ia])+"\n"
        return string

    def read_dyn_output(self, filename):
        ifile = open(filename,'r')
        self.file_lines=ifile.readlines()
        lines = iter(self.file_lines)
        ifile.close()
        #
        # Find the number of q-points
        #
        qpattern="q ="
        self.nqpoints=sum(qpattern in line for line in lines)
        #
        g
#        for line in lines:
#            if qpattern in line:
#                qline, vec = line.split()
#                qvec.append(vec)
#                next(lines) # the line with asterix
#                for im in range(0,nmodes):
#                    for ia in range(self.nat):
#                        tmp_vec     =re.findall(r_or_d, next(lines))
#                        vec=tmp_vec[0::2]+1j*tmp_vec[1::2]


