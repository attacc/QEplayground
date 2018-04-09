#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
#
from itertools import product


class map_phonons():
    """
    Map phonon in a supercell
    """
    def __init__(self,qe_input, qe_dyn, q_grid):
        self.qe_input=qe_input
        self.qe_dyn  =qe_dyn
        self.q_grid  =q_grid
        print(" \n\n\n * * * Map phonons in a supercell * * *\n")
        print(" WARNING! This code works only without symmetries!!! nosym=.true. and noinv=.true. \n")

    def build_mapping(self):
        self.n_qpoints=self.q_grid[0]*self.q_grid[1]*self.q_grid[2]
        print(" Q grid "+str(self.q_grid))
        print(" Number of q-points "+str(self.n_qpoints))

        if(self.qe_dyn.nqpoints != self.n_qpoints):
            print("\n ERROR! number of q-points in matdyn.modes different from the input q-grid!! \n")
            exit(0)

