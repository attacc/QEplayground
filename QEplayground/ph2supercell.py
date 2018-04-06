#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
#

class map_phonons():
    """
    Map phonon in a supercell
    """
    def __init__(self,qe_input, qe_dyn, Q_grid):
        self.qe_input=qe_input
        self.qe_dyn  =qe_dyn
        print("Warning this code works only without symmetries!!! nosym=.true. and noinv=.true. ")
    
    def build_mapping(self):
        n_qpoints=product(Q_
