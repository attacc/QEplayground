#!/usr/bin/python3
# Copyright (C) 2018 C. Attaccalite and E. Cannuccia
#


from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.ph2supercell import *

scf_filename  ="diamond.scf.in"
dynmat_filename ="matdyn.modes"

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

ph_map=map_phonons(qe_input, qe_dyn,q_grid=[2,2,2])
ph_map.build_mapping()

