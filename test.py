#!/usr/bin/python3

from pwscf import *
from matdyn import *


qe_input=Pwscf("C.scf")
qe_input.write("test.in")
#qe_dyn=Matdyn(qe_input,"dynmat.out")

#qe_dyn.write_modes("test.modes")
#qe_new_input=qe_dyn.generate_displacement(0,3,0.01)

#qe_new_input.run("TEST1")
