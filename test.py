#!/usr/bin/python3

from pwscf import *
from matdyn import *


qe_input=Pwscf("diamond.scf.in")
qe_input.write("test.in")

print(qe_input)
