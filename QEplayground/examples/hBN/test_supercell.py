from QEplayground.pwscf  import *
from QEplayground.supercell  import *

scf_filename    ="hBN.scf.in"
qe_input =Pwscf(scf_filename)

super_c = supercell(qe_input,[6,6,2])

qe_s = super_c.write()

qe_s.write("hBN.supercell.scf.in")
