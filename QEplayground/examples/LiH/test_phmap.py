from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.ph2supercell import *


scf_filename    ="lih.scf.in"
dynmat_filename ="joint_dynmat.out"

q_grid=[2,1,1]

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

ph_map=map_phonons(qe_input, qe_dyn, q_grid)

ph_map.build_mapping()


