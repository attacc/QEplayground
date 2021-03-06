from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.map_phononsPierre import *

scf_filename    ="hBN.scf.in"
dynmat_filename ="dynmat_all.eig"

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"
#pw="/home/elena/sources/qe-6.2/bin/pw.x"

#Serial job
qe_input.set_run_options(pw=pw)

#Parallel job
qe_input.set_run_options(pw=pw, nprocs=2, npool=2)


# Pseudo-potential directory
qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"
#qe_input.control['pseudo_dir']="'/home/elena/Research/pseudo'"

ff=[2,2,1]
new_supercell_name='hBN.supercell2x2.scf.in'
new_dynmat_name   ='dynmat_2x2.eig'

my_map = map_phonons(qe_input, qe_dyn, ff, new_supercell_name, new_dynmat_name)

translation_vectors = my_map.get_translation_vectors()

my_map.build_mapping(sort_ph=False,print_eig=False,norm_eig=True)
