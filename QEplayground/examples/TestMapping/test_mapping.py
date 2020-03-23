from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.thermal_lines import *

scf_filename    ="lih.scf.in"
dynmat_filename ="lih.eigv"

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

single_mode_thermal_line(qe_input, qe_dyn, modes=[3,4,5])


