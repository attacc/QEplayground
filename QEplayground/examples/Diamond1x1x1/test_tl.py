from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.thermal_lines import *

scf_filename  ="diamond.scf.in"
eigv_filename ="diamond.eigv"

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,eigv_filename)

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"

#Serial job
qe_input.set_run_options(pw=pw)

#Parallel job
qe_input.set_run_options(pw=pw, nprocs=2, npool=2)


# Pseudo-potential directory
qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"

single_mode_thermal_line(qe_input, qe_dyn, modes=[3,4,5])

