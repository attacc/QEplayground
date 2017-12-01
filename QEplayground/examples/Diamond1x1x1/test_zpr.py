from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.pwout  import *
from QEplayground.units  import ha2ev
from QEplayground.zpr    import *
import math

r_order=3    # Richardson extrapolation order 
delta=0.0172 # Displacement in a.u.

scf_filename    ="diamond.scf.in"
dynmat_filename ="dynmat.out"

qe_input  = Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

#pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"
pw="/home/elena/sources/qe-6.2/bin/pw.x"

#Serial job
#qe_input.set_run_options(pw=pw)

#Parallel job
qe_input.set_run_options(pw=pw, nprocs=4, npool=2)


# Pseudo-potential directory
#qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"
qe_input.control['pseudo_dir']="'/home/elena/Research/pseudo'"


# Setup the band gap calculation
qe_input.electrons['diago_full_acc']='.true.'
qe_input.system['nbnd']             =7

zpr(qe_input, qe_dyn, delta, r_order=r_order)
