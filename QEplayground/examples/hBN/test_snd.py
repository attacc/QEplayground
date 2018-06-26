from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.second_deriv import *

r_order=2    # Richardson extrapolation order 
delta=0.02  # Displacement in a.u.

scf_filename    ="hBN.scf.in"
dynmat_filename ="dynmat.eig_1"

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"

qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"

second_deriv(qe_input, qe_dyn, delta, r_order=r_order)

