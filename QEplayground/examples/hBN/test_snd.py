from QEplayground.pwscf  import *
from QEplayground.matdyn import *
from QEplayground.second_deriv import *
import math

r_order=2    # Richardson extrapolation order 
#delta=math.sqrt(1737.7564395947854)/math.sqrt(amu2au) #0.5  #sqrt(amu2au)=42.695298174 and so for delta=0.5
delta=1.0

scf_filename    ="hBN.scf.in"
dynmat_filename ="dynmat.eig_1"

qe_input =Pwscf(scf_filename)
qe_dyn=Matdyn(qe_input,dynmat_filename)

pw="/home/attacc/SOFTWARE/qe-6.1/bin/pw.x"

qe_input.control['pseudo_dir']="'/home/attacc/SOFTWARE/PSEUDO_PWSCF'"

second_deriv(qe_input, qe_dyn, delta, r_order=r_order,modes=[3,4])

