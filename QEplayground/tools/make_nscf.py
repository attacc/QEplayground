#!/usr/bin/python3
#
# Copyright (c) 2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
import glob
from QEplayground.pwscf  import *

#path='hBN.supercellx2.scf.in_TL*'
path='hBN.scf.in_M*'

for QEfile in glob.glob(path):
    qe_input =Pwscf(QEfile)
    #
    # Setupt the NSCF
    #
    qe_input.control['restart_mode']=None
    qe_input.control['calculation']="'nscf'"

    qe_input.system['nbnd']=40
    qe_input.system['force_symmorphic']='.true.'

    qe_input.electrons['mixing_beta']=None
    qe_input.electrons['mixing_mode']=None
    qe_input.electrons['conv_thr']   =None
    qe_input.electrons['diago_thr_init']='1e-6'

    qe_input.kpoints[2]=2

    ofile=QEfile.replace("scf","nscf")

    qe_input.write(filename=ofile)
