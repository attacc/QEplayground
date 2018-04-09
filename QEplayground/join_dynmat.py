#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
#
import argparse
from QEplayground.pwscf  import *
from QEplayground.matdyn import *

"""
Joint dynmat output files for different q-vectors in a single one
Author:  Claudio Attaccalite
"""
#
# parse command line
#
parser = argparse.ArgumentParser(prog='join_dynmat',description='Joint dynmat output files for different q-vectors in a single one',epilog="Copyright Claudio Attaccalite 2018")

parser.add_argument('--scf' ,   type=str,  default=None)
parser.add_argument('--out' ,   type=str,  default="joint_dynmat.out")
parser.add_argument('--files' , nargs='*', default=None)
args = parser.parse_args()

print("\n Joint dynmat.x output files \n\n")

if (args.files == None or args.scf == None):
    print('type "join_dynmat.py --help" for help')
    exit(1)

print(" Self-consistent input: %s  \n" % str(args.scf))
print(" Files to join: %s  \n" % str(args.files))

qe_input =Pwscf(args.scf)

create_file=True

for dynmat_file in args.files:
    try:
        qe_dyn=Matdyn(qe_input,dynmat_file)

        if create_file:
           wmode='w'
           create_file=False
        else:
            wmode='a'

        qe_dyn.write_modes(args.out, wmode=wmode)
        del qe_dyn
    except:
        print("Error opening file: "+dynmat_file)
        exit(0)







