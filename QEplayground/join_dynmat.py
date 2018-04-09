#!/usr/bin/python3
#
# Copyright (c) 2018, Claudio Attaccalite and Elena Cannuccia
# All rights reserved.
#
#
#
import argparse

"""
Joint dynmat output files for different q-vectors in a single one
Author:  Claudio Attaccalite
"""
#
# parse command line
#
parser = argparse.ArgumentParser(prog='join_dynmat',description='Joint dynmat output files for different q-vectors in a single one',epilog="Copyright Claudio Attaccalite 2018")

parser.add_argument('--files' , nargs='*', default=None)
args = parser.parse_args()

print("\n Joint dynmat.x output files \n\n")

if (args.files == None): exit_error('type "join_dynmat.py --help" for help ',)

print(" Files to join: %s  \n" % str(args.files))






