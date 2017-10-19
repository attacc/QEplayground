#
# Copyright (c) 2015, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy projec
#
# Modified by C. Attaccalite and E. Cannuccia
#
#
import os
import re
import math
import numpy as np
from pwscf import Pwscf

class supercell():
    """A class to generate custom supercells from a quantum espresso input file
    """

