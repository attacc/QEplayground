#!/usr/bin/python3
#
# Copyright (c) 2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from QEplayground.utilities import r_or_d
import pylab as pl


pattern=r'!\s*total energy\s*=\s*('+r_or_d+') Ry'

path = "output_*"

nfiles=len(glob.glob(path))
en=np.zeros(nfiles,dtype=float)


for efile in enumerate(glob.glob(path)):
    f = open(efile[1],'r')
    txt = f.read()
    match = re.search(pattern, txt, re.MULTILINE)
    if not match:
        print(" Error reading "+str(efile[1]))
        exit(0)
    en[efile[0]]=float(match.group(1))


print("Average Energy %lf " % np.average(en))
err=float(np.std(en))/float(nfiles)
print("Error mean std %lf " % err)

# Histogram
#heights,bins = np.histogram(en,bins=20)

# Normalize
pl.hist(en,normed=True)
pl.show()
