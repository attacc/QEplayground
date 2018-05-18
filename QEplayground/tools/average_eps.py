#!/usr/bin/python
import argparse
import numpy as np
import sys

"""
Average dielectric constant on different configurations
Author:  Claudio Attaccalite
"""

def exit_error(stringa):
    print(stringa)
    sys.exit(1)
#
# parse command line
#
parser = argparse.ArgumentParser(prog='average_eps',description='Average dielectric constant on different configurations',\
         epilog="Copyright Claudio Attaccalite 2018")
parser.add_argument('--files' , nargs='*', default=None)
args = parser.parse_args()

print("\n Average EPS files \n\n")


if (args.files == None): exit_error('type "average_eps.py --help" for help ',)

print(" Files to join: %s  \n" % str(args.files)) 

# get the number of column from the first file

data=np.genfromtxt(args.files[0],comments="#")

result=np.zeros([data.shape[0],5],float)

new_shape=[len(args.files),data.shape[0],data.shape[1]]

alldata=np.zeros(new_shape,float)

for ifile,filename in enumerate(args.files):
    try:
        data=np.genfromtxt(filename,comments="#")
    except:
        exit_error(" Error opening file : " + filename)
    alldata[ifile,:,:]=data

for idata in range(2):
    result[:,idata+1] =np.mean(alldata[:,:,idata+1],axis=0)
    result[:,idata+3] =np.std( alldata[:,:,idata+1],axis=0)

result[:,0]=alldata[0,:,0]

np.savetxt('average_xhi.out',result,fmt='%2.15e')


