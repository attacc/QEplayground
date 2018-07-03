#!/usr/bin/python3
#
# Copyright (c) 2018, E. Cannuccia and C. Attaccalite
# All rights reserved.
#
import numpy as np
import os
from math import sqrt

amu2au  = 1822.888486192

eps_name="o.eps_q1_ip"
noph_folder="NOPH"
prefix="bn"
sigma_fname="sigma.dat"
delta=1.0
#delta=41.68640593/sqrt(amu2au)

sigma_file=open(sigma_fname,'r')

modes=np.genfromtxt(sigma_fname,delimiter=",")
print("\n Analize second derivatives\n")
print("Phonon modes : "+str(modes[:,0]))
print("Phonon sigma^2 : "+str(modes[:,1]))

eps_noph=np.genfromtxt(os.path.join(noph_folder+"/"+prefix+".save",eps_name))

for mode,sigma2 in zip(modes[:,0],modes[:,1]):

    folder_l="M"+str(int(mode)).zfill(3)+"_L1/"+prefix+".save"
    folder_r="M"+str(int(mode)).zfill(3)+"_R1/"+prefix+".save"

#    folder_l=noph_folder+"/"+prefix+".save"
#    folder_r=noph_folder+"/"+prefix+".save"

    file_r=os.path.join(folder_r,eps_name)
    file_l=os.path.join(folder_l,eps_name)

    if not os.path.exists(file_l) or not os.path.exists(file_r):
        continue

    eps_l=np.genfromtxt(os.path.join(folder_l,eps_name))
    eps_r=np.genfromtxt(os.path.join(folder_r,eps_name))

    eps_snd_der=(eps_r+eps_l-2.0*eps_noph)/(delta**2)/amu2au
    eps_snd_der=eps_snd_der*sigma2/2.0  # 2 comes from Taylor expantion f''/2 * q^2
    eps_snd_der[:,0]=eps_noph[:,0] # restore energy line

    np.savetxt("eps_der_"+str(mode)+".dat",eps_snd_der)
