#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 10:46:19 2020

@author: moyuan
"""

import gsd
import gsd.hoomd
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir('/home/moyuan/simulation/hoomd/units')

with open('data/langevin_msd.txt',mode='r+b') as handle:
    l_msds = pickle.load(handle)
    
with open('data/brownian_msd.txt',mode='r+b') as handle:
    b_msds = pickle.load(handle)
    
# lmsd = np.zeros(len(l_msds[0]))
# for msd in l_msds:
#     lmsd+=msd
# lmsd = lmsd/len(l_msds)

# bmsd = np.zeros(len(b_msds[0]))
# for msd in b_msds:
#     bmsd+=msd
# bmsd = bmsd/len(b_msds)

t = np.arange(0,200)
plt.figure()
plt.title('Surface Fraction 0.67 Colloidal Simulation Langevin and Brownian Dynamics MSD')
plt.plot(t, l_msds[0],'x',label='langevin dt=1e-4',c='b')
plt.plot(t, l_msds[1],'x',c='b')
plt.plot(t, l_msds[2],'x',c='b')
plt.plot(t,b_msds[0][:len(t)],'x',label='brownian dt=1e-3',c='r')
plt.plot(t,b_msds[1][:len(t)],'x',c='r')
plt.plot(t,b_msds[2][:len(t)],'x',c='r')
plt.legend()
plt.ylabel('MSD (sigma^2)')
plt.xlabel('t (s)')
plt.show()
