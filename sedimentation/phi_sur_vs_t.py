#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:10:47 2020

@author: moyuan
"""

import gsd
import gsd.hoomd
import numpy as np
import pandas as pd

nl = [3,4,5,6,7,8,9]
for n in nl:
    name='a_Da0{0}_alpha_05'.format(n)
    
    data = gsd.hoomd.open(
        '/home/moyuan/simulation/hoomd/sedimentation/data/active_size/{0}.gsd'.format(name))
    
    Ly = data[0].configuration.box[1]
    Lx = data[0].configuration.box[0]
    a_ind = []
    p_ind=[]
    for p in range(0,data[0].particles.N):
        if data[0].particles.typeid[p]==1:
            a_ind.append(p)
        else:
            p_ind.append(p)
    
    z = np.zeros(int(len(data)))
    t = np.zeros(int(len(data)))
    i=0
    tt=0
    for frame in data:
        for ind in a_ind:
            z[i]+=frame.particles.position[ind][1]
        t[i]=tt
        z[i]=z[i]/len(a_ind)
        i+=1
        tt+=1
    
    z = z+Ly/2
    
    spacing = 2
    phi = np.zeros(len(z))
    for z_i in range(0, int(len(z))):
        frame = data[z_i]
        z_phi_low = z[z_i]-spacing/2
        z_phi_high = z[z_i] + spacing/2
        count=0
        for p in p_ind:
            pos = frame.particles.position[p][1]
            if pos>z_phi_low and pos < z_phi_high:
                count+=1
                
        phi_temp = count*(np.pi*0.25)/(spacing*Lx)
        phi[z_i] = phi_temp
        
    df = pd. DataFrame(data={'t':np.arange(0,len(phi)), 'phi':phi})
    df.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/phivt_sigma_a_0{0}.csv'.format(n))
#%%
import matplotlib.pyplot as plt
plt.figure()
t = np.arange(0,20002)
for n in nl:
    df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/phivt_sigma_a_0{0}.csv'.format(n))
    
    
    plt.plot(df['t'], df['phi'],'.',label='Da = 0.{0}'.format(n))
    

plt.title('sf v t')
plt.xlabel('time(s)')
plt.ylabel('sufrace fraction around active particles')
plt.legend(bbox_to_anchor=(1,1))
plt.show()
# plt.savefig('/home/moyuan/simulation/hoomd/sedimentation/figure/sf_arnd_v_t/da0{0}.png'.format(n))
#%%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(t, np.exp(phi),'.')

plt.xlabel('log time')
plt.ylabel('exp(sf)')
plt.title('exp(sf) v t sigma a = 0.{0} sigma'.format(n))

plt.savefig('/home/moyuan/simulation/hoomd/sedimentation/figure/sf_arnd_v_t/exp_da0{0}.png'.format(n))