#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:01:00 2020

@author: moyuan
"""

import gsd
import gsd.hoomd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def repo20(pos, Ly):
    zero_loc= -Ly/2
    for i in pos:
        i[1]=i[1]-zero_loc
        
        
#%%avg sfvz
names=['th056a','th1a','th25a','th5a']
for name in names:
    data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
    frame=data[-1]
    N=frame.particles.N
    Ly = frame.configuration.box[1]
    da=0.8
    
    spacing=2
    z = np.arange(1,frame.configuration.box[1],spacing)
    avg_n=10000
    count = np.zeros(len(z)-1)
    a_count = np.zeros(len(z)-1)
    a_ind=[]
    p_ind=[]
    
    for p in range(0,N):
        if frame.particles.typeid[p]==1:
            a_ind.append(p)
        else:
            p_ind.append(p)
    for frame in data[-avg_n:]:
    
    
        pos = []
        a_pos=[]
        for p in p_ind:
            pos.append(frame.particles.position[p])
        for p in a_ind:
            a_pos.append(frame.particles.position[p])
        repo20(pos,Ly)
        repo20(a_pos,Ly)
        for p in range(0,len(pos)):
            for i in range(0,len(z)-1):
                if pos[p][1]>z[i] and pos[p][1]<z[i+1]:
                    count[i] += 1
                    break;
        for p in range(0,len(a_pos)):
            for i in range(0,len(z)-1):
                if a_pos[p][1]>z[i] and a_pos[p][1]<z[i+1]:
                    a_count[i] += 1
                    break;
    
        
    A_rec = spacing*frame.configuration.box[0]
    phi = []
    for i in count:
        phi.append(i*np.pi*0.25/(A_rec*avg_n))
    for i in range(0, len(a_count)):
        phi[i]+=a_count[i]*np.pi*0.25*da**2/(A_rec*avg_n)
            
    z_avg=np.arange(spacing/2,frame.configuration.box[1]-spacing/2,spacing)
    
    df = pd.DataFrame({'z':z_avg, 'phi':phi})
    df.to_csv('phivz_active{0}.csv'.format(name))
os.system('spd-say "your program has finished"')

#%%
name = 'th056a'
df=pd.read_csv('phivz_active{0}.csv'.format(name))
plt.plot(df['z'],df['phi'])

#%%
# #%%
# f=plt.plot(title='phi vs z averaged over last {0} seconds {1}'.format(avg_n, name))
# f.circle(z_avg, phi, fill_color='tomato',size=8, legend_label='l_c = 2 sigma')
# f.xaxis.axis_label='z (sigma)'
# f.yaxis.axis_label='phi (surface fraction)'
# # f.circle(phi, z_avg, fill_color='tomato',size=8,legend='g=0.05')
# # f.xaxis.axis_label='phi (surface fraction)'
# # f.yaxis.axis_label='z (sigma)'



# # show(f)
# export_png(f, filename='figure/sed_active{0}.png'.format(name))

# df = pd.DataFrame({'z':z_avg, 'phi':phi})
# df.to_csv('phivz_active{0}.csv'.format(name))