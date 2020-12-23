#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:15:52 2020

@author: moyuan
"""

import hoomd
from hoomd import md
import numpy as np
import gsd
import gsd.hoomd
from bokeh.plotting import figure, show
from bokeh.io import export_png
import matplotlib.pyplot as plt

m_u=17.05
sigma_u=2.79
e_u = 4.11
tau_u = np.sqrt(m_u*sigma_u**2/e_u)
#basic parameters
kT=1
sigma=1
# gamma=23400*tau_u/m_u
gamma=100
time_real = 200e3#miliseconds
dt = 1e-3
tsteps = time_real/(tau_u*dt)
epsilon=10
seed=np.random.randint(1000,9999)
a=2
box_dim=[3,6]
fg=-1100


hoomd.context.initialize()
#%%
uc = hoomd.lattice.sq(a = a, type_name='A')
sys = hoomd.init.create_lattice(uc, n=box_dim)

nl=md.nlist.cell()
wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
wca.pair_coeff.set('A', 'A', epsilon=0, sigma=0)


integrator = md.integrate.brownian(group=hoomd.group.all(), kT=kT, seed=seed)
integrator.set_gamma('A', gamma=gamma)
md.integrate.mode_standard(dt)

#active rotation
activity = []
A=100
theta = np.random.rand(box_dim[0]*box_dim[1])*2*np.pi
activity=[(np.cos(i)*A, np.sin(i)*A,0) for i in theta]
force = md.force.active(seed= seed, group= hoomd.group.all(),f_lst = activity,rotation_diff=0.01,
                        orientation_reverse_link=True, orientation_link=False)


#Walls at bottom and top
walls=md.wall.group()

walls.add_plane((0,-a*box_dim[1]/2,0), (0,1,0))
walls.add_plane((0,a*box_dim[1]/2,0),(0,-1,0))

lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
lj_wall.force_coeff.set('A', epsilon=epsilon, sigma=sigma)

#gravity
const_f = md.force.constant(fvec=(0,fg,0),group=hoomd.group.type('A'))
#%%
hoomd.dump.gsd(filename='/home/moyuan/simulation/hoomd/sedimentation/tests/ori_dist_{0}.gsd'.format(seed),
               period=int(1/(tau_u*1e-3*dt)), group=hoomd.group.all(), overwrite=True, dynamic=['momentum'])

hoomd.run(tsteps = tsteps)
#%%
import os
def quat_ro(q, vec):
    qr,qi,qj,qk = q
    mat = [[1-2*(qj**2+qk**2),2*(qi*qj-qk*qr),2*(qi*qk+qj*qr)],
            [2*(qi*qj+qk*qr),1-2*(qi**2+qk**2),2*(qj*qk-qi*qr)],
            [2*(qi*qk-qj*qr),2*(qj*qk+qi*qr),1-2*(qi**2+qj**2)]
            ]
    return np.matmul(mat, vec)

data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/ori_dist_{0}.gsd'.format(seed),mode='rb')
#%%
v=[[],[]]
for frame in data:
    for i in frame.particles.velocity:
        v[0].append(i[0])
        v[1].append(i[1])
plt.hist2d(v[0],v[1],bins=10)
plt.axes().set_aspect('equal')
#%%
plt.scatter([i**2 for i in v[0]],[i**2 for i in v[1]])
#%%
plt.figure()
for i in range(0,len(v[0])):
    th = np.arctan2(v[1][i],v[0][i])
    plt.scatter(th,np.sqrt(v[0][i]**2+v[1][i]**2))
#%%

# import pandas as pd
# a=0
ol = []
for frame in data:
    ori_q=frame.particles.orientation
    ori = np.zeros(len(ori_q)).tolist()
    v=frame.particles.velocity
    i=0
    for q in ori_q:
        ori[i]=(quat_ro(q,[0,0,1])[:2])
        i+=1
    i=0
    for o in ori:
        ol.append([np.arctan2(o[1],o[0]),v[i]])
        i+=1
# #%%
# theta = [i[0] for i in ol]
# v = [np.sqrt(i[1][0]**2+i[1][1]**2) for i in ol]

# dth = np.arange(-np.pi,np.pi,0.3)
# count = np.zeros(int(len(dth)-1))
# vavg = np.zeros(int(len(dth)-1))
# for i in range(0,len(dth)-1):
#     for j in range(0,len(theta)):
#         if theta[j]>dth[i] and theta[j]<dth[i+1]:
#             vavg[i]+=v[j]
#             count[i]+=1
#             # break;
# for i in range(0,len(vavg)):
#     vavg[i]=vavg[i]/count[i]
# import matplotlib.pyplot as plt
# plt.figure()
# plt.plot(dth[:len(dth)-1],vavg,'x')
#%%
################################################################
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/ori_dist_{0}.gsd'.format(seed),mode='rb')
try:
    os.mkdir('/home/moyuan/simulation/hoomd/sedimentation/tests/snap_{0}'.format(seed))
except:
    pass;
j=0
ind = np.arange(0,len(data))
for fr in range(0,len(ind)):
    frame = data[int(ind[fr])]
    Lx, Ly = frame.configuration.box[:2]
    r = sigma/2
    offset =r
    p = frame.particles.position
    
    f = figure(x_range=[-Lx/2-offset,Lx/2+offset],y_range=[-Ly/2-offset,Ly/2+offset])
    f.circle([i[0] for i in p], [i[1] for i in p],radius=r, fill_alpha=0.1)

    ori_q=frame.particles.orientation
    ori = np.zeros(len(ori_q)).tolist()
    i=0
    for q in ori_q:
        ori[i]=(quat_ro(q,[0,0,1])[:2])
        i+=1
    i=0
    for o in ori:
        loc = [p[i][0],p[i][1]]
        f.line([loc[0],loc[0]+o[0]],[loc[1],loc[1]+o[1]])
        i+=1
    for n in range(0,frame.particles.N):
        tra_p=[]
        for c in range(0,fr):
            tra_p.append(data[int(ind[c])].particles.position[n])
        f.line([t[0] for t in tra_p],[t[1] for t in tra_p],color='tomato')
    j+=1
    export_png(f, filename='/home/moyuan/simulation/hoomd/sedimentation/tests/snap_{0}/track_o_{1}.png'.format(seed,j))