#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 21:36:15 2020

@author: moyuan
"""

import numpy as np
import hoomd
from hoomd import md
import gsd.hoomd


v0=4*(1/(1.1/(0.669e-3))) #4 micron / second
vs = v0/2               #vs/v0=0.5

gamma=4130
fr=gamma*(4/3)*(1/8)
Dr=1/fr

sigma=1
kT=0
TAU=0.669
t_real = 3600e3 #millisecond
dt=0.1
tsteps=int(t_real/(TAU*dt))

epsilon=10

a=20
box_dim=[5,30]
fg = -vs*gamma
fa = v0*gamma

seed=np.random.randint(0,99999)
hoomd.context.initialize()
#%%
def c_lambda(vs,v0,Dr):
    return (v0**2/(2*Dr*vs))*(1-(7/4)*(vs/v0)**2)
#%%
lambda_theoretical=c_lambda(vs,v0,Dr)

uc = hoomd.lattice.sq(a = a, type_name='A')
sys = hoomd.init.create_lattice(uc, n=box_dim)

nl=md.nlist.cell()
wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
wca.pair_coeff.set('A', 'A', epsilon=0, sigma=0)
i_mode = md.integrate.mode_standard(dt)
i_mode.set_params(dt=dt, aniso=False) 

integrator = md.integrate.brownian(group=hoomd.group.all(), kT=kT, seed=seed)
integrator.set_gamma('A', gamma=gamma)

#active rotation
activity = []
A=fa
theta = np.random.rand(box_dim[0]*box_dim[1])*2*np.pi
activity=[(np.cos(i)*A, np.sin(i)*A,0) for i in theta]
force = md.force.active(seed= seed, group= hoomd.group.all(),f_lst = activity,rotation_diff=Dr,
                        orientation_reverse_link=True, orientation_link=False)

#Walls at bottom and top
walls=md.wall.group()

walls.add_plane((0,-a*box_dim[1]/2,0), (0,1,0))
walls.add_plane((0,a*box_dim[1]/2,0),(0,-1,0))

lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
lj_wall.force_coeff.set('A', epsilon=epsilon, sigma=sigma/2)

#gravity
const_f = md.force.constant(fvec=(0,fg,0),group=hoomd.group.type('A'))

hoomd.dump.gsd(filename='/home/moyuan/simulation/hoomd/sedimentation/data/lambda_check_vsdv01d{0}.gsd'.format(int(v0/vs)),
               period=int(1/(TAU*1e-3*dt)), group=hoomd.group.all(), overwrite=True, dynamic=['momentum'])

hoomd.run(tsteps = tsteps)

#%%
import matplotlib.pyplot as plt
import pandas as pd
def sedimentation_profile(frames, spacing):
    frame0=frames[0]
    box = frame0.configuration.box[:2]
    N=frame0.particles.N
    Lx = box[0]
    Ly = box[1]
    
    z = np.arange(0, Ly, spacing)
    count=np.zeros(len(z)-1)
    for frame in frames:
        y = [p[1]+Ly/2 for p in frame.particles.position]
        for p in y:
            for i in range(0,len(z)-1):
                if p>=z[i] and p<=z[i+1]:
                    count[i]+=1
                    break;
    de = count*(4*np.pi*(0.5)**3/3)/(spacing*Lx*len(frames))
    return de
eq_t=3000#seconds
max_length=100
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/data/lambda_check_vsdv01d{0}.gsd'.format(int(v0/vs)))
de = sedimentation_profile(data[-eq_t:],1)

z0=np.arange(0,max_length,1)
plt.plot(z0,de[0:max_length],'o',label='simulation')

zfit=np.linspace(0,max_length,1000)
A_fit = de[0]/np.exp(-0.5/lambda_theoretical)
rho_fit = [A_fit*j for j in [np.exp(-i/lambda_theoretical) for i in zfit]]


# from scipy.optimize import curve_fit
# def rho(z, A):
#     return A*np.exp(-z/lambda_theoretical)
# popt, pcov = curve_fit(rho,z0,de[:20])

plt.plot(zfit,rho_fit,label=r'${0}exp(-z/{1})$'.format(round(A_fit,2),round(lambda_theoretical,2)))
plt.xlabel('z (sigma)')
plt.ylabel('surface fraction')
plt.title(r'$v_s/v_0 = 0.25$, no translational noise') 
plt.legend()
#%%log lin
plt.figure()
log_rho = [np.log(i) for i in de[0:max_length]]
plt.plot(z0,log_rho,'o')

log_rho_fit=[np.log(i) for i in rho_fit]
plt.plot(zfit[:290],log_rho_fit[:290])

plt.xlabel('z (sigma)')
plt.ylabel(r'log($\rho$)')
plt.title(r'$\frac{v_s}{v_0}=0.33$, no translation noise, with interaction')
#%%
frame = data[-1]
pos=frame.particles.position
plt.scatter([i[0] for i in pos],[i[1] for i in pos])
#%%polarization
def quat_ro(q, vec):
    qr,qi,qj,qk = q
    mat = [[1-2*(qj**2+qk**2),2*(qi*qj-qk*qr),2*(qi*qk+qj*qr)],
            [2*(qi*qj+qk*qr),1-2*(qi**2+qk**2),2*(qj*qk-qi*qr)],
            [2*(qi*qk-qj*qr),2*(qj*qk+qi*qr),1-2*(qi**2+qj**2)]
            ]
    return np.matmul(mat, vec)
l=[]
Ly=data[0].configuration.box[1]
for frame_ind in range(-600,-1):
    # for o in frame.particles.orientation:
    #     tmp = quat_ro(o,[0,0,1])
    #     l.append(tmp)
    
    # for v in frame.particles.velocity:
    #     l.append(v)
    f0 = data[frame_ind]
    f1 = data[frame_ind+1]
    for pa in range(0,20):
        pos0=f0.particles.position[pa]
        if pos0[1]>=-Ly/2+1:
            pos1=f1.particles.position[pa]
            img0 = f0.particles.image[pa]
            img1=f1.particles.image[pa]
            
            pos0=pos0+img0*f0.configuration.box[:3]
            pos1=pos1+img1*f1.configuration.box[:3]
            
            dr = (pos1-pos0)
            
            l.append(dr)

plt.hist2d([i[0] for i in l],[i[1] for i in l],bins=20)
plt.axes().set_aspect('equal')
#%%
import bokeh
from bokeh.plotting import *
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/data/lambda_check_vs{0}_v0{1}_Dr{2}_kT{3}.gsd'.format(vs,v0,Dr,kT))
frame = data[-500]
# for frame in data[0]:
f = figure(match_aspect=True)
f.circle(frame.particles.position[:,0],frame.particles.position[:,1],radius=0.5)
show(f)

def quat_ro(q, vec):
    qr,qi,qj,qk = q
    mat = [[1-2*(qj**2+qk**2),2*(qi*qj-qk*qr),2*(qi*qk+qj*qr)],
            [2*(qi*qj+qk*qr),1-2*(qi**2+qk**2),2*(qj*qk-qi*qr)],
            [2*(qi*qk-qj*qr),2*(qj*qk+qi*qr),1-2*(qi**2+qj**2)]
            ]
    return np.matmul(mat, vec)
l=[]
for o in frame.particles.orientation:
    l.append(quat_ro(o,[0,0,1]))
# import matplotlib.pyplot as plt
for i in l:
    plt.plot(i[0],i[1],'o')

# show(f)
