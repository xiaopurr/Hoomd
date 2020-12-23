#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 10:02:28 2020

@author: moyuan
"""


import hoomd
from hoomd import md
import numpy as np

import os


def simulate(path,name, time,alpha,theta):
    hoomd.context.initialize()
    seed=np.random.randint(0,9999)
    try:
        os.chdir(str(path))
    except:
        print('path error: check if directory exist')  
    
    name=name
    #unit parameters
    m_u=17.05

    sigma_u=2.79
    e_u = 4.11
    tau_u = np.sqrt(m_u*sigma_u**2/e_u)
    g = 9.8#m/s
    f_u=m_u*sigma_u/(tau_u**2)
    #basic parameters
    kT=1
    sigma=1
    gamma=23400*tau_u/m_u
    time_real = 100e3#480e3#miliseconds
    time_relax = time*1e3#miliseconds
    dt = 0.05
    tsteps = time_real/(tau_u*dt)
    tsteps_relax = time_relax/(tau_u*dt)
    #activity parameters
    # alpha=2/100
    # #force field constant
    # alpha = 0.56*360/(2*np.pi)
    # fg=-0.116
    fg=-20
    epsilon=10
    theta = theta
    fg_real = m_u*(g)*np.sin(theta*np.pi/180)
    fg_relax=-fg_real/f_u
    #make the unit box
    L=1.3*sigma
    box = hoomd.data.boxdim(Lx = L, Ly= L, Lz=1, xy = 0, xz = 0, yz = 0,
                            dimensions=2)
    
    snapshot_init = hoomd.data.make_snapshot(N = 1, box = box, 
                                             particle_types=['A','active'])
    
    snapshot_init.particles.position[:] = [[0,0,0]]
    #replicate it 100 times
    box_dim = [20,40]
    snapshot_init.replicate(box_dim[0],box_dim[1],1)
    N=int(snapshot_init.particles.N)
    N_a = int(N*alpha)
    # snapshot_init.particles.typeid[0:N_a] = 1;
    # np.random.shuffle(snapshot_init.particles.typeid)
    sys = hoomd.init.read_snapshot(snapshot_init)
    
    groupA = hoomd.group.type('A', update=True)
    groupactive = hoomd.group.type('active',update=True)
    nl = md.nlist.cell()
    wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
    wca.pair_coeff.set('A', 'A', epsilon=epsilon, sigma=sigma)
    wca.pair_coeff.set('A', 'active', epsilon=epsilon, sigma=sigma)
    wca.pair_coeff.set('active', 'active', epsilon=epsilon, sigma=sigma)
    # wca.pair_coeff.set('A', 'A', epsilon=10.0, sigma=1.0)
    wca.set_params(mode= 'shift')
    
    walls=md.wall.group()
    
    walls.add_plane((0,-L*box_dim[1]/2,0), (0,1,0))
    walls.add_plane((0,L*box_dim[1]/2,0),(0,-1,0))
    
    lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
    lj_wall.force_coeff.set('A', epsilon=10, sigma=sigma/2)
    lj_wall.force_coeff.set('active', epsilon=10, sigma=sigma/2)
    # lj_wall.force_coeff.set('i', epsilon=0, sigma=0)
    
    const_f = md.force.constant(fvec=(0,fg,0),group=groupA)
    const_a = md.force.constant(fvec=(0,fg,0),group=groupactive)
    integrator = md.integrate.brownian(group=hoomd.group.all(), kT=kT, seed=seed)
    integrator.set_gamma('A', gamma=gamma)
    integrator.set_gamma('active', gamma=gamma)
    md.integrate.mode_standard(dt)
    
    hoomd.dump.gsd('data/'+name+'.gsd', period=int(1/(tau_u*1e-3*dt)), 
                    group=hoomd.group.all(),overwrite = True,dynamic=['momentum'])
    
    hoomd.run(tsteps=tsteps)
    const_f.set_force(fvec=(0,fg_relax,0),group=groupA)
    const_a.set_force(fvec=(0,fg_relax,0),group=groupactive)
    hoomd.run(tsteps=tsteps_relax)

#%%
# import sys
# path = input('path:')
path = '/home/mike_chen/data/sed'
name = input('name: ')
time = float(input('time(seconds): '))
# alpha = float(input('alpha: '))
alpha=0.05
theta = float(input('sed angle (degree): '))
simulate(path, name, time, alpha,theta)

