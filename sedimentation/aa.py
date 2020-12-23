#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 19:06:36 2020

@author: moyuan
"""

import hoomd
from hoomd import md
import numpy as np

import os

def simulate(path, restart_name, new_name, alpha, time,fga,da,theta):
    hoomd.context.initialize()
    seed=np.random.randint(0,9999)
    try:
        os.chdir(path+'/data/')
    except:
        os.mkdir(path+'/data/')
    os.chdir(path)
    name = restart_name
    active_alpha = alpha
    D_a = da
    name_a=new_name
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
    time_real = time*1e3#miliseconds
    dt = 2e-2
    tsteps = time_real/(tau_u*dt)
    #force field constant
    # alpha = 0.56*360/(2*np.pi)
    # fg=-0.116
    theta = theta
    fg_real = m_u*(g)*np.sin(theta*np.pi/180)
    fg=-fg_real/f_u
    fg_a = fg*fga
    epsilon=10
    L=1.3*sigma
    box_dim = [20,40]

    system = hoomd.init.read_gsd('data/'+name+'.gsd',frame=-1)
    for i in system.particles.types:
        if i == 'active':
            flag = 1
    if not flag:
        system.particles.types.add('active')
    
    snap = system.take_snapshot()
    Lx=snap.box.Lx
    Ly=snap.box.Ly
    N=snap.particles.N
    
    N_a = int(N*active_alpha)
    
    loc = snap.particles.position.tolist()
    loc.sort(key=lambda x: x[1])
    top_particles_ind=[]
    for i in range(0, len(loc)):
        if snap.particles.position[i][1] >= loc[-N_a][1]:
            top_particles_ind.append(i)
    
    for i in range(0,N_a):
        snap.particles.typeid[top_particles_ind[i]]=1
        
    system.restore_snapshot(snap)

    r_cut = np.power(2,1/6)
    nl = md.nlist.cell()
    wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
    wca.pair_coeff.set('A', 'A', epsilon=epsilon, sigma=sigma,r_cut=r_cut)
    wca.pair_coeff.set('A', 'active', epsilon=epsilon, sigma=(D_a+sigma)/2,
                       r_cut=r_cut*(D_a+sigma)/2)
    wca.pair_coeff.set('active', 'active', epsilon=epsilon, sigma=sigma*D_a,
                       r_cut=r_cut*D_a)
    # wca.pair_coeff.set('A', 'A', epsilon=10.0, sigma=1.0)
    wca.set_params(mode= 'shift')
    #Walls at bottom and top
    walls=md.wall.group()
    
    walls.add_plane((0,-L*box_dim[1]/2,0), (0,1,0))
    walls.add_plane((0,L*box_dim[1]/2,0),(0,-1,0))
    
    lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
    lj_wall.force_coeff.set('A', epsilon=10, sigma=sigma/2)
    lj_wall.force_coeff.set('active', epsilon=10, sigma=D_a/2,
			    r_cut=r_cut*D_a/2)
    # lj_wall.force_coeff.set('i', epsilon=0, sigma=0)
    #gravity
    const_f = md.force.constant(fvec=(0,fg,0),group=hoomd.group.type('A'))
    const_a = md.force.constant(fvec=(0,fg_a,0), group=hoomd.group.type('active'))
    
    #activity
    A = 1 #amplitude factor
    angles = np.random.rand(N_a)*2*np.pi
    activity = [(np.cos(i)*A, np.sin(i)*A,0) for i in angles]
    Dr = 0.1
    md.force.active(seed=seed, group = hoomd.group.type('active'), f_lst = activity, rotation_diff=Dr,
                    orientation_link=False)
    
    integrator = md.integrate.langevin(group=hoomd.group.all(), kT=kT, seed=seed)
    integrator.set_gamma('A', gamma=gamma)
    integrator.set_gamma('active',gamma=gamma*D_a)
    
    md.integrate.mode_standard(dt)

    hoomd.dump.gsd('data/'+name_a+'.gsd', period=int(1/(tau_u*1e-3*dt)), 
                    group=hoomd.group.all(),overwrite = True,dynamic=['momentum'])
    
    hoomd.run(tsteps=tsteps)
#%%
#u_path = input('path = ')
u_path = '/home/mike_chen/data/sed'
u_rname = input('non active gsd file name:')
# u_rname = 'relaxed'
u_aname = input('new gsd file name:')
# u_time = float(input('time(s): '))
u_time=float(2e5)
# u_alpha = float(input('percentage of active particles: '))
u_alpha = 0.05
#u_fga = float(input('mass ratio of active to passive: '))
u_fga = 1
u_Da = float(input('diameter of the active particles: '))
# u_Da = 0.8
u_theta = float(input('sedimentation angle( degree ): '))

simulate(u_path,u_rname,u_aname,u_alpha,u_time, u_fga, u_Da,u_theta)
