#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:53:40 2020

@author: moyuan
"""

import hoomd
from hoomd import md
import numpy as np
import os

os.chdir('/home/moyuan/simulation/hoomd/units/')

hoomd.context.initialize()
seed=np.random.randint(1000,9999)
##modes explained
#0: old parameters
#1: old parameters with converted units (kT=1)
#2: new parameters with new units
mode=4
#%%
#parameters
mass_u = 0.131
d_u = 5


T_u_0 = 298.15*2e3
T_u = 298.15
kB = 13.8

tau_u_0=np.sqrt(mass_u*d_u**2/(kB*T_u_0))
tau_u_1 = 0.0282 #seconds

time_length_real=200#seconds
gamma_real=8.38e3

lattice_dim = [10,10]
surface_fraction = 0.5
#%%
if mode==0:

    kT = 5e-4  #room temperature at 5e-4
    gamma=1    #arbitary gamma
    epsilon=10 #interaction energy
    dt = 1e-3
    tsteps=int((time_length_real/tau_u_0)/dt)
    active=0
    
if mode==1:
    kT=1
    gamma=(mass_u/tau_u_0)*tau_u_1/mass_u
    epsilon=10*2e3
    dt = 1e-3*tau_u_0/tau_u_1
    tsteps=int((time_length_real/tau_u_1)/dt)
    active=0
    
if mode==2:
    kT=1
    gamma=gamma_real*tau_u_1/mass_u
    epsilon=10*2e3
    dt = 1e-3*tau_u_0/tau_u_1
    tsteps=int((time_length_real/tau_u_1)/dt)
    active=0
    
if mode==3:
    kT=1
    gamma=gamma_real*tau_u_1/mass_u
    epsilon=10*2e3
    dt = 1e-3*tau_u_0/tau_u_1
    tsteps=int((time_length_real/tau_u_1)/dt)
    active=1
    alpha = (2)/100 #percent of active particles
    
if mode==4:
    kT=1
    gamma=gamma_real*tau_u_1/mass_u
    epsilon=10*2e3
    dt = 1e-3*tau_u_0/tau_u_1
    tsteps=int((time_length_real/tau_u_1)/dt)
    active=0
    name='sf0001'
    surface_fraction=0.001
#area of two colloids
pa = ((1/2)**2)*np.pi*2
a_x = np.sqrt((pa/surface_fraction)/np.sqrt(3))
a_y=a_x*np.sqrt(3)

uc = hoomd.lattice.hex(a=a_x, type_name='P')

system = hoomd.init.create_lattice(uc, lattice_dim)
nl = md.nlist.cell()
wca = md.pair.lj(r_cut = np.power(2, 1/6), nlist=nl)
wca.pair_coeff.set('P', 'P', epsilon=epsilon, sigma=1.0)
wca.set_params(mode= 'shift')

if active == 1:
    r_cut = np.power(2, 1/6)
    system.particles.types.add('a')
    active_id = 1 #the particle id of the active particle. Passive id is 0.
    snap = system.take_snapshot()
    N = snap.particles.N #number of total particles
    N_a = int(alpha*N) #number of active particles
    snap.particles.typeid[0:N_a]=active_id
    np.random.shuffle(snap.particles.typeid)
    
    system.restore_snapshot(snap)
    
    D_P = 1
    D_a = 0.2
    wca.pair_coeff.set('P','P',epsilon = 10, sigma = D_P, r_cut = r_cut*D_P)
    wca.pair_coeff.set('a','a',epsilon = 10, sigma = D_a, r_cut = r_cut*D_a)
    wca.pair_coeff.set('P','a',epsilon = 10, sigma = (D_P+D_a)/2, 
                       r_cut = r_cut*(D_P+D_a)/2)
    wca.set_params('shift')
    #active force
    A = 1 #amplitude factor
    angles = np.random.rand(N_a)*2*np.pi
    activity = [(np.cos(i)*A, np.sin(i)*A,0) for i in angles]
    Dr = 0.1
    md.force.active(seed=8244, group = hoomd.group.type('a'), f_lst = activity, rotation_diff=Dr,
                    orientation_link=False)


integrator = md.integrate.langevin(group=hoomd.group.all(), kT=kT, seed=seed)
integrator.set_gamma(a='P', gamma=gamma)
if active==1:
    integrator.set_gamma(a='a', gamma=(gamma*D_a))
md.integrate.mode_standard(dt)
#%%
hoomd.dump.gsd('data/{0}_long.gsd'.format(name), period=1e3, 
               group=hoomd.group.all(),overwrite = True,dynamic=['momentum'])

hoomd.run(tsteps);