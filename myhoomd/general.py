#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:47:45 2021

@author: moyuan
"""

import numpy as np

import hoomd
from hoomd import md
def savedata(l,name):
    import pandas as pd
    df = pd.DataFrame(l)
    df.to_csv(name)


def make_mov(data_name,filename,skip=25):

    from bokeh.plotting import figure
    from bokeh.io import export_png
    import gsd.hoomd
    
    data=gsd.hoomd.open(data_name)
    rad_dict = np.load(data_name[:len(data_name)-4]+'_radius_dict.npy', allow_pickle='TRUE').item()
    s=0
    frame0=data[0]
    act_ind = [i for i,x in enumerate(frame0.particles.typeid) if x==max(frame0.particles.typeid)]
    rad_list=[]
    for r in rad_dict:
        if r != 'active':
            rad_list.append([i for i,x in enumerate(frame0.particles.typeid) if x == r])
    x_lim = frame0.configuration.box[0]
    for frame in data[::skip]:
        f = figure(x_range=[-x_lim/2,x_lim/2], y_range=[-x_lim/2,x_lim/2],match_aspect=True)
        key=0
        for r in rad_list:
            f.circle(frame.particles.position[r,0],frame.particles.position[r,1],radius=rad_dict[key])
            key+=1
        f.circle(frame.particles.position[act_ind,0],frame.particles.position[act_ind,1],radius=0.5,color='red')
        export_png(f, filename=filename+'{0}.png'.format(str(s).zfill(6)))
        s+=1


def simulate(snap, r_dict, save_file, va=4, t=5, dt=0.1, fps=25,null_lattice=False):
    
    seed=np.random.randint(0,9999)
    kT=1
    v0=va*(1/(1/(0.669e-3))) #4 micron / second
    gamma=4130
    fr=gamma*(4/3)*(1/4)

    
    Dr=1/fr
    fa = v0*gamma
    tau = 0.669
    t_real = t*3600e3 #millisecond = 5 hr
    dt=dt
    tsteps=int(t_real/(tau*dt))
    sigma=1
    # epsilon=10
    if null_lattice:
        epsilon=0
    else:
        epsilon=10
    # typeid = snap.particles.typeid
    print(snap.box)
    sys = hoomd.init.read_snapshot(snap)
    
    nl = md.nlist.cell()
    wca = md.pair.lj(r_cut = sigma, nlist=nl)
    #set all interaction potential between each type
    print('the radii corresponding to the typeid is: ')
    print(r_dict)
    for i in sys.particles.types:
        if i != 'active':
            s = r_dict[int(i)]+r_dict['active']
            # wca.pair_coeff.set(i, 'active', epsilon=epsilon, sigma=s,r_cut=s*(2**(1/6)))
            wca.pair_coeff.set(i, 'active', epsilon=epsilon, sigma=s,r_cut=s)
            for j in sys.particles.types:
                if j != 'active':
                    wca.pair_coeff.set(i,j, epsilon=0, sigma=0)
    wca.pair_coeff.set('active', 'active', epsilon=0, sigma=0) #changed after two tests
    wca.set_params(mode= 'shift')
    
    #set activity
    a_ind = max(snap.particles.typeid)
    Na = len([i for i in snap.particles.typeid if i == a_ind])
    print(Na)
    activity = []
    A=fa
    theta = np.random.rand(Na)*2*np.pi
    activity=[(np.cos(i)*A, np.sin(i)*A,0) for i in theta]
    force = md.force.active(seed= seed, group= hoomd.group.type('active'),f_lst = activity,rotation_diff=Dr,
                            orientation_reverse_link=True, orientation_link=False)
    
    integrator = md.integrate.brownian(group=hoomd.group.all(), kT=kT, seed=seed)
    integrator.set_gamma('active', gamma=gamma)
    for i in sys.particles.types:
        if i != 'active':
            integrator.set_gamma(i, gamma=gamma*1e40)
    
    
    md.integrate.mode_standard(dt)
    
    hoomd.dump.gsd(filename=save_file,
                   period=int(1e3/(fps*tau*dt)), group=hoomd.group.all(), 
                   overwrite=True, dynamic=['momentum'])
    metadata_file_name = save_file[:len(save_file)-4]+'_radius_dict.npy'
    np.save(metadata_file_name, r_dict)
    hoomd.run(tsteps = tsteps)