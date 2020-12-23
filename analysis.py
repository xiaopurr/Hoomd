#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 15:30:33 2020

@author: moyuan
"""

import gsd
import gsd.hoomd
import numpy as np
import os
name = 'sf01'
os.chdir('/home/moyuan/simulation/hoomd/units/')
mode=2


tra = gsd.hoomd.open('data/{0}.gsd'.format(name), mode='rb')
tra= tra[:317028]
box_x = tra[-1].configuration.box[0]
box_y = tra[-1].configuration.box[1]
#number of particles
p_num = tra[0].particles.N;

def f_len(a):
    l=0
    for i in a:
       l=l+i**2
    return l

msds=[]
for l in range(0,3):
    total_steps = len(tra)
    max_step_size = int(0.25* total_steps)
    msd = np.array([])
    counter = 0
    #step_skip = 1250
    step_skip = 1585 #
    p_N = 30;
    p_list = np.random.choice(np.arange(0,p_num),size=p_N)
    t0_skip = 200
    for step_size in range(0, max_step_size,step_skip):
        d_avg_t0 = 0;
        active_count=0
        for p in [int(p) for p in p_list]:
            if tra[0].particles.typeid[p] == 0: #check if the particle is passive
                d_sum = 0;
                for t0 in range(0,total_steps - step_size, t0_skip):
                    t0_comp = [tra[t0].particles.image[p][0]*box_x, # putting All images into a list in the beginning
                               tra[t0].particles.image[p][1]*box_y, 0]
                    
                    t0_pos = tra[t0].particles.position[p]+t0_comp
                    #step forward in time and grad position of particle #p
                    snap_comp = [tra[t0+step_size].particles.image[p][0]*box_x, 
                                 tra[t0+step_size].particles.image[p][1]*box_y, 0];
                    snap_pos = tra[t0+step_size].particles.position[p]+snap_comp
                    displacement = snap_pos - t0_pos
        
                    distance = f_len(displacement)
        
                    d_sum = d_sum + distance
                #divide to get the average over all t0
                d_avg_t0 = d_avg_t0 + d_sum*t0_skip/(total_steps - step_size + 1)
            else:
                active_count +=1
        completion_percentage = float(step_size/max_step_size)
        print(completion_percentage/3+l/3)
        msd = np.append(msd, d_avg_t0/(p_N-active_count))

    msds.append(msd)

import pickle
with open('data/MSD/msd_{0}_long.txt'.format(name),'w+b') as handle:
    pickle.dump(msds, handle)
    

#%% Ensemble Avg

import pickle
from scipy import stats
from bokeh.plotting import figure, show
from bokeh.io import export_png
from bokeh.models import Label, Legend
from bokeh.palettes import Category20_13
t_real = 200
max_tstep=t_real/4
names = ['001','01','02','03','04','05','06','07','08','09','1']
f = figure(title='Stoke-Einstein dilute limit correction')
palette=0
legend_item=[]
for name in names:
    with open('data/MSD/msd_sf0{0}_long.txt'.format(name),'r+b') as handle:
        msds=pickle.load(handle)
        

    msd_ensemble = np.zeros(len(msds[0]))
    for each in msds:
        msd_ensemble = msd_ensemble+each
    
    msd_ensemble = msd_ensemble/len(msds)
    
    #unit conversion to real units
    msd_ensemble = msd_ensemble*25#diameter=5, d^2=5^2=25
    #
    t = np.linspace(0,max_tstep,len(msd_ensemble))
    linreg = stats.linregress(t, msd_ensemble)
    t_fit = np.linspace(0,max(t),int(1e5))
    msd_fit = t_fit*linreg[0]+linreg[1]
    ps = f.circle(t,msd_ensemble,fill_color=Category20_13[palette],line_alpha=0)
    pf = f.line(t_fit,msd_fit,color=Category20_13[palette],line_dash='dashed')
    
    legend_item.append(('sf={0}%'.format(name),[ps]))
    legend_item.append(('sf{0}% - fit'.format(name), [pf]))
    
    f.xaxis.axis_label='time(s)'
    f.yaxis.axis_label='MSD (micron^2)'
    f.legend.location='top_left'
    
    # slope_label = Label(x=t[1],y=msd_ensemble[-1],text = 'slope = {0} (micron^2/s)'.format(round(linreg[0],3)))
    # f.add_layout(slope_label)
    # show(f)
    palette+=1
    
    print('sf{0}, slope={1}'.format(name, linreg[0]))
    
legend = Legend(items=legend_item, location='center')
f.add_layout(legend,'right')
export_png(f, filename='sf_long.png')

#%%Plotting D v sf
sf = [0.1e-2,1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,0.1]
D = np.array([1.98345,1.89890,1.75650,1.8221,1.84822,1.73382,1.89113,1.58419,1.7805,1.58629,1.57495])
D = D/4
from bokeh.plotting import figure, show
from bokeh.models import Span
f = figure(title='Diffusion coeifficient vs surface fraction')
f.circle(sf,D)
hline = Span(location=0.494, dimension='width',line_color='orange',line_dash='dashed',
             line_width=2)
f.xaxis.axis_label='surface fraction'
f.yaxis.axis_label='D (micron^2/s)'
f.renderers.extend([hline])
show(f)