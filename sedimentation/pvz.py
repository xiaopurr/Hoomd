#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 10:51:30 2020

@author: moyuan
"""

import numpy as np
from bokeh.plotting import figure,show
from bokeh.models import Label
import bokeh
from bokeh.io import export_png
from bokeh.palettes import Inferno256
import freud
import gsd
import gsd.hoomd
import os

os.chdir('/home/moyuan/simulation/hoomd/sedimentation/')
name='sed_fg05_1HR_biglong'
#%%repositioning data
def repo20(pos):
    zero_loc= min([i[1] for i in pos])
    for i in pos:
        i[1]=i[1]-zero_loc
#%%


data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
frame=data[-1]
spacing=5
z = np.arange(0,frame.configuration.box[1],spacing)
A_rec = spacing*frame.configuration.box[0]

t=0
for frame in data:
    count = np.zeros(len(z)-1)
    pos = frame.particles.position
    repo20(pos)
    for p in range(0,len(pos)):
        for i in range(0,len(z)-1):
            if pos[p][1]>z[i] and pos[p][1]<z[i+1]:
                count[i] += 1
                break;


    phi = []
    for i in count:
        phi.append(i*np.pi*0.25/(A_rec))
        
    z_avg=np.arange(spacing/2,frame.configuration.box[1]-spacing/2,spacing)
    
    f=figure(title='phi vs z, time={0}s'.format(t),x_range=[0,0.8])
    # f.circle(z_avg, phi, fill_color='tomato',size=8, legend_label='l_c = 2 sigma')
    # f.xaxis.axis_label='z (sigma)'
    # f.yaxis.axis_label='phi (surface fraction)'
    f.circle(phi, z_avg, fill_color='tomato',size=8,legend_label='l_c = 2 sigma')
    f.xaxis.axis_label='phi (surface fraction)'
    f.yaxis.axis_label='z (sigma)'
    
    #show(f)
    export_png(f, filename='snaps/sfvz/fg05_s/s{0}.png'.format(t))
    t+=1
