#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:34:37 2021

@author: moyuan

this package
"""
import gsd.hoomd
import numpy as np

def load_data(file):
    data=gsd.hoomd.open(file)
    rad_dict = np.load(file[:len(file)-4]+'_radius_dict.npy', allow_pickle='TRUE').item()
    # print(rad_dict)
    return data, rad_dict

def binning_test(file, file_fps=25, bin_size = [1/25,1/5,1/2,1,2,5,10,20,40]):
    # data, rad_dict = load_data(file)
    out=[]
    for b in bin_size:
        veff=calculate_effective_velocity_for_uniform_pillars(file, fps=int(file_fps*b))
        out.append([b,veff])
    return out


def calc_D(file,fps=25,tau_r = 1):
    from scipy.stats import linregress
    s=[120,170,220,370,320,370,420]
    s = [tau_r*i for i in s]
    msd=[]
    for i in s:
        msd.append(calc_MSD(file,fps=fps,st=100,skip=i))
    output=linregress(s,msd)
    return output[0], output[-1]


def calculate_effective_velocity_for_uniform_pillars(file,fps=25,st=100,skip=40):
    data, rad_dict = load_data(file)
    data = data[int(st*fps):]
    frame0 = data[0]
    actInd = [i for i,x in enumerate(frame0.particles.typeid) if x == max(frame0.particles.typeid)]
    box = frame0.configuration.box[:3]
    
    output= []
    for f in range(0,len(data)-int(skip*fps),int(skip*fps)):
        f0 = data[f]
        f1 = data[f+int(fps*skip)]
        
        p0 = f0.particles.position[actInd]
        p1 = f1.particles.position[actInd]
        
        i0 = f0.particles.image[actInd]
        i1 = f1.particles.image[actInd]
        
        realP0 = np.array(p0+i0*box)
        realP1 = np.array(p1+i1*box)
        
        veff = realP1-realP0
        
        veff = [np.sqrt(i[0]**2+i[1]**2)/skip for i in veff]
        
        output.append(np.mean(veff))
        
    return np.mean(output)


def calc_MSD(file,fps=25,st=100,skip=40):
    data, rad_dict = load_data(file)
    data = data[int(st*fps):]
    frame0 = data[0]
    actInd = [i for i,x in enumerate(frame0.particles.typeid) if x == max(frame0.particles.typeid)]
    box = frame0.configuration.box[:3]
    
    output= []
    for f in range(0,len(data)-int(skip*fps),fps):
        f0 = data[f]
        f1 = data[f+int(fps*skip)]
        
        p0 = f0.particles.position[actInd]
        p1 = f1.particles.position[actInd]
        
        i0 = f0.particles.image[actInd]
        i1 = f1.particles.image[actInd]
        
        realP0 = np.array(p0+i0*box)
        realP1 = np.array(p1+i1*box)
        
        veff = realP1-realP0
        
        veff = [i[0]**2+i[1]**2 for i in veff]
        
        output.append(np.mean(veff))
        
    return np.mean(output)


# def v_profile(rad_dict, axis=0):
#     if axis == 0:
#         return [x,v]
#     elif axis == 1:
#         return [y,v]
#     else:
#         return 'Error: invalid axis, 0 for x and 1 for y'
# def pdf(data, normalized_by_free_space=True):
    
#     return [x,phi]