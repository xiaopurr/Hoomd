#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:37:27 2020

@author: moyuan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gsd
import gsd.hoomd

# en1s = 9
# for n in en1s:
# n=6

name='th056a1'

data = gsd.hoomd.open(
    '/home/moyuan/simulation/hoomd/sedimentation/data/{0}.gsd'.format(name))
Ly = data[0].configuration.box[1]
a_ind = []

for p in range(0,data[0].particles.N):
    if data[0].particles.typeid[p]==1:
        a_ind.append(p)

z = np.zeros(int(len(data)*len(a_ind)))
t = np.zeros(int(len(data)*len(a_ind)))
i=0
tt=0
minz=np.zeros(len(data))
medz = np.zeros(len(data))
maxz = np.zeros(len(data))
meanz = np.zeros(len(data))
for frame in data:
    temp=np.zeros(len(a_ind))
    j=0
    for ind in a_ind:
        z[i]=frame.particles.position[ind][1]
        temp[j]=z[i]
        t[i]=tt
        i+=1
        j+=1
    minz[tt]=min(temp)+Ly/2
    medz[tt]=np.median(temp)+Ly/2
    meanz[tt]=np.mean(temp)+Ly/2
    maxz[tt]=max(temp)+Ly/2
    tt+=1

z = z+Ly/2

dfz = pd.DataFrame({'t':t,'z':z})
dfzs = pd.DataFrame({'t':np.arange(0,len(minz)),'min':minz,'med':medz, 'max':maxz,
                     'mean':meanz})
dfz.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfz_{0}.csv'.format(name))
dfzs.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfzs_{0}.csv'.format(name))
    # log_z = np.log(z)
    # #%%
    # t_s = 20000
    # t_s = t_s*len(a_ind)
    # plt.figure()
    # plt.hist2d(t[0:t_s], log_z[0:t_s], bins=[50,50])
    # plt.xlabel('time(s)')
    # plt.ylabel('log(z)')
    # plt.colorbar()
    # plt.title('time distribution of log(z) w/ sigma_a = 0.{0} sigma'.format(n))
    # plt.show()
    #plt.savefig('/home/moyuan/simulation/hoomd/sedimentation/figure/2dhist/sigma_a_0{0}_log.png'.format(n))
#%%
name='th056a1'
data = gsd.hoomd.open(
    '/home/moyuan/simulation/hoomd/sedimentation/data/{0}.gsd'.format(name))
Ly = data[0].configuration.box[1]

dfz = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfz_{0}.csv'.format(name))
dfzs = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfzs_{0}.csv'.format(name))

minz=dfzs['min']
maxz=dfzs['max']
medz=dfzs['med']
meanz=dfzs['mean']


t_real = int(max(dfzs['t']))
N_bin=100
t_bin = int(t_real/N_bin)
plt.figure()
hist = plt.hist2d(dfz['t'], dfz['z'], bins=[range(0,t_real, t_bin),range(0,54)])

plt.xlabel('time(s)')
plt.ylabel('z')
plt.colorbar()
# plt.xlim(0,40000)
# plt.ylim(0,52)
plt.title('time distribution of z w/ {0}'.format(name))
plt.clim(0,12000)
plt.savefig('/home/moyuan/simulation/hoomd/sedimentation/figure/1204analysis/{0}.png'.format(name))

fig, ax = plt.subplots()
# ax.plot(np.arange(0,len(minz)),minz,color='r',label='min')
# ax.plot(np.arange(0,len(maxz)),maxz,color='r',label='max')
ax.plot(np.arange(0,len(medz)),medz,label = 'median')
ax.plot(np.arange(0,len(meanz)),meanz,label = 'mean')

ax.fill_between(np.arange(0,len(medz)),maxz,minz,alpha=0.3)

plt.ylim(0,Ly)
plt.xlim(0,len(minz))
plt.title('Minimum and Median active particles {0}'.format(name))
plt.xlabel('time(s)')
plt.ylabel('z(sigma)')
plt.legend(bbox_to_anchor=(1,1))
plt.savefig('/home/moyuan/simulation/hoomd/sedimentation/figure/1204analysis/active_z_0{0}.png'.format(name))
#%%
angles = ['056','044','035','025','0083','0067']
stat_f = pd.DataFrame({'angles':angles})
steady_time = 50000
ml=[]
al=[]
for th in angles:
    name='th{0}a08'.format(th)
    
    dfz = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfz_{0}.csv'.format(name))
    dfzs = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/dfzs_{0}.csv'.format(name))
    
    avgz = dfzs['mean'][-steady_time:]
    avgz = avgz.mean()
    
    minz = dfzs['min'][-steady_time:]
    minz = minz.mean()
    
    al.append(avgz)
    ml.append(minz)
    
stat_f['avgz']=al
stat_f['minz']=ml

#%%
stat_f['dif'] = stat_f['avgz']-stat_f['minz']
#%%
angles = ['056','044','035','025','0083','0067']
t=10000
spacing = 2
for th in angles:
    name='th{0}a1'.format(th)
    data = gsd.hoomd.open(
        '/home/moyuan/simulation/hoomd/sedimentation/data/{0}.gsd'.format(name))
    Ly = data[0].configuration.box[1]
    zs = np.arange(0,Ly+spacing,spacing)
    N=data[0].particles.N
    active_index=[]
    for p in range(0, N):
        if data[0].particles.typeid[p]==1:
            active_index.append(p)
    Na = len(active_index)

    for frame in data[-t:]:
        pl_tmp = []

        for p in active_index:
            pos = frame.particles.position[p][1]
            for j in range(1,len(zs)):
                if pos>zs[j-1] and pos<zs[j]:
                    count[j]+=1
        
        
        
        


        
        
#%%
th = [0.56,1,2.5,5]
mind1 = [19.61,24.28,27.12,26.61]
mind08=[19.45,24.12,26.42,27.43]

plt.scatter(th,mind1, c='r',label='sigma=1')
plt.scatter(th,mind08,c='b',label='sigma=0.8')
plt.legend()
plt.xlabel('Sed angle(degree)')
plt.ylabel('minimum z of active (sigma)')
plt.show()
#%%
en1s = range(3,10)
plt.figure()
alpha=1
c_offset = 0
for n in en1s:
    dfzs = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/zs_da0{0}.csv'.format(n))
    plt.plot(dfzs['t'],dfzs['min'],label='min Da0.{0}'.format(n),c=(1-c_offset,80/255+c_offset,80/255+c_offset,alpha))

    alpha = alpha-0.01
    c_offset = c_offset+0.09
plt.legend(bbox_to_anchor=(1,1))
plt.xlabel('time (s)')
plt.ylabel('z (sigma)')
plt.title('minimum z reached by different size active particles')

plt.figure()
alpha=1
c_offset = 0
for n in en1s:
    dfzs = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/zs_da0{0}.csv'.format(n))
    # plt.plot(dfzs['t'],dfzs['min'],label='min Da0.{0}'.format(n),c=(1-c_offset,80/255+c_offset,80/255+c_offset,alpha))
    plt.plot(dfzs['t'],dfzs['med'], label = 'median Da0.{0}'.format(n),
             c = (51/255+c_offset,204/255-c_offset,204/255-c_offset,alpha))
    alpha = alpha-0.01
    c_offset = c_offset+0.09
plt.legend(bbox_to_anchor=(1,1))
plt.xlabel('time (s)')
plt.ylabel('z (sigma)')
plt.title('median z of different size active particles')