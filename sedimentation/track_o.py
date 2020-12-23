#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:19:43 2020

@author: moyuan
"""

import hoomd
from hoomd import md
import numpy as np
import gsd
import gsd.hoomd
from bokeh.plotting import figure, show
from bokeh.io import export_png
import os
m_u=17.05
sigma_u=2.79
e_u = 4.11
tau_u = np.sqrt(m_u*sigma_u**2/e_u)
#basic parameters
kT=1
sigma=1
gamma=23400*tau_u/m_u
# gamma=1
time_real = 20e3#miliseconds
dt = 1e-3
tsteps = time_real/(tau_u*dt)
epsilon=1000
seed=np.random.randint(10000,99999)
a=2
box_dim=[5,5]
hoomd.context.initialize()

# def quat_ro(q, vec):
#     qr,qi,qj,qk = q
#     mat = [[1-2*(qj**2+qk**2),2*(qi*qj-qk*qr),2*(qi*qk+qj*qr)],
#            [2*(qi*qj+qk*qr),1-2*(qi**2+qk**2),2*(qj*qk-qi*qr)],
#            [2*(qi*qk-qj*qr),2*(qj*qk+qi*qr),1-2*(qi**2+qj**2)]
#            ]
#     return np.matmul(mat, vec)

uc = hoomd.lattice.sq(a = a, type_name='A')
sys = hoomd.init.create_lattice(uc, n=box_dim)

nl=md.nlist.cell()
wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
wca.pair_coeff.set('A', 'A', epsilon=0, sigma=0)

wca.set_params(mode= 'shift')

integrator = md.integrate.langevin(group=hoomd.group.all(), kT=kT, seed=seed)
# integrator.set_gamma('A', gamma=gamma)
md.integrate.mode_standard(dt)

# #Walls at bottom and top
# walls=md.wall.group()

# walls.add_plane((0,-a*box_dim[1]/2,0), (0,1,0))
# walls.add_plane((0,a*box_dim[1]/2,0),(0,-1,0))
# walls.add_plane((-a*box_dim[0]/2,0,0), (1,0,0))
# walls.add_plane((a*box_dim[0]/2,0,0),(-1,0,0))

# lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
# lj_wall.force_coeff.set('A', epsilon=epsilon/2, sigma=sigma/2)


#active rotation
activity = []
A=30
Dr=10
for i in range(0,25):
    theta = np.random.random()*2*np.pi
    activity.append((np.cos(theta)*A,np.sin(theta)*A,0))
force = md.force.active(seed= seed, group= hoomd.group.all(),f_lst = activity,rotation_diff=Dr,
                        orientation_reverse_link=True, orientation_link=False)

#%%
hoomd.dump.gsd(filename='/home/moyuan/simulation/hoomd/sedimentation/tests/track_o_{0}.gsd'.format(seed),
               group=hoomd.group.all(),period=250,overwrite=True, dynamic=['momentum'])
hoomd.run(tsteps = tsteps)
#%%
import os
# os.system('spd-say "your program has finished"')

def calc_cos_theta(s0,s1,f0,f1):
    x1=s1[0]-s0[0]
    y1=s1[1]-s0[1]
    x2=f1[0]-f0[0]
    y2=f1[1]-f0[1]
    ang1 = np.arctan2(y1,x1)
    ang2 = np.arctan2(y2,x2)
    theta = max(ang1,ang2)-min(ang1,ang2)
    cos_th = np.cos(theta)
    return cos_th

def geo_len(p0,p1):
    x = (p1[0]-p0[0])**2;
    y = (p1[1]-p0[1])**2;
    total = np.sqrt(x+y)
    return total

cos_th = []
iteration = 0
nl = ['th056a1','th044a1','th035a1','th025a1','th0083a1','th0067a1']
for name in nl:
    tra = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/data/{0}.gsd'.format(name),
                         mode='rb')
    box_x = tra[0].configuration.box[0]
    box_y = tra[0].configuration.box[1]
    position=[]
    for snap in tra:
        raw = snap.particles.position
        for m in range(0, len(raw)):
            raw[m][0] = raw[m][0]+snap.particles.image[m][0]*box_x
            raw[m][1] = raw[m][1]+snap.particles.image[m][1]*box_y
        position.append(raw)
    #position is a list of 400 snapshots of 30 particles over 100 time unit
    cos_th.append([])
    for s_lim in np.linspace(0,10,10):
        cos_th_tmp = 0
        for p_i in range(0, len(position[0])):
            max_t0=len(position)-1
            arc_length_lim = 0;
            n=0;
            while arc_length_lim <= s_lim:
                n= n+1
                arc_length_lim += geo_len(position[max_t0-n][p_i],position[max_t0-n+1][p_i])
            max_t0 = max_t0-n
            for t0 in range(0,max_t0):
                
                    j=1
                    arc_length=geo_len(position[t0][p_i],position[t0+j][p_i])
                    while arc_length < s_lim:
                        j = j+1
                        arc_length += geo_len(position[t0+j][p_i],
                                              position[t0+j-1][p_i])
                    cos_th_tmp += calc_cos_theta(position[t0][p_i], position[t0+1][p_i],
                                                 position[t0+j-1][p_i], position[t0+j][p_i])
                    
        cos_th_tmp = cos_th_tmp/(len(position)*len(position[0]))
        cos_th[iteration].append([s_lim,cos_th_tmp])
    iteration=iteration+1
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    def func(x,a,b,c):
        return a*np.exp(-b*x)
    
    xdata = [i[0] for i in cos_th[0]]
    ydata = [i[1] for i in cos_th[0]]
    para, popt = curve_fit(func, xdata,ydata)
    x_fit = np.linspace(0,10,100)
    y_fit = func(x_fit, *para)
    plt.figure()
    plt.plot(x_fit,y_fit,label='y={0}*EXP(-{1}*x)'.format(round(para[0],2),round(para[1],2)))
    plt.plot(xdata,ydata,'x',label='data')
    plt.xlabel('arc length (sigma)')
    plt.ylabel('<cosθ>')
    plt.legend()
os.system('spd-say "your program has finished"')
#%%
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def func(x,a,b,c):
    return a*np.exp(-b*x)

xdata = [i[0] for i in cos_th[0]]
ydata = [i[1] for i in cos_th[0]]
para, popt = curve_fit(func, xdata,ydata)
x_fit = np.linspace(0,10,100)
y_fit = func(x_fit, *para)
plt.figure()
plt.plot(x_fit,y_fit,label='y={0}*EXP(-{1}*x)'.format(round(para[0],2),round(para[1],2)))
plt.plot(xdata,ydata,'x',label='data')
plt.xlabel('arc length (sigma)')
plt.ylabel('<cosθ>')
plt.legend()

##%%Calculate average active velocity
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/track_o_{0}.gsd'.format(seed),
                     mode='rb')
va = []
stdl = []
for frame in data:
    temp= np.mean(frame.particles.velocity)
    std = np.std(frame.particles.velocity)
    va.append(np.abs(temp))
    stdl.append(std/np.sqrt(frame.particles.N))
se = np.sqrt(np.sum(stdl))/len(data)
print('A = ', A)
print('Dr = ', Dr)
print('lp=',round(1/(2*para[1]),2))
print('va=',np.mean(va),'+-',se)


print('seed=',seed)
#%%make movie
name=5602
j=0
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/track_o_{0}.gsd'.format(name), mode = 'rb')
for frame in data[-100:]:
    pos = frame.particles.position
    g = figure()
    g.circle([i[0] for i in pos],[i[1] for i in pos], radius=0.5)
    export_png(g, filename = '{0}_{1}.png'.format(name, j))
    j+=1

#%%
i=0
f = figure(title='Ensemble average of Cosθ vs arc length surface fraction = %s' %i)
f.xaxis.axis_label = 'arc length'
f.yaxis.axis_label = '<cos(θ)>'
x=[]
y=[]
for j in range(0,9):
    x.append(cos_th[i][j][0])
    y.append(cos_th[i][j][1])
f.circle(x,y)
f.line(x,y)
export_png(f,filename='figure/cos_th_vs_dt_sf%s.png' % i)
    

#%%
for k in range(0,5):
    xdata = [i[0] for i in cos_th[k]]
    ydata = [i[1] for i in cos_th[k]]
    para, popt = curve_fit(func, xdata,ydata)
    g = plotting.figure(title = 'exponential fit for surface density = %s' %af_l[k])
    x_fit = np.linspace(0,10,100)
    y_fit = func(x_fit, *para)
    
    label_pl = Label(x=2,y=0.8,text = 'Lp= %.2f' % (1/(2*para[1])))
    g.add_layout(label_pl)
    g.circle(xdata,ydata)
    g.line(x_fit,y_fit,color='orange',line_dash = 'dashed')
    
    g.xaxis.axis_label = 'Arc length'
    g.yaxis.axis_label = '<cosθ>'
    export_png(g, filename='/hdd/random_walk/hoomd/figure/exp/exp_fit_sf%s.png' %af_l[k])

#%%
def quat_ro(q, vec):
    qr,qi,qj,qk = q
    mat = [[1-2*(qj**2+qk**2),2*(qi*qj-qk*qr),2*(qi*qk+qj*qr)],
           [2*(qi*qj+qk*qr),1-2*(qi**2+qk**2),2*(qj*qk-qi*qr)],
           [2*(qi*qk-qj*qr),2*(qj*qk+qi*qr),1-2*(qi**2+qj**2)]
           ]
    return np.matmul(mat, vec)

def calc_cos_theta(s0,s1,f0,f1):
    x1=s1[0]-s0[0]
    y1=s1[1]-s0[1]
    x2=f1[0]-f0[0]
    y2=f1[1]-f0[1]
    ang1 = np.arctan2(y1,x1)
    ang2 = np.arctan2(y2,x2)
    theta = max(ang1,ang2)-min(ang1,ang2)
    cos_th = np.cos(theta)
    return cos_th
################################################################
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/track_o_{0}.gsd'.format(seed),mode='rb')
N = data[0].particles.N
box = data[0].configuration.box[:3]
costh = []
for df in range(2,20):

    for ind in range(0, len(data)-df-1,10):
        for p in range(0,N):
            arc_l=0

            for i in range(1,df):
                frame0 = data[ind+i-1]
                frame1 = data[ind+i]
                arc_dis=frame1.particles.position[p]-frame0.particles.position[p]
                arc_dis = arc_dis + (frame1.particles.image[p]-frame0.particles.image[p])*box
                arc_l+= np.sqrt(arc_dis[0]**2+arc_dis[1]**2)
            
            s0 = data[ind].particles.position[p]+data[ind].particles.image[p]*box
            s1 = data[ind+1].particles.position[p]+data[ind+1].particles.image[p]*box
            f0 = data[ind+df].particles.position[p]+data[ind+df].particles.image[p]*box
            f1 = data[ind+df+1].particles.position[p]+data[ind+df+1].particles.image[p]*box
            
            cos_th_tmp = calc_cos_theta(s0,s1,f0,f1)
            
            costh.append([arc_l, cos_th_tmp])

#%%
import matplotlib.pyplot as plt
maxarc = max([np.abs(i[0]) for i in costh])

arc_disc = np.arange(0, round(maxarc,2),0.1)

temp=np.zeros(len(arc_disc)-1)
i= np.zeros(len(arc_disc)-1)
for a in range(0, len(costh)):

    for l in range(0, len(arc_disc)-1):

        if costh[a][0]>arc_disc[l] and costh[a][0]<arc_disc[l+1]:
            temp[l] += costh[a][1]
            i[l]+=1
pl = []
for j in range(0, len(temp)):
    pl.append(temp[j]/i[j])
plt.plot(arc_disc[1:],pl)
plt.plot(arc_disc[1:],[i**2 for i in pl])
#%%
data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/tests/track_o_{0}.gsd'.format(seed),mode='rb')
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
    export_png(f, filename='/home/moyuan/simulation/hoomd/sedimentation/tests/snap_{1}/track_o_{0}.png'.format(j,seed))