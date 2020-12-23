#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 09:43:08 2020

@author: moyuan
"""

import hoomd
from hoomd import md
import numpy as np
from bokeh.plotting import figure,show
import bokeh
from bokeh.io import export_png
from bokeh.palettes import Inferno256
from bokeh.models import Span
import freud
import gsd
import gsd.hoomd
import os

import pandas as pd

hoomd.context.initialize()
seed=np.random.randint(0,9999)
os.chdir('/home/moyuan/simulation/hoomd/sedimentation/')

name='sed_fg05_active_2hr_ww_c2_c'
#unit parameters
m_u=17.05
sigma_u=2.79
e_u = 4.11
tau_u = np.sqrt(m_u*sigma_u**2/e_u)
#basic parameters
kT=1
sigma=1
gamma=23400*tau_u/m_u
time_real = 20e3#miliseconds
dt = 1e-1
tsteps = time_real/(tau_u*dt)
#force field constant
alpha = 0.56*360/(2*np.pi)
# fg=-0.116
fg=-20
epsilon=10

#make the unit box
L=1.3*sigma
box = hoomd.data.boxdim(Lx = L, Ly= L, Lz=1, xy = 0, xz = 0, yz = 0,
                        dimensions=2)
#create a blank snapshot, and specify Big and small particles
snapshot_init = hoomd.data.make_snapshot(N = 1, box = box, 
                                         particle_types=['A'])

snapshot_init.particles.position[:] = [[0,0,0]]
#replicate it 100 times
box_dim = [20,40]
snapshot_init.replicate(box_dim[0],box_dim[1],1)
N=int(snapshot_init.particles.N)
# snapshot_init.particles.typeid[0:int(N/2)] = 0;
# snapshot_init.particles.typeid[int(N/2):100] = 1;
sys = hoomd.init.read_snapshot(snapshot_init)

nl = md.nlist.cell()
wca = md.pair.lj(r_cut = sigma*np.power(2, 1/6), nlist=nl)
wca.pair_coeff.set('A', 'A', epsilon=epsilon, sigma=sigma)
# wca.pair_coeff.set('A', 'A', epsilon=10.0, sigma=1.0)
wca.set_params(mode= 'shift')

walls=md.wall.group()

walls.add_plane((0,-L*box_dim[1]/2,0), (0,1,0))
walls.add_plane((0,L*box_dim[1]/2,0),(0,-1,0))

lj_wall = md.wall.force_shifted_lj(walls, r_cut=sigma*np.power(2,1/6)/2)
lj_wall.force_coeff.set('A', epsilon=10, sigma=sigma/2)
# lj_wall.force_coeff.set('i', epsilon=0, sigma=0)

const_f = md.force.constant(fvec=(0,fg,0))

integrator = md.integrate.brownian(group=hoomd.group.all(), kT=kT, seed=seed)
integrator.set_gamma('A', gamma=gamma)
md.integrate.mode_standard(dt)

#%%
hoomd.dump.gsd('data/'+name+'.gsd', period=int(1/(tau_u*1e-3*dt)), 
                group=hoomd.group.all(),overwrite = True,dynamic=['momentum'])

hoomd.run(tsteps=tsteps)

#%%repositioning data
def repo20(pos):
    zero_loc= min([i[1] for i in pos])
    for i in pos:
        i[1]=i[1]-zero_loc
#%%colormap
name='th0067a1'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
# data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/HEX/crystal/data/HEX_rekt.gsd',mode='rb')

i= 0
theta_low = -np.pi/6
theta_high = np.pi/6
try:
    os.mkdir('snaps/sed/'+name)
except:
    
    pass;
for frame in data[::40000]:
    Lx = frame.configuration.box[0]
    Ly = frame.configuration.box[1]
    box = freud.box.Box(Lx = Lx, Ly = Ly, is2D = True)
    k = 6
    op = freud.order.Hexatic(k=k)
    
    position = []
    position_a = []
    for j in range(0,frame.particles.N):
        if frame.particles.typeid[j] == 0:
            position.append(frame.particles.position[j]-
                            np.array([0,min(frame.particles.position[:,1]),0]))
        else:
            position_a.append(frame.particles.position[j]-
                            np.array([0,min(frame.particles.position[:,1]),0]))
    # repo20(position)
    # repo20()
    
    phi_6 = op.compute(system = (box, position),neighbors={'num_neighbors':6}).particle_order
    # amp = np.array([])
    theta = np.array([])
    for p in phi_6:
        # amp = np.append(amp, np.absolute(p))
        theta = np.append(theta,(np.angle(p))/6)

    # theta= theta - min(theta)
    # theta_n = theta/max(theta)
    
    source = bokeh.models.ColumnDataSource({'x':[i[0] for i in position],
                                            'y':[i[1] for i in position],
                                            'theta':theta})
    source_a = bokeh.models.ColumnDataSource({'x':[i[0] for i in position_a],
                                              'y':[i[1] for i in position_a]})

    offset = 1
    f = figure(y_range=[0,Ly], x_range=[-Ly/2,Ly/2],
               title='snapshots')
    
    mapper = bokeh.transform.linear_cmap('theta', palette=Inferno256,low=theta_low,high=theta_high)
    
    f.circle(x = 'x',y = 'y',radius=0.5,fill_color = mapper, line_alpha = 0, 
             source = source)
    f.circle(x = 'x', y = 'y', radius = 0.5, fill_color = 'green', line_alpha=0,source = source_a)
    s = str(i)
    s=s.zfill(5)
    f.xaxis.axis_label='x (sigma)'
    f.yaxis.axis_label='z (sigma)'
    #export_png(f, filename = 'snaps/sed/'+name+'/s%s.png' %s)
    export_png(f, filename = 'snaps/sed/'+name+'/s%s.png' %s)
    # show(f)
    # output_file('/home/moyuan/activeStirring/thermal/images/s%s.html' % s)
    # save(f)
    i+=1
    print(i)

#%%area fraction retro freud

data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
i=0
for frame in data[::20]:
    var = np.var(frame.particles.position[:,1])
    spacing=[10,50,1]
    gd = freud.density.GaussianDensity(width=spacing,r_max=np.power(2,1/6)*2,sigma=var)
    gd.compute(frame)
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    gd.plot(ax=ax)
    fig.savefig(fname='snaps/density/ni_kt5e-4_g5e-4/s{0}.png'.format(i))
    i+=1

#%%time snap sfvz
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
df = pd.read_csv('phivz_relax.csv')
frame=data[-1]
spacing=2
z = np.arange(0,frame.configuration.box[1],spacing)
z_avg=np.arange(spacing/2,frame.configuration.box[1]-spacing/2,spacing)
# os.mkdir('snaps/sfvz/fg05_1HR_2')
# t=0
for t in range(0,len(data)):
    frame = data[t]

    count = np.zeros(len(z)-1)
    count_a = np.zeros(len(z)-1)
    pos_repo = frame.particles.position
    repo20(pos_repo)
    pos=[[0,0,0]]*frame.particles.N
    pos_a=[[0,0,0]]*frame.particles.N
    for p in range(0, frame.particles.N):
        if frame.particles.typeid[p]==0:
            pos[p]=(pos_repo[p])
        else:
            pos_a[p]=(pos_repo[p])

    for p in range(0,len(pos)):
        for i in range(0,len(z)-1):
            if pos[p][1]>z[i] and pos[p][1]<z[i+1]:
                count[i] += 1
                break;
    for p in range(0,len(pos)):
        for i in range(0,len(z)-1):
            if pos_a[p][1]>z[i] and pos_a[p][1]<z[i+1]:
                count_a[i] += 1
                break;
    A_rec = spacing*frame.configuration.box[0]
    phi = []
    for i in count:
        phi.append(i*np.pi*0.25/(A_rec))
    for i in range(0, len(count_a)):
        phi[i]+=count_a[i]*np.pi*((0.5*0.3)**2)/(A_rec)


    min_a_loc = min([i[1] for i in pos_a if i[1] != 0])
    min_a = Span(location=min_a_loc, dimension='height',line_dash='dashed',
                 line_color='orange',line_width=2)
    f=figure(title='phi vs z, t={0} mins'.format(t/60),y_range=[0,0.8])
    f.circle(z_avg, phi, fill_color='tomato',size=8, legend_label='active')
    f.xaxis.axis_label='z (sigma)'
    f.yaxis.axis_label='phi (surface fraction)'
    # f.circle(phi, z_avg, fill_color='tomato',size=8,legend='g=0.05')
    # f.xaxis.axis_label='phi (surface fraction)'
    # f.yaxis.axis_label='z (sigma)'
    
    # show(f)
    f.add_layout(min_a)
    
    
    #normal compare
    f.circle(df['z'],df['phi'],fill_color='white',size=8,legend_label='thermal')
    
    
    export_png(f, filename='snaps/sfvz/fg05_active_2hr/s{0}.png'.format(t))

#%%avg sfvz
name='th0083a1'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
frame=data[-1]
N=frame.particles.N
a_ind = []
Ly = frame.configuration.box[1]
i=0
for p in data[0].particles.typeid:
    if p==1:
        a_ind.append(i)
    i+=1
Na = len(a_ind)
spacing=2
z = np.arange(1,frame.configuration.box[1],spacing)
avg_n=50000
count = np.zeros(len(z)-1)
ind=[]

for frame in data[-avg_n:]:


    pos = []
    for p in range(0,Na):
        pos.append(frame.particles.position[a_ind[p]])
    repo20(pos)
    for p in range(0,len(pos)):
        for i in range(0,len(z)-1):
            if pos[p][1]>z[i] and pos[p][1]<z[i+1]:
                count[i] += 1
                break;


A_rec = spacing*frame.configuration.box[0]
phi = []
for i in count:
    phi.append(i*np.pi*0.25/(A_rec*avg_n))
        
z_avg=np.arange(spacing/2,frame.configuration.box[1]-spacing/2,spacing)


f=figure(title='phi vs z averaged over last {0} seconds {1}'.format(avg_n, name))
f.circle(z_avg, phi, fill_color='tomato',size=8)
f.xaxis.axis_label='z (sigma)'
f.yaxis.axis_label='phi (surface fraction)'
# f.circle(phi, z_avg, fill_color='tomato',size=8,legend='g=0.05')
# f.xaxis.axis_label='phi (surface fraction)'
# f.yaxis.axis_label='z (sigma)'



show(f)
export_png(f, filename='figure/dpa_sed_active{0}.png'.format(name))

df = pd.DataFrame({'z':z_avg, 'phi':phi})
df.to_csv('zl_dpa_phivz_active{0}.csv'.format(name))
#%%sf from csv
import matplotlib.pyplot as plt

plt.figure()
nl = ['th056a1','th044a1','th035a1','th025a1']
for name in nl:
    df = pd.read_csv('zl_dpa_phivz_active{0}.csv'.format(name))
    
    plt.plot(df['z'],df['phi'],'o-',label=name)

plt.legend()
plt.xlabel('z(sigma)')
plt.ylabel('surface fraction')
# f=figure(title='phi vs z averaged over last {0} seconds {1}'.format(avg_n, name))
# f.circle(df['z'], df['phi'], fill_color='tomato',size=8, legend_label='l_c = 2 sigma')
# f.xaxis.axis_label='z (sigma)'
# f.yaxis.axis_label='phi (surface fraction)'
# # f.circle(phi, z_avg, fill_color='tomato',size=8,legend='g=0.05')
# # f.xaxis.axis_label='phi (surface fraction)'
# # f.yaxis.axis_label='z (sigma)'



# # show(f)
# export_png(f, filename='figure/sed_active{0}.png'.format(name))
#%%phi v z rescaled
hg_p =  1.38e-23*300/((17.05e-15/3)*9.8*np.sin(0.252*np.pi/180))
df = pd.read_csv('phivz_relax.csv')
df['z']=df['z']*2.79e-6/hg_p
f = figure(title='phi v z averaged over 1000 seconds rescaled by hg||')
f.circle(df['z'],df['phi'])
f.xaxis.axis_label='z/hg||'
f.yaxis.axis_label='phi'
show(f)
export_png(f, filename='figure/phivz_rescaled.png')
#%%phi v z compare relaxation
df0 = pd.read_csv('phivz.csv')
dfr = pd.read_csv('phivz_relax.csv')
f = figure(title='phi v z normal and relaxed comparison')
f.circle(df0['z'],df0['phi'],fill_color='white',size=8,legend='normal')
f.circle(dfr['z'],dfr['phi'],fill_color='black',size=8,legend='relax')
show(f)
#%%Integrate surface fraction to get equation of state
from scipy import integrate
ES=np.array([])
df = pd.read_csv('phivz_active{0}.csv'.format(name))
df['z']=df['z']*2.79e-6
for i in range(1,len(df['z'])):
    ES=np.append((integrate.simps(df['phi'][:i], df['z'][:i]))/df['phi'][i],ES)
#calcualte hg||
hg_p =  1.38e-23*300/((17.05e-15/3)*9.8*np.sin(0.252*np.pi/180))

ES=ES/hg_p

#theory fit
def f(phi):
    return 1/((1-phi)**2)
phi_fit = np.linspace(min(df['phi']),max(df['phi']),int(1e5))
ES_fit = f(phi_fit)
g = figure(title='Equation of state vs phi')
g.circle(df['phi'][5:],ES[5:],legend='normal',size=8,fill_color='white')
g.line(phi_fit,ES_fit,line_color='orange',line_dash='dashed',line_width=2,
       legend = '1/(1-phi)^2')

g.xaxis.axis_label='phi'
g.yaxis.axis_label='Pi/rho kb T'

g.legend.location='top_left'
export_png(g,filename='figure/eqstate_v_phi_relax_compare{0}.png'.format(name))
show(g)
#%%Integrate surface fraction to get equation of state
th = ['056','044','035','025','0083','0067']
from scipy import integrate
# mind1 = [19.61,24.28,27.12,26.61]
dfl=[]
for theta in th:
    name ='th{0}a1'.format(theta)

    ES=np.array([])
    df = pd.read_csv('phivz_active{0}.csv'.format(name))
    df['z']=df['z']*2.79e-6
    
    
    
    for i in range(1,len(df['z'])):
        ES=np.append((integrate.simps(df['phi'][:i], df['z'][:i]))/df['phi'][i],ES)
    #calcualte hg||
    hg_p =  1.38e-23*300/((17.05e-15/3)*9.8*np.sin(0.252*np.pi/180))
    
    ES=ES/hg_p
    
    #theory fit
    def f(phi):
        return 1/((1-phi)**2)
    phi_fit = np.linspace(min(df['phi']),max(df['phi']),int(1e5))
    ES_fit = f(phi_fit)
    ES=np.append(ES, 0)
    df['ES'] = ES
    
    pdf = df
    
    pdf['z']=pdf['z']/2.79e-6

    dfl.append(pdf)

#%%
value=1.31
pdf=dfl[5]
z_fit = np.linspace(0,50,5001)
df_inter = pd.DataFrame(data={'z':z_fit})
df_inter['phi']=np.nan
df_inter['ES']=np.nan
for j in df_inter.index:
    for i in pdf.index:
        if pdf['z'][i]==df_inter['z'][j]:
            df_inter['phi'][j]=pdf['phi'][i]
            df_inter['ES'][j]=pdf['ES'][i]
inter = df_inter.interpolate(method='linear',axis=0).ffill().bfill()
inter = inter[inter['ES']!=np.inf]
import matplotlib.pyplot as plt
plt.figure()
plt.plot(inter['z'],inter['ES'])
plt.figure()
plt.plot(inter['z'],inter['phi'])
print(inter[inter['z']==value])
#%%positional order
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
frame = data[-1]

po = freud.order.Translational()
po.compute(frame)
po_list = po.particle_order;

op = np.abs(po_list)
qx = [i.real for i in po.particle_order]
qy = [i.imag for i in po.particle_order]
plt.scatter(qx,qy,color='red')
plt.show()
#%%Order parameter
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
# data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/HEX/crystal/data/HEX_rekt.gsd',mode='rb')

i= 0
theta_low = -np.pi/6
theta_high = np.pi/6
f = figure(title='average order parameter as a function of sedimentation time')

op_list=[]
for frame in data:
    opsum=0
    Lx = frame.configuration.box[0]
    Ly = frame.configuration.box[1]
    box = freud.box.Box(Lx = Lx, Ly = Ly, is2D = True)
    k = 6
    op = freud.order.Hexatic(k=k)
    
    position = []
    position_a = []
    for j in range(0,frame.particles.N):
        # if frame.particles.typeid[j] == 0:
        position.append(frame.particles.position[j])
        # else:
            # position_a.append(frame.particles.position[j])
    repo20(position)
    
    phi_6 = op.compute(system = (box, position),neighbors={'num_neighbors':6}).particle_order
    for p in phi_6:
        opsum+=np.abs(p)
    opsum=opsum/len(phi_6)
    op_list.append(opsum)

t=np.linspace(0,len(op_list),len(op_list))
f.circle(t,op_list,fill_alpha=0)
f.xaxis.axis_label='time(s)'
f.yaxis.axis_label='order parameter'
show(f)
#%%Order parameter as a function of height
name='sed_fg05_active_2hr_ww'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')

i= 0
theta_low = -np.pi/6
theta_high = np.pi/6

frame=data[-1]
spacing=2
z = np.arange(0,frame.configuration.box[1],spacing)
t=0
for frame in data:
    f = figure(title='average order parameter as a function of height, t={0}s'.format(t),
               y_range=[0,1])

    opsum=0
    Lx = frame.configuration.box[0]
    Ly = frame.configuration.box[1]
    box = freud.box.Box(Lx = Lx, Ly = Ly, is2D = True)
    k = 6
    op = freud.order.Hexatic(k=k)
    
    position = []
    position_a = []
    for j in range(0,frame.particles.N):
        if frame.particles.typeid[j] == 0:
            position.append(frame.particles.position[j])
        # else:
            # position_a.append(frame.particles.position[j])
    # repo20(position)
    
    phi_6 = op.compute(system = (box, position),neighbors={'num_neighbors':6}).particle_order
    op_list=np.zeros(len(z))
    count_list=np.zeros(len(z))
    for p in range(0,len(position)):
        for i in range(0,len(z)-1):
            if position[p][1]>z[i] and position[p][1]<z[i+1]:
                op_list[i]+=np.abs(phi_6[p])
                count_list[i]+=1
                break;
    for i in range(0,len(op_list)):
        if count_list[i]!=0:
            op_list[i]=op_list[i]/count_list[i]
    f.circle(z,op_list)
    f.xaxis.axis_label='z(sigma)'
    f.yaxis.axis_label='<|Phi_6|>'
    export_png(f,filename='snaps/opvz/fg05_active_ww/opvzt{0}.png'.format(t))
    t+=1
#%%Plotting active particle trajectories
nl = ['th056a1','th044a1','th035a1','th025a1','th0083a1','th0067a1']

for name in nl:
    
    data = gsd.hoomd.open('/home/moyuan/simulation/hoomd/sedimentation/data/{0}.gsd'.format(name))
    a_ind=[]
    i=0
    for p in data[0].particles.typeid:
        a_ind.append(i)
        i+=1
    pos = []
    for frame in data[8000:8501]:
        pos.append(frame.particles.position[a_ind[20]])
    plt.figure()
    plt.plot([i[0] for i in pos],[i[1] for i in pos],label=name)
    plt.title('500 second')
    plt.axes().set_aspect('equal')
    plt.legend()
    plt.savefig('/home/moyuan/presentations/1211/{0}_{1}.png'.format(name,500))
#%%z direction PDF
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
bin_size=2
frame = data[0]
Ly = frame.configuration.box[1]
z = np.arange(0, Ly, bin_size)
z_avg = np.arange(1,50,bin_size)
t=3600
for frame in data[::10]:

    count = np.zeros(len(z)-1)
    for p in range(0, (frame.particles.N)):
        if frame.particles.typeid[p]==1:
            pos = frame.particles.position[p][1]+Ly/2
            for i in range(0,len(z)-1):
                if pos>z[i] and pos<z[i+1]:
                    count[i]+=1
                    break;
                    
    f = figure(title='PDF for active particles z, t={0}s'.format(t),x_range=[0,0.5])
    mcount = sum(count)
    count = [i/mcount for i in count]
    f.circle(count,z_avg,size=8,fill_color='white')
    f.line(count, z_avg)
    f.xaxis.axis_label='PDF'
    f.yaxis.axis_label='z(sigma)'
    export_png(f, filename='snaps/PDF/one/s{0}.png'.format(int(t/10)))
    t+=10
    
#%%pdf v t graph
name = 'c_series/c5'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
bin_size=2
frame = data[0]
Ly = frame.configuration.box[1]
z = np.arange(0, Ly, bin_size)
z_avg = np.arange(1,50,bin_size)
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/MPZ.csv')
avg_PDF = [i for i in df['MPZ']]
for frame in data:

    count = np.zeros(len(z)-1)
    for p in range(0, (frame.particles.N)):
        if frame.particles.typeid[p]==1:
            pos = frame.particles.position[p][1]+Ly/2
            for i in range(0,len(z)-1):
                if pos>z[i] and pos<z[i+1]:
                    count[i]+=1
                    break;
                    
    scount = sum(count)
    prob = count/scount
    avg_p = sum([z_avg[i]*prob[i] for i in range(0,len(z_avg))])
    avg_PDF.append(avg_p)

df = pd.DataFrame(data={'t': np.arange(0,len(avg_PDF)),'MPZ': avg_PDF})
df.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/MPZ.csv')
#%%
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/MPZ.csv')
f = figure(title='average position of the active particles over time',
           x_range=[0,max(df['t'])],y_range=[0,max(df['MPZ'])])

f.circle(df['t'],df['MPZ'])
f.xaxis.axis_label='time(s)'
f.yaxis.axis_label='z(sigma)'
show(f)
#%%
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/MPZ.csv')
pdf = df['MPZ']
D_pdf = []
for i in range(0,len(pdf)-1,5):
    temp = pdf[i+1]-pdf[i]
    D_pdf.append(temp)
    

#%%Variance v t
name = 'c_series/c5'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')

Ly = frame.configuration.box[1]
z = np.arange(0, Ly, bin_size)
z_avg = np.arange(1,50,bin_size)
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/s2.csv')
sigma2_l=[i for i in df['s2']]
for frame in data:
    pos=[]
    for p in range(0, (frame.particles.N)):
        if frame.particles.typeid[p]==1:
            pos.append(frame.particles.position[p][1])
    pos_mean = np.mean(pos)
    sigma2=0
    for i in pos:
        sigma2+=(i-pos_mean)**2
    sigma2=sigma2/len(pos)
    sigma2_l.append(sigma2)
    
df = pd.DataFrame(data={'t': np.arange(0,len(sigma2_l)),'s2': sigma2_l})
df.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/s2.csv')
f = figure(title='average variance of the active particles z over time',
           x_range=[0,max(df['t'])],y_range=[0,max(df['s2'])])

f.circle(df['t'],df['s2'])
f.xaxis.axis_label='time(s)'
f.yaxis.axis_label='sigma**2'
show(f)
#%%
name='c3_c'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
bin_size=2
frame = data[0]
Ly = frame.configuration.box[1]
z = np.arange(0, Ly, bin_size)
z_avg = np.arange(1,50,bin_size)
count = np.zeros(len(z)-1)
for frame in data[-1000:]:

    for p in range(0, (frame.particles.N)):
        if frame.particles.typeid[p]==1:
            pos = frame.particles.position[p][1]+Ly/2
            for i in range(0,len(z)-1):
                if pos>z[i] and pos<z[i+1]:
                    count[i]+=1
                    break;
name='sed_fg05_active_2hr_ww_c2_c'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')

for frame in data:

    for p in range(0, (frame.particles.N)):
        if frame.particles.typeid[p]==1:
            pos = frame.particles.position[p][1]+Ly/2
            for i in range(0,len(z)-1):
                if pos>z[i] and pos<z[i+1]:
                    count[i]+=1
                    break;

f = figure(title='PDF for active particles z, average over t=2600s-4600s')
mcount = sum(count)
count = [i/mcount for i in count]
f.circle(np.arange(0,24,bin_size),count,size=8,fill_color='white')
f.line( np.arange(0,24,bin_size),count)
f.yaxis.axis_label='PDF'
f.xaxis.axis_label='z(sigma)'
show(f)
export_png(f, filename= '/home/moyuan/simulation/hoomd/sedimentation/figure/z_PDF_active.png')
#%%active particle trajectory
name='c_series/c5'
data = gsd.hoomd.open('data/'+name+'.gsd',mode='rb')
p=-2

ap = []
for frame in data:
    ap.append(frame.particles.position[p])
    
#%%
df = pd.DataFrame(data={'t': np.arange(0,len(ap)),'x': [i[0] for i in ap],
                        'z':[i[1] for i in ap]})
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/tra_800.csv')
# df.to_csv('/home/moyuan/simulation/hoomd/sedimentation/data/tra_799.csv')
f = figure(title='active particle #800 trajectory',x_range=[-25,25])
f.xaxis.axis_label='x(sigma)'
f.yaxis.axis_label='z(sigma)'
f.line(df['x'],df['z']-df.z.min())
export_png(f,filename='figure/active_tra_800.png')

#%%Local density vs time
df = pd.read_csv('/home/moyuan/simulation/hoomd/sedimentation/data/MPZ.csv')
