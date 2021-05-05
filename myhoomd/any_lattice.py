#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:34:37 2021

@author: moyuan

this script attempts to creating a modular hoomd simulator with a 
functional programming approach, such that this script can be easily adopted
to simulate any porous system.
"""
import hoomd
import numpy as np
hoomd.context.initialize('')
#%%Functions
def check_lat_dim(lat):
    """
    

    Parameters
    ----------
    lat : list of lattice to check
        unknown dimension

    Returns
    -------
    list of ordered pairs
        Corrected lattice list with 3 dimensions.

    """
    if len(lat[0])>2:
        return lat
    else:
        tmp=[]
        for i in lat:
            i.append(0)
            tmp.append(i)
        print('automatically corrected lattice dimension')
        return tmp
        
def lat2snap(pl,pas_rad,act_lat):
    """
    this function converts a lattice of passive and active particles to a hoomd snapshot. 
    pas[0] is the lattice and pas[1] is the radius list, ax and ay are the basic
    lattice constants

    Parameters
    ----------
    pl : lattice obejct
        Positions of the pillar lattice defined by number of repetition in each dimension
        and the lattice constant(s)
    pas_rad : numpy array
        Ordered list of pillar radius
    act_lat : list of ordered pairs
        Initial positions of the active particles.

    Returns
    -------
    snap : Hoomd snapshot
        The snapshot that can be initialized into a hoomd simulation system.
    rad_dict : Dict
        Dictionary of the pillar radius mapping to particles type ID.

    """

    pas_types = [str(i) for i in np.unique(pas_rad)]
    
    """use a dictionary to match the radius to the particle id"""
    pid = range(len(pas_types))
    p_dict = dict(zip(pas_types,pid))
    
    """assign active particle positions"""
    Np = pl.N
    Na = len(act_lat)
    
    """make the box"""
    box = hoomd.data.boxdim(Lx=pl.Lx,Ly=pl.Ly,dimensions=2)
    
    """make the hoomd snapshot"""
    p_types = list([str(i) for i in pid])
    p_types.append('active')
    snap = hoomd.data.make_snapshot(N = Np+Na, box=box, particle_types=p_types)
    active_id = len(p_types)-1

    """parse the lattice into the snapshot with the corresponding radius types"""
    pas_lat = pl.tolist()
    pas_lat=check_lat_dim(pas_lat)
    act_lat=check_lat_dim(act_lat)
    snap.particles.position[:Np]=pas_lat
    snap.particles.position[Np:]=act_lat

    """set the typeid in the snapshot"""
    snap.particles.typeid[:Np]=[p_dict[str(i)] for i in pas_rad]
    snap.particles.typeid[Np:]=active_id
    
    """reverse map the typeid and the radius for setting interaction potential"""
    
    rad_dict = {r: round(float(i),5) for i, r in p_dict.items()}
    rad_dict['active']=0.5
    return snap, rad_dict

#%%

class pillar:
    x=0;
    y=0;
    r=0;
    def __init__(self, x, y, r):
        self.x=x
        self.y=y
        self.r=r
    
    def set_x(self,x):
        self.x = x
    def set_y(self,y):
        self.y=y
    def set_r(self,r):
        self.r=r

class sq_lattice:
    """
        

    Parameters
    ----------
    Nx : int
        Number of repetition in the x direction.
    Ny : int
        Number of repetition in the y direction.
    ax : float
        lattice constant in x direction.
    ay : float
        lattice constant in y direction.

    Returns
    -------
    None.

    """
    
    Nx=None;Ny=None;N=None
    ax=None;ay=None;
    Lx=None;Ly=None;
    N_cell_x =None;
    N_cell_y=None;
    
    max_rho = 0.22
    
    lattice_type = 'sq'
    radius=[]
    def __init__(self, Nx, Ny, ax, ay,radius=[], N_cell_x=10,N_cell_y=10):

        self.Nx=Nx
        self.Ny=Ny
        self.ax=ax
        self.ay=ay
        
        self.Lx = Nx*ax
        self.Ly = Ny*ay
        self.N = int(Nx*Ny)
        
        self.radius=radius
        self.N_cell_x = N_cell_x
        self.N_cell_y=N_cell_y
        self.w=0.25
    def tolist(self):
        """
        

        Returns
        -------
        list
            The lattice positions in a list of ordered pairs [[x00, y00]...[xij,yij]].

        """
        return [[self.ax*i-self.Lx/2, self.ay*j-self.Ly/2] for i in range(self.Nx) for j in range(self.Ny)]
    
    def getx(self):
        return [i[0] for i in self.tolist()]
        
    def gety(self):
        return [i[1] for i in self.tolist()]
    
    def generateR(self,r_type='updown',packing_fraction=0,func=None):
        Nx=self.Nx
        Ny=self.Ny
        radius=[]
        if r_type == 'updown':
            radius1 = [i/int(self.Nx) for i in range(int(self.Nx/2)) for j in range(self.Ny)]
            radius2 = list(radius1)
            radius2.reverse()
            radius=list(np.concatenate((radius1,radius2)))
        if r_type == 'updown_var':
            r_max = np.sqrt((packing_fraction*self.ax*self.ay)/np.pi)
            radius1 = [(i/int(self.Nx))*2*r_max for i in range(int(self.Nx/2)) for j in range(self.Ny)]
            radius2 = list(radius1)
            radius2.reverse()
            radius=list(np.concatenate((radius1,radius2)))
        if r_type =='updown_var_rho':
            rho_max = packing_fraction
            rho1 = [(i/int(self.Nx))*2*rho_max for i in range((int(self.Nx/2))) for j in range(self.Ny)]
            rho2 = list(rho1)
            rho2.reverse()
            rho = list(np.concatenate((rho1,rho2)))
            radius = [np.sqrt(i*self.ax*self.ay/np.pi) for i in rho]
        if r_type == 'uniform':
            pillar_area = packing_fraction*self.ax*self.ay
            radius = [np.sqrt(pillar_area/np.pi) for i in range(self.Nx) for j in range(self.Ny)]
        if r_type == 'random':
            Nx_sub = int(self.Nx/self.N_cell_x)
            Ny_sub = int(self.Ny/self.N_cell_y)
            r_cell = np.random.rand(self.N_cell_x,self.N_cell_y)/2
            radius = [r_cell[xcn][ycn] for xcn in range(self.N_cell_x) for xnc in range(Nx_sub) for ycn in range(self.N_cell_y) for ync in range(Ny_sub)]
        if r_type == 'r':
            radius=[(np.sqrt((i-Nx/2)**2+(j-Ny/2)**2))/(2*np.sqrt((Nx/2)**2+(Ny/2)**2)) for i in range(Nx) for j in range(Ny)]
        if r_type =='sin_1d':
            w = self.w
            radius=[abs(np.sin(w*(i-Nx/2)))/2 for i in range(Nx) for j in range(Ny)]
        if r_type == 'sin_2d':
            w=self.w
            radius=[abs(np.sin(w*np.sqrt((i-Nx/2)**2+(j-Ny/2)**2)))/2 for i in range(Nx) for j in range(Ny)] 
        if r_type =='1/r':
            radius=[3/(np.sqrt((i-Nx/2)**2+(j-Ny/2)**2+6**2)) for i in range(Nx) for j in range(Ny)] 
        self.radius=radius
        return radius
    # def randomizeR(self):
    #     slices = [self.radius[x:x+self.Ny] for x in range(0, len(self.radius), self.Ny)]
    #     A = self.ax*self.ay*self.Ny
    #     randomized_radius=[]
    #     for s in slices:
    #         rand_s=np.zeros(Ny)
    #         sf = sum([np.pi*(i**2) for i in s])/A
    #         for i in range(Ny):
    #             rand_s[i]=0
    
class hex_lattice:
    """
    

    Parameters
    ----------
    Nx : int
        Number of unit cells in the x direction, each unit cell has two columns.
    Ny : TYPE
        Number of unit cells in the y direction, each unit cell has two columns.
    a : TYPE
        Hexagon side length.

    Returns
    -------
    None.

    """
    Nx=None;Ny=None;N=None
    a=None;
    Lx=None;Ly=None;
    lattice_type='hex'
    radius=None
    def __init__(self,Nx,Ny,a):

        self.Nx=Nx
        self.Ny=Ny
        self.a=a
        
        self.Lx = Nx*a
        self.Ly = Ny*a*np.sqrt(3)
        self.N = int(Nx*Ny*2)
    def tolist(self):
        """
        

        Returns
        -------
        list
            The lattice positions in a list of ordered pairs [[x00, y00]...[xij,yij]].

        """
        ay = self.a*np.sqrt(3)
        return np.concatenate(([[self.a*i-self.Lx/2,ay*j-self.Ly/2] for i in range(self.Nx) for j in range(self.Ny)],
                               [[self.a*(i+1/2)-self.Lx/2, ay*(j+1/2)-self.Ly/2] for i in range(self.Nx) for j in range(self.Ny)])).tolist()

    def getx(self):
        return [i[0] for i in self.tolist()]
        
    def gety(self):
        return [i[1] for i in self.tolist()]
    def generateR(self,r_type='updown',packing_fraction=0):
        radius=[]
        if r_type == 'updown':
            radius1 = [i/(int(self.Nx)*(2**(1/6))) for i in range(int(self.Nx/2)) for j in range(self.Ny)]
            radius2 = list(radius1)
            radius2.reverse()
            radius=list(np.concatenate((radius1,radius2)))
            radius = list(np.concatenate((radius,radius)))
        if r_type == 'updown_var':
            r_max = np.sqrt((packing_fraction*(self.a**2)*np.sqrt(3))/(2*np.pi))
            radius1 = [(i/int(self.Nx))*2*r_max for i in range(int(self.Nx/2)) for j in range(self.Ny)]
            radius2 = list(radius1)
            radius2.reverse()
            radius=list(np.concatenate((radius1,radius2)))
            radius=list(np.concatenate((radius,radius)))
        if r_type =='updown_var_rho':
            rho_max = packing_fraction
            rho1 = [(i/int(self.Nx))*2*rho_max for i in range((int(self.Nx/2))) for j in range(self.Ny)]
            rho2 = list(rho1)
            rho2.reverse()
            rho = list(np.concatenate((rho1,rho2)))
            radius = [np.sqrt(i*(self.a**2)*np.sqrt(3)/(2*np.pi)) for i in rho]
            radius = list(np.concatenate((radius,radius)))
        if r_type == 'uniform':
            # pillar_area = packing_fraction*(self.a**2)*np.sqrt(3) changed 04142021 by mike
            pillar_area = packing_fraction*(self.a**2)*np.sqrt(3)/2
            radius = [np.sqrt(pillar_area/np.pi) for i in range(self.Nx) for j in range(self.Ny)]*2
        self.radius=radius
        return radius


#%%
def add_active(pl,Na):
    offset = 1/2
    """check to see if active particles would overlap initially"""
    if Na > pl.N:
        print('too many active particles')
        return None

        
    if pl.lattice_type == 'sq':
        slots = range(pl.N)
        insert = np.random.choice(slots,Na)
        act = [[(i%pl.Nx+offset)*pl.ax-pl.Lx/2, ((np.floor(i/pl.Nx)+offset)*pl.ay)-pl.Ly/2] for i in insert]
        
    elif pl.lattice_type == 'hex':
        slots = range(int(pl.N/2))
        insert = np.random.choice(slots,Na)
        act = [[(i%pl.Nx)*pl.a-pl.Lx/2, (np.floor(i/pl.Nx)+offset)*pl.a*np.sqrt(3)-pl.Ly/2] for i in insert]
    else:
        print('lattice type unknown')
        return 'praise the sun!'
    return act