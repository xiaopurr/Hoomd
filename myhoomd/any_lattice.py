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
    p_types = list(pas_types)
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
    rad_dict = {r: float(i) for i, r in p_dict.items()}
    
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
    lattice_type = 'sq'
    def __init__(self, Nx, Ny, ax, ay):

        self.Nx=Nx
        self.Ny=Ny
        self.ax=ax
        self.ay=ay
        
        self.Lx = Nx*ax
        self.Ly = Ny*ay
        self.N = int(Nx*Ny)
        
    def tolist(self):
        """
        

        Returns
        -------
        list
            The lattice positions in a list of ordered pairs [[x00, y00]...[xij,yij]].

        """
        return [[self.ax*i, self.ay*j] for i in range(self.Nx) for j in range(self.Ny)]
    
    def getx(self):
        return [i[0] for i in self.tolist()]
        
    def gety(self):
        return [i[1] for i in self.tolist()]
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
        return np.concatenate(([[self.a*i,ay*j] for i in range(self.Nx) for j in range(self.Ny)],
                               [[self.a*(i+1/2), ay*(j+1/2)] for i in range(self.Nx) for j in range(self.Ny)])).tolist()

    def getx(self):
        return [i[0] for i in self.tolist()]
        
    def gety(self):
        return [i[1] for i in self.tolist()]
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
        act = [[(i%pl.Nx+offset)*pl.ax, ((np.floor(i/pl.Nx)+offset)*pl.ay)] for i in insert]
        
    elif pl.lattice_type == 'hex':
        slots = range(int(pl.N/2))
        insert = np.random.choice(slots,Na)
        act = [[(i%pl.Nx)*pl.a, (np.floor(i/pl.Nx)+offset)*pl.a*np.sqrt(3)] for i in insert]
    else:
        print('lattice type unknown')
        return 'praise the sun!'
    return act