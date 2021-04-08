#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:47:45 2021

@author: moyuan
"""

import numpy as np

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
    
    Nx=None;Ny=None;
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
        
    def tolist(self):
        """
        

        Returns
        -------
        list
            The lattice positions in a list of ordered pairs [[x00, y00]...[xij,yij]].

        """
        return [[self.ax*i, self.ay*j] for i in range(self.Nx) for j in range(self.Ny)]
        
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
    Nx=None;Ny=None;
    a=None;
    Lx=None;Ly=None;
    lattice_type='hex'
    def __init__(self,Nx,Ny,a):

        self.Nx=Nx
        self.Ny=Ny
        self.a=a
        
        self.Lx = Nx*a
        self.Ly = Ny*a*np.sqrt(3)

    def tolist(self):
        """
        

        Returns
        -------
        list
            The lattice positions in a list of ordered pairs [[x00, y00]...[xij,yij]].

        """
        ay = self.a*np.sqrt(3)
        return np.concatenate(([[self.a*i,ay*j] for i in range(self.Nx) for j in range(self.Ny)],
                               [[self.a*(i+1/2), ay*(j+1/2)] for i in range(self.Nx) for j in range(self.Ny)]))

