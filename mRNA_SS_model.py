#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:33:54 2023

@author: surbhitwagle
"""


import json
import matplotlib.pyplot as plt
import math
import numpy as np
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
from scipy.integrate import solve_bvp
import scipy.optimize as opt       # import root-finding algorithm
# import sympy as sp                 # Python toolbox for symbolic maths
from Utility import *


fig_w, fig_h = 12, 4.5
my_fontsize = 16
my_params = {'axes.labelsize': my_fontsize,
          'axes.titlesize': my_fontsize,
          'figure.figsize': [fig_w, fig_h],
          'font.size': my_fontsize,
          'legend.fontsize': my_fontsize-4,
          'lines.markersize': 8.,
          'lines.linewidth': 2.,
          'xtick.labelsize': my_fontsize-2,
          'ytick.labelsize': my_fontsize-2}

plt.rcParams.update(my_params)

color_surf = '#005f73'
color_cyto = '#9b2226'
color_spine = '#CA6702'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
color_list = [color_surf, color_cyto, color_spine, CB91_Amber, CB91_Purple, CB91_Violet]


def CreateFolderRecursive(folder):
        Path(folder).mkdir(parents=True, exist_ok=True)

L= 100
class mRNA_model():
    """
        Class for the steady state solution of  model equations:
            \frac{\partial P_c}{dt} = D_c \frac{\partial^2 P_c}{\partial x^2} - V_p\frac{\partial P_c}{\partial x} 
            - \lambda_c P_c - \beta P_c + \alpha P_s
            
            \frac{\partial P_s}{dt} &= D_s \frac{\partial^2 P_s}{\partial x^2} - \lambda_s P_s -\alpha P_s 
            + \beta P_c -\eta P_s (\omega- P_{spine}) + \gamma P_{spine}  
             
            \frac{\partial P_{spine}}{dt} &= \eta P_s (\omega- P_{spine}) - \gamma P_{spine} 
    """

    def __init__(self, D_R, v_R, half_life_mRNA,JRin, dx):
        self.D_R = D_R
        self.v_R = v_R
        self.t_half_mRNA = half_life_mRNA
        self.k_R = np.log(2)/(self.t_half_mRNA*24*60*60); 
        self.JRin = JRin
        self.dx = dx
        self.x_grid = np.arange(0,L,dx)
        
    
    def GetSSDist(self):
        lambda_r1 = (-self.v_R + np.sqrt(self.v_R**2 + (4*self.D_R*self.k_R)))/(2*self.D_R)
        lambda_r2 = (-self.v_R - np.sqrt(self.v_R**2 + (4*self.D_R*self.k_R)))/(2*self.D_R)
        C1 = (self.JRin*np.exp(-1.*lambda_r2*L))/((self.v_R*self.D_R*lambda_r1)*(np.exp(-1.*lambda_r1*L) - np.exp(-1.*lambda_r2*L)))
        C2 = (self.JRin*np.exp(-1.*lambda_r1*L))/((self.v_R*self.D_R*lambda_r2)*(np.exp(-1.*lambda_r2*L) - np.exp(-1.*lambda_r1*L)))
        print(lambda_r1,lambda_r2,C1,C2)
    
    def fun(self, x, y):
        """
        function that defines the darivatives
        """
        r,dr = y
        return [
                dr,
                (self.v_R*dr + self.k_R*r)/self.D_R
               ]
                             
        
    def bc(self,ya,yb):
        """
            function for boundary condition values
        """
        return np.array([ya[1] - self.v_R*ya[0]/self.D_R + self.JRin/self.D_R,yb[1]- self.v_R*yb[0]/self.D_R])
    
    def SolveNumerical(self):
        y = np.zeros((2,self.x_grid.size))
        soln = solve_bvp(self.fun, self.bc, self.x_grid, y,max_nodes=1e+9, verbose=0,tol=1e-3, bc_tol=1e-8)
        # breakpoint()
        r_dist = soln.y[0]
        dr_dist = soln.y[1]
        self.x_grid = soln.x
        # breakpoint()
        # print(self.bc(soln.y[:,0],soln.y[:,1]))
        return self.x_grid,r_dist

def SaveFigures(filename,ext_list = [".png",".svg",".pdf"]):
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=300)
        

def RunSS(D_R,v_R):
    mRNA_ss = mRNA_model(D_R, v_R, 0.416, 0.001,0.012 )
    x_grid,r_dist = mRNA_ss.SolveNumerical()
    return x_grid,r_dist
    # plt.plot(mRNA_ss.x_grid,r_dist)




