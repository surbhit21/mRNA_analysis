#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:35:18 2022

@author: surbhitwagle
"""

import numpy as np

def GetmRNADist(x,lam_r,lam_r2):
    
    # lam_r = (np.sqrt(V_r**2 + 4*D_r*K_r) - V_r)/(2*D_r)
    y = (lam_r/lam_r2)*np.exp(-lam_r*x)
    
    return y

