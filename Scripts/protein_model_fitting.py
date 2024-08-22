#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:34:57 2023

@author: surbhitwagle
"""

import CNIH2_SS_model
from lmfit import Minimizer
import matplotlib.pyplot as plt
import Total_AMPAR
from Utility import *

Lengths = np.array([25, 50, 75, 100, 150, 200, 550])
scale = 0.240
COLORS = ["#005f73", "#9b2226", "#CA6702", "#337357"]
# bin_size = 3
bins = np.arange(0, Lengths.max(), bin_size)


def FitModelProtein(x, data, sigmas, pars=[],molecule="CNIH2"):
    if pars == []:
        fit_paramas = Parameters()
        # np.random.seed(2024)
        D_P_min = -1.
        D_P_max = 0
        D_P_init = np.random.uniform(D_P_min, D_P_max)
        v_P_min = -5.
        v_P_max = 0.
        v_P_init = np.random.uniform(v_P_min, v_P_max)

        # bp_min = -10
        # bp_max = 0.
        # bp_init = np.random.uniform(bp_min, bp_max)
        print("initial parameters = \n"+"*"*20)
        print("D_p = {}, V_p = {} \n".format(10**D_P_init,10**v_P_init))
        print("*" * 20)
        #  parameters to fit diffusion and net drift
        fit_paramas.add('D_P', D_P_init, min=D_P_min, max=D_P_max)
        fit_paramas.add('v_P', v_P_init, min=v_P_min, max=v_P_max)
        # fit_paramas.add("bp",bp_init,min = bp_min,max= bp_max)
    else:
        fit_paramas = pars

    resudual(fit_paramas, molecule,x, data)
    mini = Minimizer(resudual, params=fit_paramas, fcn_kws={'x': x, 'data': data,'molecule':molecule})
    out2 = mini.minimize(method='leastsq')
    report_fit(out2.params)

    return FittedCalculation(out2.params, x, data, sigmas, mini, out2,molecule)


def resudual(paras, molecule, x=None, data=None):
    print(".",end="")
    r_needed = GetRequiredDist(paras, x, data,molecule,off_set=1)
    resd = data - r_needed
    # breakpoint()
    return resd  # resd.flatten()


# def GetSlidingWindowMeanMatrix(data,window_len,mode='same'):
#
#     if len(data.shape) != 2:
#         return ("data is not matrix ")
#     print("here")
#     op_matrix = []
#     # op_matrix = np.ones((data.shape[0],op.shape[0]))
#
#     for d in data:
#         # breakpoint()
#         op_matrix.append(GetSlidingWindowMean(d,window_len,mode))
#     op_matrix = np.asarray(op_matrix)
#     return op_matrix
def GetRequiredDist(paras, x, data,molecule,off_set=0):
    x1, r_dist = GetParamAndModelDist(paras,molecule)
    r_binned = BinnedSum(np.column_stack((x1, r_dist/r_dist[0])), bins, 0)[off_set:data.shape[0]+off_set, 1]
    # breakpoint()
    # taking the first N bins
    r_needed = r_binned[0:x.shape[0]]

    # normalizing with the first bin / same as soma
    r_norm_needed = r_needed / r_needed[0]
    return r_norm_needed


def GetParamAndModelDist(paras,molecule):
    # reading parameters to fit
    D_P = 10 ** paras['D_P'].value
    v_P = 10 ** paras['v_P'].value

    # bp = 10 ** paras['bp'].value

    # return model distribution
    if molecule == "CNIH2":
        return CNIH2_SS_model.RunSSProtein(D_P, v_P)
    elif molecule == "GluA2":
        return Total_AMPAR.RunSSProtein(D_P,v_P)


def FittedCalculation(paras, x, data, sigmas, mini, out2,molecule):
    x1, r_dist = GetParamAndModelDist(paras,molecule)
    # GetParamAndModelDist
    r_needed = GetRequiredDist(paras, x, data,molecule,off_set=1)
    chi_squ = RedChisq(data, r_needed, sigmas,data.shape[0]-2)
    # delta_x = paras['dx'].value
    x_n = int(np.ceil(x[-1] / delta_x))
    return x1[0:x_n+1], (r_dist / r_dist[0])[0:x_n+1], chi_squ, paras, mini, out2


def R_seq(ydata, y_fit):
    residuals = ydata - y_fit
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

