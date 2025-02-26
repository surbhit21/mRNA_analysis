#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:34:57 2023

@author: surbhitwagle
"""


from mRNA_SS_model import *
from lmfit import Minimizer
import matplotlib.pyplot as plt
from Utility import *
Lengths = np.array([25,50,75,100,150,200,550])
scale = 0.240
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
bin_size = 5
bins = np.arange(0, Lengths.max(), bin_size)

def FitModel(x,data,sigmas,j_r_fcator=1,pars=[]):
    if pars == []:
        fit_paramas = Parameters()
        # np.random.seed(2022)
        D_R_min = 1e-3
        D_R_max = 1.
        D_R_init = np.random.uniform(D_R_min,D_R_max)
        v_R_min = 1e-10
        v_R_max = 1e-5
        v_R_init = np.random.uniform(v_R_min,v_R_max)
        J_R_min = 1e-3*j_r_fcator
        J_R_max = 1*j_r_fcator
        J_R_init = np.random.uniform(J_R_min, J_R_max)
        # with open("./mRNAparamsweep_{}.csv".format(molecule), "a") as myfile:
        #     myfile.write("Initial param values D_r ={}, V_r = {}, Jr = {}\n".format(D_R_init, v_R_init, J_R_init))
        #     myfile.close()
        #  parameters to fit diffusion and net drift
        fit_paramas.add('D_R',D_R_init,min=D_R_min,max=D_R_max)
        fit_paramas.add('v_R',v_R_init,min=v_R_min,max=v_R_max)
        fit_paramas.add('J_R', J_R_init, min=J_R_min, max=J_R_max)
        
    else:
        fit_paramas = pars
    # breakpoint()
    resudual(fit_paramas,x,data)
    mini = Minimizer(resudual,params=fit_paramas,fcn_kws={'x':x, 'data':data})
    out2 = mini.minimize(method='leastsq')
    report_fit(out2.params)
    
    return FittedCalculation(out2.params,x,data,sigmas,mini,out2)
def resudual(paras,x=None,data=None):
    
    r_needed= GetRequiredDist(paras,x,data)
    resd = data - r_needed
    return resd #resd.flatten()



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
def GetRequiredDist(paras,x,data,off_set=1):
    
    x1,r_dist = GetParamAndModelDist(paras)
    # breakpoint()
    r_binned = BinnedSum(np.column_stack((x1,r_dist)), bins,0)[1:data.shape[0]+1,1]
    # taking the first N bins
    r_needed = r_binned
    # breakpoint()
    # normalizing with the first bin / same as soma
    r_norm_needed = r_needed
    return r_norm_needed

def GetParamAndModelDist(paras):
    # reading parameters to fit
    D_R = paras['D_R'].value
    v_R = paras['v_R'].value
    jrin = paras['J_R'].value
    # with open("test.txt", "a") as myfile:
    #     myfile.write("current param values D_r ={}, V_r = {}, Jr = {}\n".format(D_R, v_R, jrin))
    #     myfile.close()


    # return model distribution
    return RunSS(D_R,v_R,jrin)
    
def FittedCalculation(paras,x,data,sigmas,mini,out2):
    x1,r_dist = GetParamAndModelDist(paras)
    # GetParamAndModelDist
    r_needed = GetRequiredDist(paras,x,data)
    chi_squ = ChiSq(data,r_needed,sigmas)
    # delta_x = paras['dx'].value
    x_n = int(np.ceil(x[-1]/delta_x))
    return x1[0:x_n],r_dist[0:x_n],chi_squ,paras,mini,out2

def R_seq(ydata,y_fit):
    residuals = ydata- y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def plotFittedModelWithData(x1,mean_surf,mean_int,fit_surf,fit_int):
    fig,ax = plt.subplot(nrows=2,ncols=1)
    ax[0].plot