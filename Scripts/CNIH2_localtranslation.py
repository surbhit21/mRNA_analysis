#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:01:02 2023

@author: surbhitwagle
"""
import matplotlib.pyplot as plt
import numpy as np
from CNIH2_SS_model import RunSSProtein
# RunSim1(0.24,1e-5,0.1,0.1)
from CNIH2_TemporalIntegration import *
from datetime import datetime


# t_max = 2000
# times = 100
# t = np.arange(0, t_max, times * dt)  # times for which to save the y values
# pc = y_init_orig[::3]
# ps = y_init_orig[1::3]
# pspine = y_init_orig[2::3]
# x_grid = np.arange(0, L, dx)
# y_init = y_init_orig





# print("BC initial = ", CheckBC(P_s_init, P_c_init))
# print("Integal condition = ", SP_model1.IntegralBC(P_s_init, P_c_init))


# breakpoint()

def PlasticityExperiment(beta,alpha,gamma,dc,vp, ds, locns, x_sdx, up_or_down_factor,a_factor,g_factor,dc_factor,vp_factor,ds_factor,step_num,op_dir,dt):

    # change_profile = np.arange(lo-x_sdx,lo+x_sdx,1)
    print(up_or_down_factor,g_factor,dc_factor,vp_factor)
    beta_final = beta.copy()
    alpha_final = alpha.copy()
    gamma_final = gamma.copy()
    dc_final = dc.copy()
    vp_final = vp.copy()
    ds_final = ds.copy()
    for lo in locns:
        beta_updated = beta.copy()
        alpha_updated = alpha.copy()
        gamma_updated = gamma.copy()
        dc_updated = dc.copy()
        vp_updated = vp.copy()
        ds_updated = ds.copy()
        if not up_or_down_factor == 1:
            beta_updated[lo-x_sdx:lo+x_sdx] *= (up_or_down_factor)
            beta_final[lo-x_sdx:lo+x_sdx] += beta_updated[lo-x_sdx:lo+x_sdx]
        if not a_factor == 1:
            alpha_updated[lo - x_sdx:lo + x_sdx] *= (a_factor)
            alpha_final[lo - x_sdx:lo + x_sdx] += alpha_updated[lo - x_sdx:lo + x_sdx]
        if not g_factor == 1:
            gamma_updated[lo] *= g_factor
            gamma_final[lo] = gamma_updated[lo]
        # breakpoint()
        if not dc_factor == 1:
            dc_updated[lo] *= (dc_factor)
            dc_final[lo] = dc_updated[lo]
        if not vp_factor == 1:
            vp_updated[lo] *= (vp_factor)
            vp_final[lo] = vp_updated[lo]
        if not ds_factor == 1:
            ds_updated[lo] *= (ds_factor)
            ds_updated[lo] = ds_updated[lo]

    # vp_updated[lo + x_sdx] *= vp_factor
    plt.plot(x_grid, beta_final / beta, label=r"$\beta$")
    plt.plot(x_grid, alpha_final / alpha, label=r"$\alpha$")
    plt.plot(x_grid, gamma_final / gamma, label=r"$\gamma$")
    # plt.plot(x_grid, dc_final / dc, label=r"$D_c$")
    # plt.plot(x_grid, vp_final / vp, label=r"$V_p$")
    # plt.plot(x_grid, ds_final / ds, label=r"$D_s$")
    plt.title("at time step {}".format(step_num))
    plt.legend()
    SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dt,step_num),dpi=300)
    plt.show()
    return beta_final,alpha_final, gamma_final, dc_final, vp_final, ds_final

def _1gaussian(x,cen1,sigma1):

    return (np.exp((-((x-cen1)/sigma1)**2)))
# def getNet_interval(locations,span):
def Getchange_profile(x_grid,locns,sigma,factor):
    change_profile = np.ones(x_grid.shape)
    for l1 in locns:
        print(l1,sigma,factor)
        change_profile += (factor*_1gaussian(x_grid,l1,sigma))
    # plt.plot(x_grid,change_profile)
    # plt.show()
        # breakpoint()
    return change_profile

def cLTPExperiment(params, factors,param_names,op_dir,step_num):
    # breakpoint()
    new_pa = params.copy()
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    plt.yscale("log")
    for pdx, p in enumerate(params):
        new_pa[pdx] = params[pdx]*factors[pdx]
        print(factors[pdx])
        if not factors[pdx] == 1:
            ax.plot(x_grid, new_pa[pdx] / params[pdx], label=param_names[pdx])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(r"Distance from soma ($\mu$m)")
    ax.set_ylabel(r"$Log_{10}$[Fold change]")
    ax.tick_params(axis='both', which='major', labelsize=18)
    # plt.xlabel("Simulation time in mins")
    # plt.ylabel("Fold change")
    plt.legend(frameon=False, fontsize=18)

    SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dt, step_num), dpi=300)
    plt.show()

    return new_pa
def ParamChangeGauss(x_grid,param, locns, sigma, factor,uod):
    bf = factor
    b_sig= sigma
    beta= param
    beta_final = beta
#     for each parameter, if the change factor is not = 1, we put a gaussian change at
#     each location in locns array with mu = location, sigma = sigma and amp = factor
    if not bf ==1:
        # beta_final = beta
        b_change_profile = Getchange_profile(x_grid,locns,b_sig,bf)
        if uod == 1:
            beta_final  = beta*b_change_profile
        elif uod == -1:
            beta_final = beta / b_change_profile
        else:
            return beta

    return beta_final

def Singledxchange(pa,locns,factor,uod,dx):
    p = pa.copy()
    for l1 in locns:
        p[int(l1/dx)]  *= ((factor)**uod)
    return p
def PlasticityExperimentGauss(x_grid,params,param_names, locns, sigmas, factors,up_or_down,step_num,op_dir,dt,dx):
    new_pa = params.copy()
    fig,ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    plt.yscale("log")
    for pdx,p in enumerate(params):
        if sigmas[pdx] == dx:
            # print("modifying ",param_names[pdx],((factors[pdx])**up_or_down[pdx]))
            # for l1 in locns:
            new_pa[pdx] = Singledxchange(params[pdx],locns,factors[pdx],up_or_down[pdx],dx)
            # for l1 in locns:
            #     print( new_pa[pdx][int(l1/dx)], new_pa[pdx][int(l1/dx)-1],params[pdx][int(l1/dx)], params[pdx][int(l1/dx)-1])
        else:
            new_pa[pdx] = ParamChangeGauss(x_grid,params[pdx],locns,sigmas[pdx],factors[pdx],up_or_down[pdx])
        if not factors[pdx] == 1:
            ax.plot(x_grid, new_pa[pdx] / params[pdx], label=param_names[pdx])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(r"Distance from soma ($\mu$m)")
    ax.set_ylabel(r"$Log_{10}$[Fold change]")
    ax.tick_params(axis='both', which='major', labelsize=18)
    # plt.xlabel("Simulation time in mins")
    # plt.ylabel("Fold change")
    plt.legend(frameon=False, fontsize=18)

    SaveFigures("{0}/fig_protocol_{1}_step_{2}".format(op_dir, dt, step_num), dpi=300)
    plt.show()
    return new_pa
def CNIH2Plasticity(dp_arr,vp_arr,kp_arr,betap_arr,mrna_arr,jpin,dx,x_grid,location):
    """
    step zero, we simulate the model in steady state for 10 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds in 50 seconds
    step Second, we step-increase the exocytosis rate (beta) to f2 folds for 10 seconds
    step last, we change back the parameters to basal level and run simulation for 30 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/../CNIH2-translation/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    sim_time = 0
    time_steps = []

    dp_factors = []
    vp_factors = []
    kp_factors = []
    betap_factors = []
    # mrna_factors = []
    # ds_factors = []

    l_scales0 = [dx,dx,dx,2]
    # breakpoint()
    """
    Step 0
    """
    dp0 = 1  # ramped up by a factor of 1
    vp0 = 1
    kp0 = 1
    bp0 = 2e+3
    # vp0 = 1
    # ds0 = 1
    P_init = protein_ss_dist
    t_step0 = 120*60  # running for t_step0 secs
    dp_factors.append(dp0)
    vp_factors.append(vp0)
    kp_factors.append(kp0)
    betap_factors.append(bp0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)

    dp_step0, vp_step0, kp_step0, betap_step0 = \
        PlasticityExperimentGauss(x_grid,
                                  [dp_arr,
                                   vp_arr,
                                   kp_arr,
                                   betap_arr],
                              [r"$D_p$", r"$V_p$", r"$K_p$", r"$\beta_p$"], location,
                              l_scales0, [dp0,vp0,kp0,bp0], [1, 1, 1, 1 ], 0, op_dir,date_time,dx)
    model_params = [dp_step0 , vp_step0, kp_step0, betap_step0, jpin,mRNA_arr, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, P_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)

    """
    Step 1
    """
    dp1 = 1  # ramped up by a factor of 1
    vp1 = 1
    kp1 = 1
    bp1 = 1
    l_scales0 = [dx, dx, dx, dx]
    # vp0 = 1
    # ds0 = 1
    P_init = data_mat[:, -1]
    t_step1 = 180 * 60  # running for t_step0 secs
    dp_factors.append(dp1)
    vp_factors.append(vp1)
    kp_factors.append(kp1)
    betap_factors.append(bp1)
    time_steps.append(t_step1)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)

    dp_step1, vp_step1, kp_step1, betap_step1 = \
        PlasticityExperimentGauss(x_grid,
                                  [dp_arr,
                                   vp_arr,
                                   kp_arr,
                                   betap_arr],
                                  [r"$D_p$", r"$V_p$", r"$K_p$", r"$\beta_p$"], location,
                                  l_scales0, [dp1, vp1, kp1, bp1], [1, 1, 1, 1 ], 1, op_dir, date_time, dx)
    # breakpoint()
    model_params = [dp_step1, vp_step1, kp_step1, betap_step1, jpin, mRNA_arr, dx]
    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1 = DynamicSimRun(model_params, t_range, t_eval, P_init, max_step=100 * dt, method='RK45')
    total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)

    saveoutput(op_dir, date_time, data_mat, total_tps, 100,baseline_param_file)


def CNIH2cLTP(dp_arr, vp_arr, kp_arr, betap_arr, mrna_arr, jpin, dx, x_grid, location):
    """
    step zero, we simulate the model in steady state for 10 secs
    step First, we step-increase the exocytosis rate (beta) to f1 folds in 50 seconds
    step Second, we step-increase the exocytosis rate (beta) to f2 folds for 10 seconds
    step last, we change back the parameters to basal level and run simulation for 30 mins
    we integrate for a total of 30 mins to see the GluA1 dynamics
    """
    now = datetime.now()
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    op_dir = os.getcwd() + "/CNIH2-translation/" + date_time
    os.makedirs(op_dir, exist_ok=True)
    print("date and time:", date_time)
    sim_time = 0
    time_steps = []

    dp_factors = []
    vp_factors = []
    kp_factors = []
    betap_factors = []
    # mrna_factors = []
    # ds_factors = []

    l_scales0 = [dx, dx, dx, 2]
    # breakpoint()
    """
    Step 0
    """
    dp0 = 1  # ramped up by a factor of 1
    vp0 = 1
    kp0 = 1
    bp0 = 2e+3
    # vp0 = 1
    # ds0 = 1
    P_init = protein_ss_dist
    t_step0 = 2*60 * 60  # running for t_step0 secs
    dp_factors.append(dp0)
    vp_factors.append(vp0)
    kp_factors.append(kp0)
    betap_factors.append(bp0)
    time_steps.append(t_step0)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)

    dp_step0, vp_step0, kp_step0, betap_step0 = \
        cLTPExperiment([dp_arr,vp_arr,kp_arr,betap_arr],
                       [dp0, vp0, kp0, bp0],
                       [r"$D_p$", r"$V_p$", r"$K_p$", r"$\beta_p$"],
                       op_dir,
                       0)
    model_params = [dp_step0, vp_step0, kp_step0, betap_step0, jpin, mRNA_arr, dx]
    t_range = [0, t_step0]
    t_eval = np.arange(0, t_step0, dt)
    soln0 = DynamicSimRun(model_params, t_range, t_eval, P_init, max_step=100 * dt, method='RK45')
    data_mat = soln0.y
    total_tps = soln0.t
    sim_time += t_step0
    print("Step 0 finished at simulation time  = ", sim_time)

    """
    Step 1
    """
    dp1 = 1  # ramped up by a factor of 1
    vp1 = 1
    kp1 = 1
    bp1 = 1
    l_scales0 = [dx, dx, dx, dx]
    # vp0 = 1
    # ds0 = 1
    P_init = data_mat[:, -1]
    t_step1 = 1*60 * 60  # running for t_step0 secs
    dp_factors.append(dp1)
    vp_factors.append(vp1)
    kp_factors.append(kp1)
    betap_factors.append(bp1)
    time_steps.append(t_step1)
    # beta_step0,alpha_Step0,gamma_step0,dc_step0, vp_step0,ds_step0 = PlasticityExperiment(beta,alpha,gamma,dc,vp,ds, lo, x_sdx, f0,a0, g0,dc0,vp0,ds0,0,op_dir,date_time)

    dp_step1, vp_step1, kp_step1, betap_step1 = \
        cLTPExperiment([dp_arr,
                                   vp_arr,
                                   kp_arr,
                                   betap_arr],
                                   [dp1, vp1, kp1, bp1],
                       [r"$D_p$", r"$V_p$", r"$K_p$", r"$\beta_p$"],
                       op_dir,
                       1)
    # breakpoint()
    model_params = [dp_step1, vp_step1, kp_step1, betap_step1, jpin, mRNA_arr, dx]
    t_range = [0, t_step1]
    t_eval = np.arange(0, t_step1, dt)
    soln1 = DynamicSimRun(model_params, t_range, t_eval, P_init, max_step=100 * dt, method='RK45')
    total_tps = np.concatenate((total_tps, (soln1.t + sim_time)))
    sim_time += t_step1
    data_mat = np.concatenate((data_mat, soln1.y), axis=1)
    print("Step 1 finished at simulation time  = ", sim_time)

    saveoutput(op_dir, date_time, data_mat, total_tps, 100, baseline_param_file)

dP_array = ss_model.D_P * np.ones(protein_ss_dist.shape)
vP_array = ss_model.v_P * np.ones(protein_ss_dist.shape)
kP_array  = ss_model.k_P * np.ones(protein_ss_dist.shape)
betaP_array = ss_model.beta_p * np.ones(protein_ss_dist.shape)
mRNA_arr = ss_model.r_dist
jpin = ss_model.JPin
x_grid = ss_model.x_grid
loc = [250]
loc_arr = [int(l) for l in loc]
# x_span = 3
# x_span_dx = int(x_span )
# breakpoint()
# CNIH2Plasticity(dP_array,
#                              vP_array,
#                              kP_array,
#                              betaP_array,
#                              mRNA_arr,
#                              jpin,
#                              ss_model.dx,
#                              x_grid,
#                              loc_arr)
CNIH2cLTP(dP_array,
                             vP_array,
                             kP_array,
                             betaP_array,
                             mRNA_arr,
                             jpin,
                             ss_model.dx,
                             x_grid,
                             loc_arr)