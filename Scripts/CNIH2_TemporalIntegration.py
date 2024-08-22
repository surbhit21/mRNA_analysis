from CNIH2_SS_model import RunSSProtein, RunModelWithFile

import json
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.integrate import solve_ivp
import shutil
from Utility import model_L

dt = 1e-2
baseline_param_file = "./CNIH2SSModelParamsTemporal.json"
protein_ss_dist, ss_model = RunModelWithFile(baseline_param_file)

def Phi(t,t0):
    if t >= t0:
        return 1
    return 0
def local_translation(mRNa_dist,trans_rate):
    return trans_rate*mRNa_dist

def Degradation(protein_dist,deg_rate):
    return protein_dist*deg_rate

def get_temporal_profile(t_c,t):
    t0, a1, a2, b1, b2, T_max = t_c
    return Phi(t,t0)*(t-t0)*(a1*np.exp(-a2*(t-t0)) - b1*np.exp(-b2*(t-t0)))/T_max
    # return t_profile

def get_spatial_profile(x_c,x_grid):
    beta_b, x_0, sigma = x_c
    x_prof = beta_b*np.exp(-((x_grid-x_0)**2)/sigma)
    return x_prof

def get_burst_profile(betap,t_c,t,x_c,x_grid):
    t_profile = get_temporal_profile(t_c,t)
    x_profile = get_spatial_profile(x_c,x_grid)
    breakpoint()
    return x_profile

def CNIH2TimeDynamics(t, y, Dp,vP,kP,betaP,Jpin,mRNA, dx):
    """
   Differential equations for the 1-D coupled diffusion-advection-reaction-trapping equations.

   The ODEs are derived using the method of lines.
   """
    # The vectors pc, ps and pspineare interleaved in y.  We define
    # views of pc,ps and pspine by slicing y.
    if t % 2 == 0:
        print("time = ", t, )
    # print(".",)
    P = y
    # dydt is the return value of this function.
    dydt = np.empty_like(y)
    # updated after calculations on 22 Aug 2024
    dydt[0] =  (Dp[0]/dx**2 - vP[0]/dx)* P[1] - \
               (Dp[0] / dx**2) * P[0] - \
               Degradation(P[0],kP[0]) + \
               (Jpin / dx) +\
               local_translation(mRNA[0],betaP[0])

    dydt[1:-1] = Dp[1:-1]* np.diff(P,2) / dx**2 - \
                 vP[1:-1] * np.diff(P,1) [1:] / dx - \
                 Degradation(P[1:-1],kP[1:-1]) + \
                 local_translation(mRNA[1:-1],betaP[1:-1])

    dydt[-1] = (Dp[-1] / dx**2)*P[-2] - \
               (Dp[-1]/dx**2 + vP[-1]**2/Dp[-1] - vP[-1]/dx) * P[-1] - \
               Degradation(P[-1],kP[-1]) + \
               local_translation(mRNA[-1], betaP[-1])



    return dydt

def DynamicSimRun(model_params, t_range, t_eval, y_init, method="LSODA", dense_op=True, vectorize=True, lband=2,
                  uband=2, rtol=1e-5, atol=1e-11, max_step=100):
    soln = solve_ivp(CNIH2TimeDynamics, t_range, y_init, args=(model_params), \
                     method=method, dense_output=dense_op, vectorize=vectorize, lband=lband, uband=uband, rtol=rtol,
                     atol=atol, max_step=max_step, t_eval=t_eval)
    # breakpoint()
    return soln

def SaveFigures(filename,ext_list = [".png",".svg",".pdf"],dpi=300):
    """
        function to save figures
        required arguments:
            filename
    """
    print("saving as" ,filename)
    for ext in ext_list:
        plt.savefig(filename+ext,dpi=dpi)
# def Bleach(y, lo, x_sdx):
#     pc_l, ps_l, psp_l = GetProteindata(y)
#     ps_l[lo - x_span_dx:lo + x_span_dx] = 0
#     psp_l[lo - x_span_dx:lo + x_span_dx] = 0
#     return np.vstack((np.vstack((pc_l, ps_l)), psp_l)).T.flatten()

def saveoutput(op_dir, date_time, soln, tps, percent, orig_params_file):
    p_cnih2= soln
    with open('{0}/CNIH2_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, p_cnih2)
    f.close()
    t_points = tps
    with open('{0}/timepoints_{1}_{2}_percent.npy'.format(op_dir, date_time, percent), 'wb') as f:
        np.save(f, t_points)
    f.close()
    shutil.copyfile(orig_params_file, '{0}/baseline_parameters_{1}.json'.format(op_dir, date_time))


def GetProteindata(soln):
    pc_t = soln[::3]
    ps_t = soln[1::3]
    p_sp_t = soln[2::3]
    return pc_t, ps_t, p_sp_t

def savesimsettings(num_steps, time_steps, locations, protocol_file,
                    dp_factors=[],vp_factors=[],kp_factors=[],betap_factors=[],jp_factors=[]):
    # assert num_steps == len(time_steps) - 1
    # assert num_steps == len(factors) - 1
    protocol_details = {}
    protocol_details["num_steps"] = num_steps
    protocol_details["time_steps"] = time_steps
    protocol_details["locations"] = locations
    # protocol_details["x_span"] = x_span

    protocol_details["dp_factors"] = dp_factors
    protocol_details["vp_factors"] = vp_factors
    protocol_details["betap_factors"] = betap_factors
    protocol_details["kp_factors"] = kp_factors
    protocol_details["jp_factors"] = jp_factors

    with open (protocol_file, 'a') as fp:
        json.dump(protocol_details,fp)
    fp.close()
