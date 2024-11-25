#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
import math
import os.path

import CNIH2_TemporalIntegration
from CNIH2_SS_model import RunModelWithFile

import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from Utility import *
from matplotlib.animation import FuncAnimation, PillowWriter

date_time = "11_08_2024_16_13_52"
per = "100"
labels = ["CNIH2"]
file_names = ["{}_{}_{}_percent.npy".format(i, date_time, per) for i in
              labels]  # ["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
input_folder = os.getcwd() + "/CNIH2-translation/{}/".format(date_time);
data = {}
op_folder = os.path.join(input_folder + "figures_{}".format(date_time))
os.makedirs(op_folder, exist_ok=True)

plt.rc('font', family='Arial')
def SaveFigures(filename, ext_list=[".png", ".svg", ".pdf"], dpi=300):
    """

        function to save figures
        required arguments:
            filename
    """
    for ext in ext_list:
        plt.savefig(filename + ext, dpi=dpi)


for idx, fname in enumerate(file_names):
    data[labels[idx]] = np.load(input_folder + fname)
timepoints = np.load(input_folder + "timepoints_{0}_{1}_percent.npy".format(date_time, per))
# plt.xlabel('x');
# plt.ylabel('concentration')
# breakpoint()
dpi = 300
interval = 500
locn = 250
beta_span = 2
alpha_span = 10
fps = 25
f_size= 30
ss_dist, ss_model = RunModelWithFile(CNIH2_TemporalIntegration.baseline_param_file)
L = 500.0
dx = ss_model.dx
x_grid = np.arange(0, L, dx)
x_points = x_grid.shape[0]
num_frames = int(2*60 / CNIH2_TemporalIntegration.dt)
f_rames = [i for i in range(0, len(timepoints), num_frames)]
on_time = 0
off_time = 2*60*60
# T = 20
# dt = 0.002*10
# t_grid =  np.arange(0,T+dt,dt)
# t_points = t_grid.shape[0]
# breakpoint()
"""
Plotting the complete profile dynamics over simulated time
"""
# fig,ax = plt.subplots()
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# line, = ax.plot(x_grid,  data[labels[0]][:,0], color = COLORS_dict["shaft_i"], lw=1,label=labels[0])
#
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# plt.legend(loc="upper right")
# # plt.ylim([0,6000])
# def animate(i):
#     line.set_data(x_grid,  data[labels[0]][:,i])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     return  line, title,txt1
#
# ani = FuncAnimation(fig, animate, interval=interval, blit=True, repeat=True, frames=f_rames)
# ani.save("{}/{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()

"""
Plotting a time gif of ratio of change from initial state
"""

# fig,ax = plt.subplots()
# line, = ax.plot(x_grid,  data[labels[0]][:,0]/data[labels[0]][:,0], color = COLORS_dict["shaft_i"], lw=1,label=labels[0])
# txt1 = ax.text(x=0.3, y=0.8, s=r"stim location = {} $\mu m$".format(locn), transform=ax.transAxes)
# txt2 = ax.text(x=0.1, y=0.5, s=r"Translation \n up-regulation ={}".format("ON"), transform=ax.transAxes)
# plt.legend(loc="upper left")
# plt.plot(x_grid,  np.ones(x_grid.shape), color = 'black', lw=1,label=labels[0])
# plt.ylim([0.75,2])
# title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
#                 transform=ax.transAxes, ha="center")
# def animateRatio(i):
#     line.set_data(x_grid,  data[labels[0]][:,i]/data[labels[0]][:,0])
#     title.set_text('time = %.3f mins'  % (timepoints[i]/60))
#     txt1.set_text(r"stim location = {} $\mu m$".format(locn))
#     if timepoints[i] > on_time and timepoints[i] < off_time:
#         txt2.set_text(r"Translation \n up-regulation ={}".format("ON"))
#
#     else:
#         txt2.set_text(r"Translation \n up-regulation ={}".format("OFF"))
#     return line, title, txt1, txt2
#
# aniratio = FuncAnimation(fig, animateRatio, interval=interval, blit=True, repeat=True, frames=f_rames)
# aniratio.save("{}/ratio_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
# plt.close()


"""
Plotting a time gif of zoomed at location of stimlation loc with a surrounding area of span 
"""
fig,ax = plt.subplots(figsize=(6,4))
# loc =  50
span = 100
index_min = int((locn - span)/dx)
index_max = int((locn + span)/dx)
beta_minx = int((locn-beta_span)/dx)
beta_maxx = int((locn+beta_span)/dx)
# alpha_minx = int((locn-alpha_span)/dx)
# alpha_maxx = int((locn+alpha_span)/dx)
# range_new = [index_min:index_max]
x_grid_new = x_grid[index_min:index_max]
title = ax.text(0.5,0.95, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center",fontsize=8)
line, = ax.plot(x_grid_new,  data[labels[0]][index_min:index_max,0], color = 'blue', lw=1,label=labels[0])
ax.plot(x_grid_new,np.ones(x_grid_new.shape),"k--")
# ax.set_ylabel("Normalized [CNIH2]",fontsize=10)
# ax.set_xlabel(r"Distance from soma ($\mu$m)",fontsize=10)
txt1 = ax.text(x=0.7, y=0.8, s=r"stim location \n = {} $\mu m$".format(locn), transform=ax.transAxes,fontsize=8)
txt2 = ax.text(x=0.1, y=0.8, s=r"Translation \n up-regulation ={}".format("ON"), transform=ax.transAxes,fontsize=8)
# min_line = ax.vlines(x= x_grid[beta_minx],ymin=0.9,ymax=1.1,linestyles="dashed")
# max_line = ax.vlines(x= x_grid[beta_maxx],ymin=0.9,ymax=1.1,linestyles="dashed")
# min_line = ax.vlines(x= x_grid[alpha_minx],ymin=0.9,ymax=1.1,linestyles="-",color='k')
# max_line = ax.vlines(x= x_grid[alpha_maxx],ymin=0.9,ymax=1.1,linestyles="-",color="k")
ax.tick_params(labelsize=8)
plt.legend(loc="upper right")
plt.ylim([0.5,3])

def animateZoomed(i):
    line.set_data(x_grid_new, data[labels[0]][index_min:index_max, i]/data[labels[0]][index_min:index_max, 0])
    title.set_text('time = %.3f mins' % ((timepoints[i])/60))
    txt1.set_text(r"stim location = {} $\mu m$".format(locn))

    if timepoints[i] > on_time and timepoints[i] < off_time:
        txt2.set_text(r"Translation \n up-regulation ={}".format("ON"))

    else:
        txt2.set_text(r"Translation \n up-regulation ={}".format("OFF"))
    return line, title, txt1, txt2

anizoomed = FuncAnimation(fig, animateZoomed, interval=interval, blit=True, repeat=True, frames=f_rames)
anizoomed.save("{}/zoomed_{}.gif".format(op_folder,date_time), dpi=dpi, writer=PillowWriter(fps=fps))
plt.close()

"""
Plotting the time evolution of stimulated locations over time
"""
# print(timepoints.shape)
# breakpoint()
# patterson_1g = patterson_data_1g()
# patterson_4c = patterson_data_4c()
# tanaka_psd, tanaka_dend = tanaka_data_s3()


def plotStimuLocation(locations, sou="s",sos=0,eos=0):
    stim_or_unstim = "ST"
    if sou == "u":
        stim_or_unstim = "UN"
    fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)
    tps = timepoints / 60
    x_tics = np.arange(0,tps[-1],30)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xticks(x_tics)
    ax.set_xlabel("Time (Minutes)")
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylabel("Normalized [CNIH2]")
    ax.legend(loc="upper right", frameon=False, fontsize=f_size)
    # xlabels = [item.get_text() for item in ax.get_xticklabels()]
    for loc in locations:
        total_surf = data[labels[0]][loc, :]  # +data[labels[2]][loc, :]
        norm_total_surf = total_surf/total_surf[0]
        ax.plot(tps, norm_total_surf, color=COLORS_dict["shaft_s"], lw=1, label=loc)
        if sou == "s":
            ax.fill_betweenx(np.arange(1, norm_total_surf.max()+0.5, 0.5), sos / 60, eos / 60, alpha=0.2, color="r")

        ax.plot(tps,  np.ones(tps.shape), 'k--')
        ax.set_xticks(x_tics)
        x_loc = math.ceil(loc * dx)
        # ax.text(x=0.3,y=0.7,s=r"Uncaging at = {} $\mu m$".format(x_loc),transform=ax.transAxes,fontsize=20)
        # for axs in ax:
        # xlabels = [int(float(i)/60) for i in xlabels]
        # ax.set_xticklabels(xlabels)
    # plt.legend()
    # ax.errorbar(data_4g[:,0]+5,data_4g[:,1]/100,data_4g[:,2]/100,color = "orange")
    plt.tight_layout()
    SaveFigures("{}/{}_Time_Evolution_{}_at_location_{}".format(op_folder, stim_or_unstim, date_time, x_loc))
    plt.show()



def Plot3DMatrixIndividual(dat,locn,off_set,inset_offset,step,stim_start,stim_end,title,lab,c_map,ax_label=0):
    min_y = locn - off_set
    max_y = locn + off_set
    x_g = [i for i in range(int(min_y / dx), int(max_y / dx), int(step / dx))]

    in_min_y = locn - inset_offset
    in_max_y = locn + inset_offset
    x_g_inset = [i for i in range(int(in_min_y / dx), int(in_max_y / dx), int(step / dx))]

    fps = 60
    num_frames = int(fps / CNIH2_TemporalIntegration.dt)
    num_in_frames =  int(2*60*60/CNIH2_TemporalIntegration.dt)
    f_rames = [i for i in range(0, len(timepoints), num_frames)]
    in_f_rames = [i for i in range(0, len(timepoints[0:num_in_frames]), num_frames)]
    pc_loc_data = dat[x_g, :]
    inpc_loc_data = dat[x_g_inset, :]

    pc_init_data = pc_loc_data[:, 0]
    inpc_init_data = inpc_loc_data[:, 0]

    pc_time_loc_data = pc_loc_data[:, f_rames]
    inpc_time_loc_data = inpc_loc_data[:, in_f_rames]

    pc_evol = (pc_time_loc_data.T / pc_init_data).T
    inpc_evol = (inpc_time_loc_data.T / inpc_init_data).T
    _min, _max = 1, np.amax(pc_evol)
    # print(_min,_max)
    asp = 'auto'
    # f_size = 14
    xlab_in, xlab, ylab = "Time (Mins)", "Time (Mins)", r"Distance from Soma ($\mu$m)"
    pos, orient, pad = 'right', 'vertical', 0.5
    fig, ax = plt.subplots(figsize=(8, 6), ncols=1)
    im1 = ax.imshow(pc_evol, aspect=asp, cmap=c_map, origin='lower',vmin = _min, vmax = _max)
    # ticks_labs = np.linspace(locn - off_set, locn + off_set, len(y_ticks))
    # print("ax1 ticks ", ticks_labs,"original ticks ",y_ticks)
    # breakpoint()
    y_ticks = np.arange(min_y,max_y+1,50)  # np.arange(1, 2 * off_set + 10, 50)
    # ax.set_yticks(y_ticks)
    ax.set_xmargin(0)
    ax.set_ymargin(0)
    # x_ticks =  np.arange(0, f_rames[-1]/60/60  * CNIH2_TemporalIntegration.dt, 1)
    # ax.set_xticklabels(x_ticks)
    # ax.set_ylim([min_y,max_y])
    left, bottom, width, height = [0.55, 0.55, 0.4, 0.4]
    ax1 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)

    # print(max_y_lim)
    im2 = ax1.imshow(inpc_evol, aspect=asp, cmap=c_map, origin='lower',vmin = _min, vmax = _max)
    # ax1.set_xticklabels(x_ticks)
    y_ticks = np.linspace(0, 2*inset_offset ,3)  # np.arange(1, 2 * off_set + 10, 50)
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels(in_min_y+y_ticks)
    # ax1.set_xmargin(0)
    # ax1.set_ymargin(0)
    x_tics = ax1.get_xticks()
    # x_labs = np.linspace(0,num_hr,len(x_tics))
    # breakpoint()
    # ax1.set_xticks([x_tics[0],x_tics[-1]])
    # ax1.set_xticklabels([0,num_hr])
    # ax1.set_yticks(y_ticks)
    # ax1.set_yticklabels(ticks_labs)
    # ax1.set_xmargin(0)
    # ax1.set_ymargin(0)
    # print("ax1 ticks ", ticks_labs,"original ticks ",y_ticks)
    ax1.text(y = inset_offset,x = f_rames[0], s='X', color="#000000", fontsize=f_size)
    # breakpoint()
    ax.text(y=locn, x=f_rames[0], s='X', color="#000000", fontsize=f_size)
    if ax_label == 1:
        ax.set_ylabel(ylab, fontsize=f_size)
        ax.set_xlabel(xlab, fontsize=f_size)
        ax1.set_xlabel(xlab_in, fontsize=f_size)
        # ax.set_title("[CNIH2]")
        ax.tick_params(labelsize=f_size)
    # ax.tick_params(labelsize=12)

    divider1 = make_axes_locatable(ax)
    cax1 = divider1.append_axes(pos, size='5%', pad=pad)
    cbar = plt.colorbar(im1, cax=cax1, orientation=orient)
    cbar.ax.yaxis.set_tick_params(labelsize=f_size)
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.set_ylabel(r'Fold $\Delta$ in concentration',fontsize=f_size)
    ax.hlines(y=20, xmin=stim_start / 60, xmax=stim_end / 60, linewidth=2, color='k')
    ax.set_xticks(np.arange(timepoints[0]/60, (timepoints[-1] + 1)/60,60))
    ax.set_yticks(np.arange(min_y, max_y+ 1, 50))
    # ax1.text(x=5, y=10, s="Translation \n up-regulation", color="k", fontsize=f_size)
    # breakpoint()
    # locs, labs = ax1.yticks()
    # ax1.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    # ax2 = fig.add_subplot(132)

    plt.tight_layout()
    SaveFigures("{}/Matrix_plot_{}_{}".format(op_folder,lab, date_time))
    print("figured saved")
    plt.show()


s_start = 0
s_end = 2*60*60
# plotStimuLocation([int(locn/dx)],"s",s_start,s_end)
# # Plot3DMatrixIndividual(data[labels[0]],locn,250,1,s_start,s_end,r"$P_c$","Pc","summer")
Plot3DMatrixIndividual(data[labels[0]],150,150,10,1,s_start,s_end,"","CNIH2","RdPu",ax_label=1)

# breakpoint()
# PlotContour(data[labels[0]],locn,50,1,s_start,s_end,"[CNIH2]","CNIH2","RdPu",ax_label=1)
# def animateflux()