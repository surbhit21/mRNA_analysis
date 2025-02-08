#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:26:13 2022

@author: surbhitwagle
"""
import argparse
import os.path

import CNIH2_TemporalIntegration
from CNIH2_SS_model import RunModelWithFile

import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from Utility import *
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rc('font', family='Arial')


              # labels]  # ["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]

def SaveFigures(filename, ext_list=[".png", ".svg", ".pdf"], dpi=300):
    """

        function to save figures
        required arguments:
            filename
    """
    for ext in ext_list:
        plt.savefig(filename + ext, dpi=dpi)



def Plot3DMatrixIndividual(dat,fname,locn,off_set,inset_offset,step,stim_start,stim_end,title,lab,c_map,ax_label=0):
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

    asp = 'auto'
    # f_size = 14
    xlab_in, xlab, ylab = "Time (Mins)", "Time (Mins)", r"Distance from Soma ($\mu$m)"
    pos, orient, pad = 'right', 'vertical', 0.5
    fig, ax = plt.subplots(figsize=(8, 6), ncols=1)
    im1 = ax.imshow(pc_evol, aspect=asp, cmap=c_map, origin='lower',vmin = _min, vmax = _max)
    # ticks_labs = np.linspace(locn - off_set, locn + off_set, len(y_ticks))
    # print("ax1 ticks ", ticks_labs,"original ticks ",y_ticks)
    y_ticks = np.arange(min_y,max_y+1,50)  # np.arange(1, 2 * off_set + 10, 50)
    ax.set_yticks(y_ticks)
    ax.set_xmargin(0)
    ax.set_ymargin(0)
    # x_ticks =  np.arange(0, f_rames[-1]/60/60  * CNIH2_TemporalIntegration.dt, 1)
    # ax.set_xticklabels(x_ticks)
    ax.set_ylim([min_y,max_y])
    # left, bottom, width, height = [0.55, 0.55, 0.4, 0.4]
    # ax1 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)

    # print(max_y_lim)
    # im2 = ax1.imshow(inpc_evol, aspect=asp, cmap=c_map, origin='lower',vmin = _min, vmax = _max)
    # ax1.set_xticklabels(x_ticks)
    # y_ticks = np.linspace(0, 2*inset_offset ,3)  # np.arange(1, 2 * off_set + 10, 50)
    # ax1.set_yticks(y_ticks)
    # ax1.set_yticklabels(in_min_y+y_ticks)
    # ax1.set_xmargin(0)
    # ax1.set_ymargin(0)
    # x_tics = ax1.get_xticks()
    # x_labs = np.linspace(0,num_hr,len(x_tics))
    # breakpoint()
    # ax1.set_xticks([x_tics[0],x_tics[-1]])
    # ax1.set_xticklabels([0,num_hr])
    # ax1.set_yticks(y_ticks)
    # ax1.set_yticklabels(ticks_labs)
    # ax1.set_xmargin(0)
    # ax1.set_ymargin(0)
    # print("ax1 ticks ", ticks_labs,"original ticks ",y_ticks)
    # ax1.text(y = inset_offset,x = f_rames[0], s='X', color="#000000", fontsize=f_size)
    # breakpoint()
    # ax.text(y=locn, x=f_rames[0], s='X', color="#000000", fontsize=f_size)
    if ax_label == 1:
        ax.set_ylabel(ylab, fontsize=f_size)
        ax.set_xlabel(xlab, fontsize=f_size)
        # ax1.set_xlabel(xlab_in, fontsize=f_size)
        ax.set_title("[CNIH2]")
        ax.tick_params(labelsize=f_size)
    # ax.tick_params(labelsize=12)

    divider1 = make_axes_locatable(ax)
    cax1 = divider1.append_axes(pos, size='5%', pad=pad)
    cbar = plt.colorbar(im1, cax=cax1, orientation=orient)
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.set_ylabel(r'% fold change',fontsize=f_size)
    # ax.hlines(y=20, xmin=stim_start / 60, xmax=stim_end / 60, linewidth=2, color='k')
    ax.set_xticks(np.arange(timepoints[0]/60, (timepoints[-1] + 1)/60,60))
    ax.set_yticks(np.arange(min_y, max_y+ 1, 50))
    # ax1.text(x=5, y=10, s="Translation \n up-regulation", color="k", fontsize=f_size)
    # breakpoint()
    # locs, labs = ax1.yticks()
    # ax1.text(x=2, y=92.0, s="Gly stimulation", fontsize=16)  # , transform=ax.transAxes)
    # ax2 = fig.add_subplot(132)

    plt.tight_layout()
    SaveFigures(fname)
    print("figured saved")
    plt.show()




# plotStimuLocation([int(locn/dx)],"s",s_start,s_end)
# # Plot3DMatrixIndividual(data[labels[0]],locn,250,1,s_start,s_end,r"$P_c$","Pc","summer")
# breakpoint()
# PlotContour(data[labels[0]],locn,50,1,s_start,s_end,"[CNIH2]","CNIH2","RdPu",ax_label=1)
# def animateflux()

per = "100"
labels = ["CNIH2"]
if __name__ == '__main__':
    """
        parsing argumenst on mRNA and widths to be analysed. Each length is analysed seperately. mRNAs can be analysed in combinations
    """
    parser = argparse.ArgumentParser(description='mRNA analysis.py file ')
    parser.add_argument('-d', "--date_time", type=str, default="11_08_2024_16_13_52",
                        help='name of the folder to read data from')
    parser.add_argument('-se', "--end_t", type=int, default=2*60*60,
                        help='end of the stimulation')

    s_start = 0
    # reading the argumenst
    args = parser.parse_args()
    date_time = args.date_time
    s_end = args.end_t
    # s_start = 0
    # s_end = 10*60
    # l_start = 0
    # l_end = 10*60

    file_names = ["{}_{}_{}_percent.npy".format(i, date_time, per) for i in
                  labels]  # ["PC_{}_{}_percent.npy".format(date_time,per),"PS_{}_{}_percent.npy".format(date_time, per),"PSPINE_{}_{}_percent.npy".format(date_time,per)]
    input_folder = os.getcwd() + "/CNIH2-translation/CNIH2_simulations/{}/".format(date_time);
    op_folder = os.path.join(input_folder + "figures_{}".format(date_time))
    os.makedirs(op_folder, exist_ok=True)
    data = {}
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
    f_size = 30
    ss_dist, ss_model = RunModelWithFile(CNIH2_TemporalIntegration.baseline_param_file)
    L = 500.0
    dx = ss_model.dx
    x_grid = np.arange(0, L, dx)
    x_points = x_grid.shape[0]
    num_frames = int(2 * 60 / CNIH2_TemporalIntegration.dt)
    f_rames = [i for i in range(0, len(timepoints), num_frames)]
    on_time = 0
    off_time = 2 * 60 * 60

    Plot3DMatrixIndividual(data[labels[0]], "{}/Matrix_plot_{}_{}".format(op_folder,labels[0], date_time),150, 150, 10, 1, s_start, s_end, "", "CNIH2", "RdPu", ax_label=1)

