#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:39:47 2022

@author: surbhitwagle
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import json
import os
from pathlib import Path
from scipy.stats import ks_2samp, kruskal
import seaborn as sb
from scipy.optimize import curve_fit
folder = os.path.abspath(os.getcwd())+"/CNIH2/"
Channel_names = {0: "Dapi", 1: "Cnih2", 2: "MAP2", 3: "CAMKII"}
scale1 = 2.41  # pixels/uM
scale2 = 7.587  # pixels/uM
channels = [1, 3]
dend_data = {}
soma_data = {}
dend_data_meta = {}
soma_data_meta = {}
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
COLORS = [CB91_Blue, CB91_Green, CB91_Pink,
          CB91_Purple, CB91_Violet, CB91_Amber]
Lengths = np.asarray([25, 50, 100, 150, 200, 250])
widths = [5.0, 10, 15.0]
Analysis_dend_width = 0


def ReadFiles(dir):
    cells = range(1, 29)
    # print(cells)
    total_d_count = 0
    for i in channels:
        soma_data['channel_'+str(i)] = {}
        dend_data['channel_'+str(i)] = {}
        for cell in cells:
            # print(cell)
            file = folder+"cell_"+str(cell)+"/"+str(i)+"/"+str(
                widths[Analysis_dend_width])+"/dend_puncta_channel{0}_{1}.json".format(i, widths[Analysis_dend_width])
            file2 = folder+"cell_"+str(cell)+"/"+str(2)+"/"+str(
                widths[Analysis_dend_width])+"/dend_puncta_channel{0}_{1}.json".format(2, widths[Analysis_dend_width])
            # print(file)
            if Path(file).is_file():
                dend_data_meta['cell_'+str(cell)], soma_data_meta['cell_'+str(
                    cell)] = LoadDend(folder+"cell_"+str(cell)+"/", scale1)
                total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                f_d = open(file)
                dend_data['channel_'+str(i)]['cell_' +
                                             str(cell)] = json.load(f_d)
                f_d.close()
                f_s = open(folder+"cell_"+str(cell)+"/"+str(i)+"/"+str(
                    widths[Analysis_dend_width])+"/soma_puncta_channel{0}_{1}.json".format(i, widths[Analysis_dend_width]))
                # print(f_d,f_s)
                soma_data['channel_'+str(i)]['cell_' +
                                             str(cell)] = json.load(f_s)
                f_s.close()
            if Path(file2).is_file():
                dend_data_meta['cell_'+str(cell)], soma_data_meta['cell_'+str(
                    cell)] = LoadDend(folder+"cell_"+str(cell)+"/", scale2)
                total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                f_d = open(file2)
                dend_data['channel_1']['cell_'+str(cell)] = json.load(f_d)
                f_d.close()
                f_s = open(folder+"cell_"+str(cell)+"/"+str(2)+"/"+str(
                    widths[Analysis_dend_width])+"/soma_puncta_channel{0}_{1}.json".format(2, widths[Analysis_dend_width]))
                # print(f_d,f_s)
                soma_data['channel_1']['cell_'+str(cell)] = json.load(f_s)
                f_s.close()
    soma_sigmas,soma_cell_sigmas = GetSigmaList(soma_data)
    dend_sigmas,dend_cell_sigmas = GetSigmaList(dend_data)
    breakpoint()


def GetSigmaList(data_dict):
    all_sigma_list = {}
    all_cell_sigma_list = {}
    for c in data_dict.keys():
        all_sigma_list[c] = []
        all_cell_sigma_list[c] = {}
        for cell in data_dict[c].keys():
            all_cell_sigma_list[c][cell] = []
            for unit in data_dict[c][cell].keys():
                for stat in data_dict[c][cell][unit]:
                    sigma = stat[2]/np.sqrt(2)
                    all_cell_sigma_list[c][cell].append(sigma)
                    all_sigma_list[c].append(sigma)
    return all_sigma_list,all_cell_sigma_list

def PlotSigmas(self,fractions,lab,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        # ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        pos = np.arange(1,3,1)
        # breakpoint()
        # for i in range(0,fractions.shape[0]):
        bp1 = ax.boxplot(np.transpose(fractions[0:2]),widths = 0.25,positions=pos,showmeans=True,meanline=True,labels=lab[0:2],showfliers=False)
        bp2 = ax.boxplot(np.transpose(fractions[2:]),widths = 0.25,positions=3+pos,showmeans=True,meanline=True,labels=lab[2:],showfliers=False)
        self.setBoxColors(bp1,COLORS[0])
        self.setBoxColors(bp2,COLORS[1])
        # ax.bars(cells,dend,label=lab2,color=CB91_Blue)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        # plt.legend(prop={'size': fsize})
        # plt.show()
        fig.tight_layout()
        folder = "."
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=150)
            plt.savefig(file_name+".svg",dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
def LoadDend(Dir, scale):
    """
    Input:
            Dir (String)   : Super directory we are looking at
    Output:
            DendArr  (numpy array) : location of dendrite

    Function:
            Read npy files to obtain dendrites
    """

    DendArr = {}
    # DendArr_names = []
    soma = []
    for x in os.listdir(Dir):
        if('AllDendrites' in x):
            DendArr = np.load(Dir+x,allow_pickle=True)
        elif("Soma" in x) and x.endswith(".npy"):
            soma.append(np.load(Dir+x))
    return DendArr, soma


def GetDendlength(xys, scale):
    length = 0
    for idx in range(0, xys.shape[0]-1):
        length += Distance(xys[idx], xys[idx+1])
    return length/scale


def Distance(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def GetPolygonArea(x, y, scale):
    return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))/(scale**2)


if __name__ == '__main__':
    ReadFiles(folder)
