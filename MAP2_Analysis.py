#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:08:54 2021

@author: surbhitwagle
"""
import csv
import h5py
import json
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import tables
import seaborn as sb
import scipy.io
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal
from SNSPlottingWidget import SNSPlottingWidget

Lengths = np.array([2,25,50,75,100,125,150,200])
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
COLORS = [CB91_Blue,CB91_Green,CB91_Pink,CB91_Purple,CB91_Violet,CB91_Amber]
# COLORS = ["#CB48B7","#6D9F71","#2E2D4D","#337357"]
# bin_size = 5
# bins = np.arange(0, Lengths.max(), bin_size)
# print(bins)

def LoadData(fpath,cells_to_exclude,molecule,width,num_cells,bins):
    len_arr = np.zeros(Lengths.shape)
    data_dict = {}
    binned_data_dict = {}
    length_dict = {}
    data_length_wise = {}
    binned_data_length_wise = {}
    for l in Lengths:
        binned_data_length_wise[l] = []
        data_length_wise[l] = []
    for i in range(1,num_cells):
        # print(i)
        if i not in cells_to_exclude:
            data_dict["cell_{}".format(i)] = {}
            length_dict["cell_{}".format(i)] = {}
            binned_data_dict["cell_{}".format(i)] = {}
            data_folder = fpath + "cell_{}/Rois/{}/{}/".format(i,molecule,width)
            All_files = files = os.listdir(data_folder)
            for fxd,fn in enumerate(All_files):
                if fn.startswith("Dendrite") and fn.endswith(".csv"):
                    # print(fn)
                    data_dict["cell_{}".format(i)][fn.split(".")[0]] = ReadCSVFull(data_folder+fn)
                    length_dict["cell_{}".format(i)][fn.split(".")[0]]  = data_dict["cell_{}".format(i)][fn.split(".")[0]][-1,0]
                    len_arr = LenCount(len_arr,length_dict["cell_{}".format(i)][fn.split(".")[0]],"cell_{}_".format(i)+fn.split(".")[0] )
                    binned_data_dict["cell_{}".format(i)][fn.split(".")[0]] = BinnedSum(data_dict["cell_{}".format(i)][fn.split(".")[0]],bins,0,"cell_{}_".format(i)+fn.split(".")[0])[:,1]
                    for le in Lengths:
                        if length_dict["cell_{}".format(i)][fn.split(".")[0]]  >= le:
                            # breakpoint()
                            # print("appending data from cell_{} dendrite {}  {}".format(i,fn.split(".")[0],le))
                            data_length_wise[le].append(data_dict["cell_{}".format(i)][fn.split(".")[0]])
                            binned_data_length_wise[le].append(binned_data_dict["cell_{}".format(i)][fn.split(".")[0]] )
    for lex in Lengths:
        binned_data_length_wise[lex] = np.asarray(binned_data_length_wise[lex])
        data_length_wise[lex] = np.asarray(binned_data_length_wise[lex])
        
    breakpoint()
    return data_dict,binned_data_dict,data_length_wise,binned_data_length_wise,length_dict,len_arr

def BinnedSum(arr,bins,num_col=-1,name= None):
    # print(name)
    return arr
    if len(arr.shape) == 2:
        rr,cc = arr.shape
        binned_sum = np.zeros((len(bins),cc))
        digitized = bins.searchsorted(arr[:,num_col])
        # breakpoint()
        digitized[0] = digitized[1]
        for c in range(0,cc):
            binned_sum[:,c] = np.bincount(digitized, weights=arr[:,c], minlength=len(bins))
        binned_sum[:,num_col] = bins
        # breakpoint()
        return binned_sum
    else:
        print("quite not the shape",arr.shape)
        return np.zeros((len(bins),arr.shape[1]))

def LenCount(len_arr,length,name=None):
    len_arr += (length >= Lengths).astype(int)
    if length > 150:
        print("{} {}".format(name,length))
    return len_arr

def GetSumNormDistribution(data):
    breakpoint()
    sums = data.sum(axis=1)
    norm_data = np.transpose(data)/sums
    return np.transpose(norm_data)

def GetSomaNormDistribution(data,index=0):
    d_vec = data[:,index]
    norm_data = data / d_vec[:,None]
    return norm_data 
def ReadCSVFull(filename):
    csv_data = []
    # print(filename)
    with open(filename) as csvDataFile:
        # read file as csv file 
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        for row in csvReader:
    
           # add cell [0] to list of dates
           csv_data.append(row)
    # breakpoint()
    op = np.asarray(csv_data[1:])
    
    return op[:,-2:].astype(float)


def exponential(x, a, b):
    return a*np.exp(b*x)

def ExpFit(xdata,ydata,Fidx,Lidx):
    param_bounds=([-np.inf,-np.inf],[np.inf,0])
    popt, pcov = curve_fit(exponential, xdata[Fidx:Lidx], ydata[Fidx:Lidx],bounds =param_bounds)
    print(popt)
    y_fit = exponential(xdata, *popt)
    residuals = ydata[Fidx:Lidx]- y_fit[Fidx:Lidx]
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata[Fidx:Lidx]-np.mean(ydata[Fidx:Lidx]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return y_fit,r_squared


def GetSlidingWindowMean(data,window_len,mode='same'):
    try:
        conv_window = np.ones(window_len)*(1/window_len)
        sw_mean = np.convolve(data, conv_window,mode=mode)
        return sw_mean
    except:
        print("exception")
def GetSlidingWindowMeanMatrix(data,window_len,mode='same'):
    
    if len(data.shape) != 2:
        return ("data is not matrix ")
    print("here")
    op_matrix = []
    # op_matrix = np.ones((data.shape[0],op.shape[0]))
    
    for d in data:
        # breakpoint()
        op_matrix.append(GetSlidingWindowMean(d,window_len,mode))
    op_matrix = np.asarray(op_matrix)
    return op_matrix

# if __name__ == '__main__':
def plotAndSaveMap2( num_cells_dict,Excluded_cells,bins,bin_size,width,mRNA,showplot=1):
    # # mRNA = "CNIH2"
    # num_cells_dict = {"CNIH2":28,"Gria1":19,"Gria2":15}
    # Excluded_cells= {"CNIH2":[],"Gria1":[],"Gria2":[3,9,10]}
    # # width = 5
    breakpoint()
    map2_5_data,binned_data_dict,data_length_wise,binned_data_length_wise,length_data,length_arr = LoadData(mRNA+"/cells/",Excluded_cells[mRNA],"MAP2",width,num_cells_dict[mRNA],bins)
    # breakpoint()
    sum_norm_MAP2 = {}
    mean_sum_norm_MAP2 = {}
    std_sum_norm_MAP2  = {}
    sem_sum_norm_MAP2 = {}
    pwmrna = SNSPlottingWidget()
    folder = mRNA+"/MAP2-Figures/{}".format(width)
    pwmrna.CreateFolderRecursive(folder)
    # breakpoint()
    for l in Lengths[1:-1]:
        sum_norm_MAP2[l] = GetSumNormDistribution(binned_data_length_wise[l])
        mean_sum_norm_MAP2[l] = sum_norm_MAP2[l].mean(axis = 0)
        std_sum_norm_MAP2[l] = sum_norm_MAP2[l].std(axis = 0)
        sem_sum_norm_MAP2[l] = std_sum_norm_MAP2[l]/np.sqrt(binned_data_length_wise[l].shape[0])
        x = np.arange(0, l, bin_size)
        np.save(folder+"/Sum_norm_MAP2_{}_{}.npy".format(width,l),sum_norm_MAP2[l])
        np.save(folder+"/mean_MAP2_{}_{}.npy".format(width,l), mean_sum_norm_MAP2[l])
        np.save(folder+"/std_MAP2_{}_{}.npy".format(width,l), std_sum_norm_MAP2[l])
        # breakpoint()
        if showplot == 1:
            pwmrna.PlotMAP2(np.asarray([x]), np.asarray([mean_sum_norm_MAP2[l][1:x.shape[0]+1]]), np.asarray([sem_sum_norm_MAP2[l][1:x.shape[0]+1]]), ["MAP2"], "Dendritic Distance", " Fluorescence intensity (a.b.u)", "N={}".format(binned_data_length_wise[l].shape[0]), "{}/Sum_norm_MAP_length_{}_width{}".format(folder,l,width), bin_size,fit_exp=0)
            
            
def sliding_window_analysis():
    bin_size = 5
    bins = np.arange(0,200,bin_size)
    mRNA = "Gria1"
    num_cells_dict = {"CNIH2":28,"Gria1":19,"Gria2":15}
    Excluded_cells= {"CNIH2":[],"Gria1":[],"Gria2":[3,9,10]}
    width = 5#
    map2_5_data,binned_data_dict,data_length_wise,binned_data_length_wise,length_data,length_arr = LoadData(mRNA+"/cells/",Excluded_cells[mRNA],"MAP2",width,num_cells_dict[mRNA],bins)
    # breakpoint()
    for l1 in Lengths[4:5]:
        x1 = np.arange(0,l1,bin_size)
        breakpoint()
        binned_data_within_len = []
        for i in range(binned_data_length_wise[l1].shape[0]):
            binned_data_within_len.append(binned_data_length_wise[l1][i][0:x1.shape[0]])
        binned_data_within_len = np.asarray(binned_data_within_len)    
        sum_norm_MAP2 = GetSomaNormDistribution(binned_data_within_len)
        breakpoint()
        mean_sum_norm_MAP2 = sum_norm_MAP2.mean(axis = 0)
        std_sum_norm_MAP2 = sum_norm_MAP2.std(axis = 0)
        sem_sum_norm_MAP2 = std_sum_norm_MAP2/np.sqrt(binned_data_length_wise[l1].shape[0])
        x1_mean =  np.vstack((x1,mean_sum_norm_MAP2));
        x1_std = np.vstack((x1,std_sum_norm_MAP2));
        x1_sem = np.vstack((x1,sem_sum_norm_MAP2));
        np.save("../../Mo/Python-code/ampa-dynamics/Final_MAP2_{}_{}.npy".format(width,l1),x1_mean)
        np.save("../../Mo/Python-code/ampa-dynamics/Final_MAP2_std_{}_{}.npy".format(width,l1),x1_std)
        np.save("../../Mo/Python-code/ampa-dynamics/Final_MAP2_sem_{}_{}.npy".format(width,l1),x1_sem)
        op_folder = mRNA+"/MAP2-Figures/{}".format(width)
        breakpoint()
        CI = 1.96 * sem_sum_norm_MAP2
        fig, ax = plt.subplots()
        for i in range(sum_norm_MAP2.shape[0]):
            plt.plot(x1,sum_norm_MAP2[i],'gray',alpha=0.1)
        ax.plot(x1,mean_sum_norm_MAP2,'k')
        ax.fill_between(x1, (mean_sum_norm_MAP2-CI), (mean_sum_norm_MAP2+CI), color='b', alpha=.1)
        ax.set_xlabel("Dendritic distance in microns",fontsize=14)
        ax.set_ylabel("Normalized Intensity",fontsize=14)
        plt.savefig(op_folder+"MAP2_{}_{}.png".format(width,l1),dpi=300)
        # np.save(folder+"/Final_MAP2_{}_{}.npy".format(width,l),sum_norm_MAP2[l])
        # np.save(folder+"/Final_Mean_MAP2_{}_{}.npy".format(width,l), mean_sum_norm_MAP2[l])
        # np.save(folder+"/Final_std_MAP2_{}_{}.npy".format(width,l), std_sum_norm_MAP2[l])
       
        plt.show()
        plt.plot(x1,mean_sum_norm_MAP2)
        ax.fill_between(x1, mean_sum_norm_MAP2-CI, mean_sum_norm_MAP2+CI, color='r', alpha=.1)
        plt.savefig(op_folder+"MAP2_{}_{}_with_CI.png".format(width,l1),dpi=300)
        plt.show()
# sliding_window_analysis()
    
    
    