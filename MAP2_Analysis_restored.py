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

# # from PlotBinned import *
# def MAP2_Analysis(data_file_path,cells_to_exclude):
#     molecule = "MAP2"
#     dendritic_lens = np.array([25,50,75,100,125]); # dendritic lengths for which data is available
#     analysis_for_lengths = np.array([0,0,0,1,0]); # dendritic lengths for which the analysis has to be done each value should be either 0 or 1
#     count = np.array([0,0,0,0,0])
#     step_size = 25;
#     bin_size = 5; # for binned analysis
#     # scale = 1/4.8177; #um/pixel
#     data_dict = {};
#     for lens in dendritic_lens:
#         data_dict[int(lens/step_size)-1] = {}
#     # data = np.array('data')
#     # print(data_file['MAP2_data'].shape[1])
#     for cells in data_file['MAP2_data'][0,3:]:
#         for dendrite in cells[0,:]:
#             dendrite_length = dendrite[0,2].shape[0]*scale;
#             for lens in dendritic_lens:
#                 if dendrite_length > lens:
#                     # print(dendrite[0,2].shape)
#                     # print(int(lens/step_size))
#                     cindex= int(lens/step_size) - 1
#                     data_dict[int(lens/step_size)-1][count[cindex]] =  dendrite[0,2];
#                     count[cindex]  += 1;
#     # print(data_dict[0].keys())
    
#     for i in range(0,len(analysis_for_lengths)):
#         file_name = "MAP2_distribution_for_%i"%dendritic_lens[i];
#         lab =  r'Mean $\pm$ S.E.M'
#         y_lab = "Normalized Intensity"
#         x_lab = r'Dendritic distance (in $\mu$M)'
#         title = 'MAP2 intensity'
#         if analysis_for_lengths[i] == 1:
#             # print(count[i])
#             print(file_name)
#             # PlotBinned(data_dict[i],dendritic_lens[i], scale, count[i], bin_size, file_name,title,lab,x_lab,y_lab, molecule,1)
Lengths = np.array([2,25,50,75,100,150,200])
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
COLORS = [CB91_Blue,CB91_Green,CB91_Pink,CB91_Purple,CB91_Violet,CB91_Amber]
# COLORS = ["#CB48B7","#6D9F71","#2E2D4D","#337357"]
bin_size = 2
bins = np.arange(0, Lengths.max(), bin_size)
print(bins)

def LoadData(fpath,cells_to_exclude,molecule,width,num_cells):
    len_arr = np.zeros(Lengths.shape)
    data_dict = {}
    binned_data_dict = {}
    length_dict = {}
    data_length_wise = {}
    for l in Lengths:
        data_length_wise[l] = []
    for i in range(1,num_cells):
        if i not in cells_to_exclude:
            data_dict["cell_{}".format(i)] = {}
            length_dict["cell_{}".format(i)] = {}
            binned_data_dict["cell_{}".format(i)] = {}
            data_folder = fpath + "cell_{}/Rois/{}/{}/".format(i,molecule,width)
            All_files = files = os.listdir(data_folder)
            for fxd,fn in enumerate(All_files):
                if fn.startswith("Dendrite") and fn.endswith(".csv"):
                    data_dict["cell_{}".format(i)][fn.split(".")[0]] = ReadCSVFull(data_folder+fn)
                    length_dict["cell_{}".format(i)][fn.split(".")[0]]  = data_dict["cell_{}".format(i)][fn.split(".")[0]][-1,0]
                    len_arr = LenCount(len_arr,length_dict["cell_{}".format(i)][fn.split(".")[0]],"cell_{}_".format(i)+fn.split(".")[0] )
                    binned_data_dict["cell_{}".format(i)][fn.split(".")[0]] = BinnedSum(data_dict["cell_{}".format(i)][fn.split(".")[0]],bins,0,"cell_{}_".format(i)+fn.split(".")[0])[:,1]
                    for l in Lengths:
                        if length_dict["cell_{}".format(i)][fn.split(".")[0]]  >= l:
                            data_length_wise[l].append(binned_data_dict["cell_{}".format(i)][fn.split(".")[0]] )
    for l in Lengths:
        data_length_wise[l] = np.asarray(data_length_wise[l])
    return data_dict,binned_data_dict,data_length_wise,length_dict,len_arr

def BinnedSum(arr,bins,num_col=-1,name= None):
    # print(name)
    if len(arr.shape) == 2:
        rr,cc = arr.shape
        binned_sum = np.zeros((len(bins),cc))
        digitized = bins.searchsorted(arr[:,num_col])
        # breakpoint()
        digitized[0] = digitized[1]
        for c in range(0,cc):
            binned_sum[:,c] = np.bincount(digitized, weights=arr[:,c], minlength=len(bins))
        binned_sum[:,num_col] = bins
        return binned_sum
    else:
        print("quite not the shape",arr.shape)
        return np.zeros((len(bins),arr.shape[1]))

def LenCount(len_arr,length,name=None):
    len_arr += (length >= Lengths).astype(int)
    return len_arr

def GetSumNormDistribution(data):
    sums = data.sum(axis=1)
    norm_data = np.transpose(data)/sums
    return np.transpose(norm_data)

def ReadCSVFull(filename):
    csv_data = []
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
class PlottingWidgetmRNA():
    def CreateFolderRecursive(self,folder):
        Path(folder).mkdir(parents=True, exist_ok=True)
        
    def PlottwoStats(self,x_data,y_data,lab,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.scatter(x_data,y_data,label=lab)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        plt.legend(prop={'size': fsize})
        # plt.show()
        folder = "."
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=150)
            plt.savefig(file_name+".svg",dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()

    def PlotFourStats(self,x_data1,y_data1,x_data2,y_data2,lab1,lab2,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        ax.scatter(x_data1,y_data1,label=lab1)
        ax.scatter(x_data2,y_data2,label=lab2)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        plt.legend(prop={'size': fsize})
        # plt.show()
        folder = "."
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=150)
            plt.savefig(file_name+".svg",dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    
    def PlotBinnedStats(self,xs,means,stds,labs,xlab,ylab,title_string,file_name,bin_size,width=10,height=8,fsize=16,save_it = 1,fit_exp =1,auto_close=1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0')
        if xs.shape[0] == means.shape[0]:
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i]+bin_size/2,means[i],stds[i],label=labs[i],color=COLORS[i],marker='o')
            ax.set_xlabel(xlab,fontsize=fsize)
            ax.set_ylabel(ylab,fontsize=fsize)
            plt.title(title_string,fontsize=fsize)
            
            folder = "."
            
            if fit_exp == 1:
                for i in range(xs.shape[0]):
                    yi_fit, ri_squared = ExpFit(xs[i],means[i],1,-1)
                    ax.plot(xs[i]+bin_size/2,yi_fit,'^--',c=COLORS[i],label=labs[i]+"-fit, $r^2$ =%0.2f" %(ri_squared))
            plt.legend(prop={'size': fsize})
            # plt.show()
            fig.tight_layout()
            if save_it == 1:
                plt.savefig(file_name+".png",dpi=150)
                plt.savefig(file_name+".svg",dpi=150)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
            if auto_close == 1:
                plt.close()
        else:
            print("Not same number of xs and means")
    
    def PlotCellFraction(self,fractions,lab,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
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
    
    def setBoxColors(self,bp,c):
        setp(bp['boxes'], color=c)
        setp(bp['caps'], color=c)
        setp(bp['caps'], color=c)
        setp(bp['whiskers'], color=c)
        setp(bp['whiskers'], color=c)
        setp(bp['fliers'], color=c)
        setp(bp['fliers'], color=c)
        setp(bp['medians'], color=c)
    
    def PlotBinnedStatsMulti(self,x,data_1,lab1,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        
        
        for d in data_1:
            plt.plot(x,d,color=CB91_Blue,alpha=0.2)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        plt.legend(prop={'size': fsize})
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
    def ViolinPlotStats(self,data1,data2,lab1,lab2,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        
        # ax.set_axis_style([lab1,lab2])
        labels=[lab1+"\n N=%d cells"%(data1.shape[0]),lab2+"\n N=%d cells"%(data2.shape[0])]
        x = np.arange(len(labels))
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        xs = np.zeros(data1.shape)
        violin_parts = ax.violinplot([data1,data2],positions=x,showextrema=True,showmeans = True,showmedians = True,points=data1.shape[0])
        ax.scatter(xs,data1,color=CB91_Blue,alpha=0.6)
        x2s = np.zeros(data2.shape)+(x[1])
        ax.scatter(x2s,data2,color=CB91_Blue,alpha=0.6)
        # print(violin_parts['bodies'])
        for pc in violin_parts['bodies']:
            pc.set_facecolor(COLORS[0])
            pc.set_edgecolor(COLORS[1])
            pc.set_alpha(0.5)
        
        violin_parts['cmeans'].set_edgecolor('red')
        violin_parts['cmedians'].set_edgecolor('black')
        violin_parts['cbars'].set_edgecolor(None)
        violin_parts['cbars'].set_alpha(0)
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        # stat,p_val = kruskal(data1, data2)
        plt.title(title_string,fontsize=fsize)#+"\n stat = %0.2e, p-value = %0.2e")%(stat,p_val),fontsize=fsize)
        fig.tight_layout()
        if save_it == 1:
            plt.savefig(file_name+".png",dpi=150)
            plt.savefig(file_name+".svg",dpi=150)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
        
    def ClosePlot(self):
        plt.close()

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



if __name__ == '__main__':
    num_cells = 17
    map2_5_data,binned_data_dict,data_length_wise,length_data,length_arr = LoadData("Gria1/cells/",[3],"MAP2",5,num_cells)
    breakpoint()
    soma_norm_MAp2 = {}
    mean_soma_norm_MAP2 = {}
    for l in Lengths[1:]:
        soma_norm_MAp2[l] = GetSumNormDistribution(data_length_wise[l])
        mean_soma_norm_MAP2[l] = 
    # breakpoint()