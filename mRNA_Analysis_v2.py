#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 13:55:38 2021

@author: surbhitwagle
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import json
import lmfit
from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
from operator import add
import os
import pandas as pd
from pathlib import Path
from functools import reduce
from scipy.stats import ks_2samp, kruskal
import seaborn as sb
from scipy.optimize import curve_fit
import scikit_posthocs as sp
import seaborn as sns
import MAP2_Analysis as mp2a

# from scipy import stats

# The type of mRNA to analyse, possbile choices: "Gria1", "Gria2", "CNIH2", "All"
mRNA = "Gria1"
# number of cells for each mRNA 
num_cells = {"CNIH2":28,"Gria1":19,"Gria2":15}

#root folder, will be changed later if mRNA type is set to all

folder = os.path.abspath(os.getcwd())+"/{}/cells/".format(mRNA)
# cells to exclude from mRNA analysis

Excluded_cells= {"CNIH2":[],"Gria1":[],"Gria2":[3,9,10]}

# channel_names
Channel_names = {0:"Dapi",1:mRNA,2:"MAP2",3:"CAMKII"}

#  two scales 
scale1 = 4.817  #pixels/uM
scale2 = 7.587 #pixels/uM

# channels to include
channels = [1,3]

dend_data = {}
soma_data = {}
dend_data_meta = {}
soma_data_meta = {}

mRNA_COLOR_code = {"Gria1":'#EF476F',"Gria2":'#FFD166',"CNIH2":'#06D6A0',"CAMKII":'#118AB2'}

# possbile lengths for dendritic distributions
Lengths = np.asarray([2,25,50,100,150,200,250])

#  dendritic widths analysis
widths = [5.0,10,15.0]
Analysis_dend_width = 0 
def ReadFiles(dir,cell_num,cells_to_exclude,width,in_set=1):
     cells = range(1,cell_num)
     total_d_count = 0;
     for i in channels:
         soma_data['channel_'+str(i)] = {}
         dend_data['channel_'+str(i)] = {}
         for cell in cells:
             if cell not in cells_to_exclude:
                 # print(cell)
                  file =  folder+"cell_{0}/{1}/{2}/dend_puncta_channel{1}_{2}.json".format(cell,i,width)
                  file2 = folder+"cell_{0}/{1}/{2}/dend_puncta_channel{1}_{2}.json".format(cell,2,width)
                  # print(file)
                  # breakpoint()
                  if Path(file).is_file():
                     # print("file exists")
                     dend_data_meta['cell_{}'.format(cell)],soma_data_meta['cell_{0}'.format(cell)] = LoadDend(folder+"cell_{0}/".format(cell),scale1)
                     total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                     f_d = open(file)
                     dend_data['channel_'+str(i)]['cell_'+str(cell)] = json.load(f_d)
                     f_d.close()
                     f_s = open(folder+"cell_{0}/{1}/{2}/soma_puncta_channel{1}_{2}.json".format(cell,i,width))
                     # print(f_d,f_s)
                     soma_data['channel_'+str(i)]['cell_'+str(cell)] = json.load(f_s)
                     f_s.close()
                  if Path(file2).is_file():
                     dend_data_meta['cell_'+str(cell)],soma_data_meta['cell_'+str(cell)] = LoadDend(folder+"cell_"+str(cell)+"/",scale2)
                     total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                     f_d = open(file2)
                     dend_data['channel_1']['cell_'+str(cell)] = json.load(f_d)
                     f_d.close()
                     f_s = open(folder+"cell_{0}/{1}/{2}/soma_puncta_channel{1}_{2}.json".format(cell,2,width))
                     # print(f_d,f_s)
                     soma_data['channel_1']['cell_'+str(cell)] = json.load(f_s)
                     f_s.close()
                 
                 
     bin_size = 5
     soma_bins = np.arange(0,1,bin_size)
     dend_bins = np.arange(0,Lengths.max(),bin_size)
    
     soma_total_count ,soma_cell_count ,soma_unit_count, soma_total_stat, soma_cell_stat, soma_unit_stat = GetPunctaDicts(soma_data,bins=soma_bins)
     dend_total_count ,dend_cell_count ,dend_unit_count, dend_total_stat, dend_cell_stat, dend_unit_stat = GetPunctaDicts(dend_data,bins=dend_bins)
     c_ls = GetLengthCounts(dend_data_meta)
     
     dend_cell_sum = {}
     # for c in dend_cell_stat.keys():
     #     dend_cell_sum[c] = {}
     #     for cell in dend_cell_stat[c].keys():
     #         dend_cell_sum[c][cell] = np.array([GetMatSum(dend_cell_stat[c][cell])])
     # breakpoint()
     
     # dend_distribution_dict = GetDendriticSpatialDistribution(soma_cell_unique_stats,dend_unit_stat,dend_data_meta) 
     # breakpoint()
     dend_sum_norm_distribution_dict = GetDendriticSpatialDistribution(None,dend_unit_stat,dend_data_meta) 
     # breakpoint()
     
     total_puncta_dist = {}
     
     for ldx,l in enumerate(Lengths):
        l1 = l #str(l)
        total_puncta_dist[l1] = {}
        for cdx,c in enumerate(dend_sum_norm_distribution_dict[l].keys()):
            if dend_sum_norm_distribution_dict[l1][c].shape[0] >0:
                total_puncta_dist[l1][c] = dend_sum_norm_distribution_dict[l1][c][:,:,5].sum(axis=0)
            else:
                total_puncta_dist[l1][c] = np.zeros((100,1))
     pw = PlottingWidgetmRNA()
     
     op_folder = folder+"../FENS/{}/".format(width)
     pw.CreateFolderRecursive(op_folder)

    
     save_it = 1
     fractions, total_cell_counts = GetSomaticDendriticFraction(soma_cell_stat,dend_cell_stat)
     
     lab1= ['{0}-Soma'.format(Channel_names[channels[0]]),'{0}-Dendrite'.format(Channel_names[channels[0]]),'{0}-Soma'.format(Channel_names[channels[1]]),'{0}-Dendrite'.format(Channel_names[channels[1]])]
     x_lab,y_lab = ["mRNA","Fraction of mRNA in Soma Vs. dendrites"]
     title = "Quatification of Dendritic and Somatic mRNA copy-number";
     stat_no = 5;
     stats_list = ["area","area X mu","area X med","mean","med","Puncta Count"]
     data_to_show = np.asarray([fractions['channel_1'][:,0,stat_no],fractions['channel_1'][:,1,stat_no],fractions['channel_3'][:,0,stat_no],fractions['channel_3'][:,1,stat_no]])
     # breakpoint()
     ax_label = 0
     pw.PlotCellFraction(data_to_show,lab1,x_lab,y_lab,title,op_folder+"soma_vs_dend_fractions_{0}_{1}_{2}_{3}".format(stats_list[stat_no],Channel_names[channels[0]],\
                                                                                                                       Channel_names[channels[1]],width)\
                         ,[Channel_names[channels[0]],Channel_names[channels[1]]],save_it = save_it,set_axis_label=ax_label)
    
     soma_dend_ratio = GetCellWiseRatio(soma_cell_stat,dend_cell_stat)
     lab1,lab2 = [Channel_names[channels[0]],Channel_names[channels[1]]]
     x_lab,y_lab = ["mRNAs",'mRNA count ratio (Dendrite/Soma)']
     title = "Ratio between dendritic and somatic mRNA copy-number"
     file_prefix = "dend_vs_soma_ratio_"
     pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,stat_no],soma_dend_ratio['channel_3'][:,stat_no],lab1,lab2,x_lab,y_lab ,title,op_folder+file_prefix+"_{0}_{1}_{2}_{3}".format(stats_list[stat_no],Channel_names[channels[0]],\
                                                                                                                        Channel_names[channels[1]],width)\
                          ,[Channel_names[channels[0]],Channel_names[channels[1]]],save_it = save_it,set_axis_label=ax_label)
     
     labs = [Channel_names[channels[0]],Channel_names[channels[1]]]
     x_lab,y_lab,y_lab_norm = ["Dendritic distance",'mRNA puncta count ',"Normalized \n mRNA count"]
     title = "Spatial distribution of mRNA copy-number"
     file_prefix = "Spatial_mRNA_distribution"
     for l in Lengths[1:-2]:
         dend_dist = {}
         xs = np.arange(0,l,bin_size)
         means = np.zeros((2,xs.shape[0]))
         stds = np.zeros((2,xs.shape[0]))
         if in_set == 1:
             mean_map2_folder = folder+"../MAP2-Figures/{0}/mean_MAP2_{0}_{1}.npy".format(int(width),l)
             std_map2_folder = folder+"../MAP2-Figures/{0}/std_MAP2_{0}_{1}.npy".format(int(width),l)
             MAP2_mean = np.load(mean_map2_folder)
             MAP2_std = np.load(std_map2_folder)
         MAP2_norm_means = np.zeros(means.shape)
         MAP2_norm_stds = np.zeros(means.shape)
         norm_density_mean = np.zeros(means.shape)
         norm_density_std = np.zeros(means.shape)
         for cdx,c in enumerate(dend_sum_norm_distribution_dict[l].keys()):
             means[cdx] = dend_sum_norm_distribution_dict[l][c][:,:,stat_no].mean(axis=0)[1:xs.shape[0]+1]
             stds[cdx] = dend_sum_norm_distribution_dict[l][c][:,:,stat_no].std(axis=0)[1:xs.shape[0]+1]
             MAP2_norm_means[cdx] = means[cdx]/MAP2_mean[1:xs.shape[0]+1]
             MAP2_norm_stds[cdx] = ((stds[cdx]/means[cdx] + (MAP2_std/MAP2_mean)[1:xs.shape[0]+1]) * MAP2_norm_means[cdx])/np.sqrt(dend_sum_norm_distribution_dict[l][c].shape[0])
             # MAP2_norm_means[cdx] = MAP2_norm_means[cdx]/MAP2_norm_means[cdx,0]
             
             norm_density_mean[cdx] = MAP2_norm_means[cdx]/MAP2_norm_means[cdx,0]
             norm_density_std[cdx,0] = 2*MAP2_norm_stds[cdx,0]/MAP2_norm_means[cdx,0]
             norm_density_std[cdx,1:] = (MAP2_norm_stds[cdx,1:]/MAP2_norm_means[cdx,1:]+MAP2_norm_stds[cdx,0]/MAP2_norm_means[cdx,0])*norm_density_mean[cdx,1:]
             # np.nan_to_num(norm_density_std)
         MAP2_norm_stds = np.nan_to_num(MAP2_norm_stds) 
         norm_density_std = np.nan_to_num(norm_density_std) 
         pw.PlotBinnedStats(np.asarray([xs,xs]), MAP2_norm_means, MAP2_norm_stds,norm_density_mean, labs, x_lab, y_lab, y_lab_norm,title, op_folder+file_prefix+"_norm_{0}_{1}_{2}_{3}_len_{4}".format(stats_list[stat_no],Channel_names[channels[0]],\
                                                                                                                         Channel_names[channels[1]],width,l)\
                           ,bin_size,save_it = save_it,fit_exp=0,in_set=in_set,set_axis_label=ax_label)
         # breakpoint()
         pw.PlotBinnedStats(np.asarray([xs,xs]), norm_density_mean, norm_density_std, norm_density_mean, labs, x_lab, y_lab, y_lab_norm, title, op_folder+file_prefix+"_norm__with_fit_{0}_{1}_{2}_{3}_len_{4}".format(stats_list[stat_no],Channel_names[channels[0]],\
                                                                                                                         Channel_names[channels[1]],width,l)\
                             , bin_size,save_it = save_it,fit_exp=1,in_set=0,set_axis_label=ax_label)
def GetCellWiseRatio(soma,dend):
    
   
    """
    function to calculate ratio of RNA between dendrites and soma
     input:  soma puncta dict type
             dendrite puncta dict
     
     returns: dict with channels and cells as key containing dend/soma ratio for each channel  
    """
    
    cell_wise_ratio = {}
    
    for channel in soma.keys():
        cell_wise_ratio[channel] = []
        count = 0
        for cell in soma[channel].keys():
            sum_dend  = GetMatSum(dend[channel][cell])
            # breakpoint()
            cell_wise_ratio[channel].append(np.divide(sum_dend,soma[channel][cell][0]))
        cell_wise_ratio[channel] = np.asarray(cell_wise_ratio[channel])
    # breakpoint()
    return cell_wise_ratio

def GetSomaticDendriticFraction(soma,dend):
    cell_wise_fraction = {}
    cells = []
    total = {}
    for channel in soma.keys():
        cell_wise_fraction[channel] = []
        total[channel] = {}
        count = 0
        for cell in soma[channel].keys():
            sum_soma = soma[channel][cell][0]
            sum_dend  = GetMatSum(dend[channel][cell])
            total[channel][cell]  = np.add(sum_soma, sum_dend)
            soma_fraction = np.divide(sum_soma,total[channel][cell])
            dend_fraction = np.divide(sum_dend,total[channel][cell])
            # breakpoint()
            # cells.append(cell)
            cell_wise_fraction[channel].append([soma_fraction,dend_fraction])
        cell_wise_fraction[channel] = np.asarray(cell_wise_fraction[channel])
        # print(cells)
    return cell_wise_fraction, total#, cells


def GetDendriticSpatialDistribution(soma_d,dend_d,dend_meta_data):
    dendritic_dist = {}
    for l in Lengths:
        dendritic_dist[l] = {}
        for channel in dend_d.keys():
            dendritic_dist[l][channel] = []
            # dendritic_dist_cell_wise[l][channel] = {}
    for channel in dend_d.keys():
        for cell in dend_d[channel].keys():
            for unit in dend_d[channel][cell].keys():
                for l in Lengths:
                    if dend_meta_data[cell][unit] >= l:
                            dendritic_dist[l][channel].append(dend_d[channel][cell][unit])                
    for l in Lengths:
        for channel in dend_d.keys():
            dendritic_dist[l][channel] = np.asarray(dendritic_dist[l][channel])
    return dendritic_dist
                        
def GetPunctaDicts(od,bins=None):
    total_count = {}
    cell_count = {}
    unit_count = {}
    total_stat = {}
    cell_stat = {}
    unit_stat = {}
    
    
    for channel in od.keys():
        total_count[channel] = 0
        cell_count[channel] = {}
        unit_count[channel] = {}
        total_stat[channel] = []
        cell_stat[channel] = {}
        unit_stat[channel] = {}
        
        for cell in od[channel].keys():
            cell_count[channel][cell] = 0
            unit_count[channel][cell] = {}
            
            cell_stat[channel][cell] = []
            unit_stat[channel][cell] = {}
            
            unique_punctas = []
            for unit in od[channel][cell].keys():
                all_punctas = od[channel][cell][unit]
                total_count[channel] += 1
                cell_count[channel][cell] += 1
                unit_count[channel][cell][unit] = len(all_punctas)
                all_stats = GetAllPunctastat(all_punctas)
                if all_stats.shape[0] != 0:
                    binned_sum = BinnedSum(all_stats,bins,-1,"{0}_{1}_{2}".format(channel,cell,unit)) # binning based on lenght  
                else:
                    # breakpoint()
                    binned_sum = np.zeros((bins.shape[0],GetRNAPunctaStat(np.zeros((10,1))).shape[0]))
                unit_stat[channel][cell][unit] = binned_sum
                
                if cell_stat[channel][cell] == []:
                    cell_stat[channel][cell] = binned_sum
                else:
                    # breakpoint()
                    cell_stat[channel][cell] = np.add(cell_stat[channel][cell],binned_sum)
                if total_stat[channel] == []:
                    total_stat[channel] = binned_sum
                else:
                    total_stat[channel] = np.add(total_stat[channel],binned_sum)
            # breakpoint()
            # all_unique_stats = GetAllPunctastat(unique_punctas)
            # cell_unique_stats[channel][cell] = BinnedSum(all_stats,bins,-1)
                # print(max_length)
    return total_count ,cell_count ,unit_count, total_stat,cell_stat,unit_stat
            
def BinnedSum(arr,bins,num_col=-1,name= None):
    # print("binned_sum for =",name)
    
    if len(arr.shape) == 2:
        rr,cc = arr.shape 
        binned_sum = np.zeros((len(bins),cc))
        arr = SortPunctas(arr,num_col)
        # max_lenght = arr[-1,-1]
        digitized = bins.searchsorted(arr[:,num_col])
        
        #     breakpoint()
        for c in range(0,cc):
            try:
                binned_sum[:,c] = np.bincount(digitized, weights=arr[:,c], minlength=len(bins))
            except:
                print("puncta is ",arr)
                breakpoint()
        binned_sum[:,num_col] = bins
        return binned_sum
    else:
        print("quite not the shape",arr.shape)
        return np.zeros((len(bins),6))
def GetUniqueRows(mat):
    return np.unique(mat,axis = 0)


    


def GetDictSum(d):
    dict_sum = np.zeros(6)
    for key in d.keys():
        dict_dum = np.sum(dict_sum, GetRNAPunctaStat(d[key]))

def GetMatSum(mat,ax=0):
    return np.sum(mat,axis=ax)

def SortPunctas(ps,column=0):
    return ps[ps[:,column].argsort()]

def GetAllPunctastat(all_punctas):
    all_stats = []
    for puncta in all_punctas:
        all_stats.append(GetRNAPunctaStat(puncta))
    return np.asarray(all_stats)
def GetRNAPunctaStat(puncta):
    #the puncta is of the form (x,y,r,max,min,mean.std,median, insoma,distance_from_soma)
    x,y,r,mx,mn,mu,delta,med,inSoma,dfs = puncta
    
    area = Area(r)
    stat1 =  area*mu
    stat2 = area*med
    stat3 = mu
    stat4 = med
    stat5 = area
    return np.array([stat1,stat2,stat3,stat4,stat5,1,dfs])

# def GetDendriteRNAPuncta(puncta):
    #the puncta is of the form (x,y,r,max,min,mean.std,median, insoma,distance_from_soma)
    # x,y,r,mx,mn,mu,delta,med,inSoma,dfs = puncta
    
    # area = Area(r)
    # stat1 =  area*mu
    # stat2 = area*med
    # stat3 = mu
    # stat4 = med
    # stat5 = area
    # return [stat1,stat2,stat3,stat4,stat5,dfs]

def Area(radius):
    return np.pi*(radius**2)

def GetLengthCounts(meta_dict):
    l_count=np.zeros(len(Lengths))
    for cell in meta_dict.keys():
        for unit in meta_dict[cell]:
            # breakpoint()
            count = (meta_dict[cell][unit] >= Lengths)
            # if (meta_dict[cell][unit] >= 75):
                # print(cell,unit,meta_dict[cell][unit])
            # l_thras = int(meta_dict[cell][unit]/50) - 1 
            # if l_thras >= 0:
            #     l_count[l_thras] += 1
            l_count = np.add(l_count,count)
    return l_count

def LoadDend(Dir,scale):

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
    soma = {}
    for x in os.listdir(Dir):
        if x.startswith('AllDendrites') and x.endswith(".npy"):
            # print(x)
            # breakpoint()
            all_dend = np.load(Dir+x,allow_pickle=True)
            for dxd,dend in enumerate(all_dend):
                DendArr["dendrite_{}".format(dxd)]= GetDendlength(dend,scale)
            # total_d_count += 1;
            # DendArr_names.append()
        elif("Soma" in x) and x.endswith(".npy"):
            s = np.load(Dir+x)
            soma[x.split(".npy")[0]] = GetPolygonArea(s[:,0],s[:,1],scale)
    return DendArr,soma


def GetDendlength(xys,scale):
    length = 0;
    for idx in range(0,xys.shape[0]-1):
        length += Distance(xys[idx],xys[idx+1])
    return length/scale

def Distance(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def GetPolygonArea(x,y,scale):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))/(scale**2)

class PlottingWidgetmRNA():
    def __init__(self,fsize=16,tsize=25,fam='Source Code Pro',pixelden = 100,lw=3.0,width=10,height=8):
        rc('font',
            family=fam,
            size=fsize)
    
        rc('figure', 
            dpi=pixelden,
            titlesize=tsize,
            titleweight='heavy',
            figsize = (width,height))
    
        rc('axes', 
            # linewidth=lw,
            titlesize=20,
            titleweight='regular',
            labelsize=fsize,
            )
    
        rc('legend',
            fontsize=fsize)
        rc('xtick',
           labelsize=fsize)
        rc('ytick',
           labelsize=fsize)
        rc('boxplot.meanprops',
           linewidth=lw,
           linestyle='--')
        rc('boxplot.boxprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.capprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.capprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.flierprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.whiskerprops',
           linewidth=lw,
           linestyle='-')
        rc('boxplot.medianprops',
           linewidth=lw,
           linestyle='-')
        rc('lines',
           linewidth=lw,
           markersize=4*lw)
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
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    
    def PlotFourStats(self,x_data1,y_data1,x_data2,y_data2,lab1,lab2,xlab,ylab,title_string,file_name,save_it = 1):
        fig, ax = plt.subplots()
        
        ax.scatter(x_data1,y_data1,label=lab1)
        ax.scatter(x_data2,y_data2,label=lab2)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        plt.title(title_string)
        plt.legend()
        # plt.show()
        folder = "."
        if save_it == 1:
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    
    def PlotBinnedStats(self,xs,means,stds,MAP2_norm,labs,xlab,ylab,ylab_norm,title_string,file_name,bin_size,save_it = 1,fit_exp =1,auto_close=1,in_set=1,lw=3.0,set_axis_label=1):
        fig, ax = plt.subplots()
        ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0' ,markersize=4*lw)
        if in_set ==1:
            left, bottom, width, height = [0.35, 0.7, 0.2, 0.2]
            ax2 = fig.add_axes([left, bottom, width, height])
            ax2.set_ylabel(ylab_norm)
        if xs.shape[0] == means.shape[0]:
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i],means[i],stds[i],label=labs[i],color=mRNA_COLOR_code[labs[i]],marker='d',linestyle='None' )
               
               if in_set==1:
                   ax2.plot(xs[i],MAP2_norm[i],color=mRNA_COLOR_code[labs[i]],marker='o',markersize=4,linestyle='dashed')
            if set_axis_label == 1:
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            # plt.title(title_string,fontsize=fsize)
            
            folder = "."
            
            if fit_exp == 1:
                for i in range(xs.shape[0]):
                    # breakpoint()
                    yi_fit, ri_squared,chi_squ = ExpFit(xs[i],means[i],stds[i,:],0,+1,labs[i])
                    ax.plot(xs[i],yi_fit,marker='None',c=mRNA_COLOR_code[labs[i]],label=labs[i]+"-fit")
            
            # plt.show()
            fig.tight_layout()
            ax.legend()
            if save_it == 1: 
                self.SaveFigures(file_name)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
            if auto_close == 1:
                plt.close()
        else:
            print("Not same number of xs and means",xs.shape,means.shape)
    def PlotCellFraction(self,fractions,lab,xlab,ylab,title_string,file_name,molecules,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        pos = np.arange(1,3,1)
      
        bp1 = ax.boxplot(np.transpose(fractions[0:2]),widths = 0.25,positions=pos,labels=lab[0:2],showmeans=True,meanline=True,showfliers=False,meanprops = dict(color=mRNA_COLOR_code[molecules[0]]))
        bp2 = ax.boxplot(np.transpose(fractions[2:]),widths = 0.25,positions=2+pos,labels=lab[2:],showmeans=True,meanline=True,showfliers=False,meanprops = dict(color=mRNA_COLOR_code[molecules[1]]))
        breakpoint()
        p_values = sp.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        x_points = np.asarray((pos,2+pos)).flatten()
        pairs = np.array([[1,2],[3,4],[1,4],[2,4],[1,3],[2,3]])
        for idx,pair in enumerate(pairs):
            txt = ''
            print(p_values[pair[0]][pair[1]],pair)
            if p_values[pair[0]][pair[1]] <= 0.05:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.01:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.001:
                txt += '*'
            # breakpoint()
            y_max = np.array([fractions[pair[0]-1].max(),fractions[pair[1]-1].max()]).max()
            self.AnnotateText(ax,x_points[pair[0]-1],x_points[pair[1]-1],y_max,0.01,txt,'k')
        self.setBoxColors(bp1,mRNA_COLOR_code[molecules[0]])
        self.setBoxColors(bp2,mRNA_COLOR_code[molecules[1]])
        if set_axis_label == 1:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        plt.ylim([0,1.2])
        # plt.title(title_string,fontsize=fsize)
        # plt.tick_params(
        #     axis='x',          # changes apply to the x-axis
        #     which='both',      # both major and minor ticks are affected
        #     bottom=False,      # ticks along the bottom edge are off
        #     top=False,         # ticks along the top edge are off
        #     labelbottom=False) # labels along the bottom edge are off
        means = []
        stds = []
        for i in range(4):
            means.append(fractions[i].mean())
            stds.append(fractions[i].std())
        # breakpoint()
        for j, line in enumerate(bp1['means']):
            x, y = line.get_xydata()[1]
            text = r'${:.2f}(\pm {:.2f})$'.format(means[j], stds[j])
            ax.annotate(text, xy=(x, y))
        for j, line in enumerate(bp2['means']):
            x, y = line.get_xydata()[1]
            text = r' ${:.2f}(\pm {:.2f})$'.format(means[2+j], stds[2+j])
            ax.annotate(text, xy=(x, y))
        fig.tight_layout()
        folder = "."
        if save_it == 1:
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    
    def AnnotateText(self,ax,x1,x2,y,h,txt,color,ha='center',va='bottom'):
        # print(x,y,txt)
        if not txt == '':
            fac = np.abs(x2-x1)*0.08
            plt.plot([x1,x1, x2,x2], [y+fac, y+h+fac, y+h+fac, y+fac], lw=1.5, c=color)
            plt.text((x1+ x2)*0.5,y+h+fac,txt, ha=ha, va=va, color=color)
        
    def setBoxColors(self,bp,c):
        setp(bp['boxes'], color=c)
        setp(bp['caps'], color=c )
        setp(bp['caps'], color=c )
        setp(bp['whiskers'], color=c )
        setp(bp['whiskers'], color=c )
        setp(bp['fliers'], color=c )
        setp(bp['fliers'], color=c )
        setp(bp['medians'], color=c )
        # setp(bp['mean'], color=c )
    
    def PlotBinnedStatsMulti(self,x,data_1,lab1,xlab,ylab,title_string,file_name,width=10,height=8,fsize=16,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        
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
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    def SaveFigures(self,filename,ext_list = [".png",".svg",".pdf"]):
        for ext in ext_list:
            plt.savefig(filename+ext,dpi=300)
            
    def ViolinPlotStats(self,data1,data2,lab1,lab2,xlab,ylab,title_string,file_name,molecules,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        
        # ax.set_axis_style([lab1,lab2])
        labels=[lab1+"\n N=%d cells"%(data1.shape[0]),lab2+"\n N=%d cells"%(data2.shape[0])]
        x = np.arange(len(labels))
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        xs = np.zeros(data1.shape)
        violin_parts = ax.violinplot([data1,data2],positions=x,showextrema=True,showmeans = True,showmedians = True,points=data1.shape[0])
        ax.scatter(xs,data1,color=mRNA_COLOR_code[molecules[0]],alpha=0.6)
        x2s = np.zeros(data2.shape)+(x[1])
        ax.scatter(x2s,data2,color=mRNA_COLOR_code[molecules[0]],alpha=0.6)
        # print(violin_parts['bodies'])
        for pc in violin_parts['bodies']:
            pc.set_facecolor(mRNA_COLOR_code[molecules[1]])
            pc.set_edgecolor(mRNA_COLOR_code[molecules[1]])
            pc.set_alpha(0.5)
        
        violin_parts['cmeans'].set_edgecolor('red')
        violin_parts['cmedians'].set_edgecolor('black')
        violin_parts['cbars'].set_edgecolor(None)
        violin_parts['cbars'].set_alpha(0)
        if set_axis_label == 1:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        # stat,p_val = kruskal(data1, data2)
        plt.title(title_string)#+"\n stat = %0.2e, p-value = %0.2e")%(stat,p_val),fontsize=fsize)
        fig.tight_layout()
        if save_it == 1:
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
        
    def ClosePlot(self):
        plt.close()

def exponential(x, b):
    return np.exp(b*x)

def ExpFit(xdata,ydata,sigmas,Fidx,Lidx,molecule):
    # breakpoint()
    param_bounds=([-np.inf],[0])
    popt, pcov = curve_fit(exponential, xdata[Fidx:], ydata[Fidx:],bounds =param_bounds)
    # breakpoint()
    print(molecule,popt,pcov,np.sqrt(np.diag(pcov)))
    y_fit = exponential(xdata, *popt)
    residuals = ydata[Fidx:]- y_fit[Fidx:]
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata[Fidx:]-np.mean(ydata[Fidx:]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = ChiSq( ydata[Fidx:], y_fit[Fidx:],sigmas[Fidx:])
    print(chi_squ)
    return y_fit,r_squared,chi_squ

def ChiSq(yd,y_fit,sigmas):
    nzs = np.nonzero(sigmas)
    print(nzs)
    r_yd = np.take(yd,nzs)
    r_yf = np.take(y_fit,nzs)
    r_sgs = np.take(sigmas,nzs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals/r_sgs)**2)
    return chi_squ
def ExpFit2(xdata,ydata,sigmas,Fidx,Lidx,molecule):
    pars = Parameters()
    pars.add('amplitude',1,vary=False)
    pars.add('decay',1,min=0)
    mod = lmfit.models.ExponentialModel()
    out = mod.fit(ydata[Fidx:], pars, x=xdata[Fidx:])
    y_fit = exponential(xdata[Fidx:],-1.0/out.params['decay'])
    residuals = ydata[Fidx:]- y_fit[Fidx:]
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata[Fidx:]-np.mean(ydata[Fidx:]))**2)
    r_squared = 1 - (ss_res / ss_tot)
    chi_squ = np.sum((residuals/sigmas)**2)
    print("here",chi_squ)
    # breakpoint()
    return y_fit,r_squared,out.chisqr
if __name__ == '__main__':
    ReadFiles(folder,num_cells[mRNA],Excluded_cells[mRNA],widths[Analysis_dend_width])