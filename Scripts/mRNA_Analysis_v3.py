#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:38:33 2022

@author: surbhitwagle
"""
import argparse
import json
import matplotlib.pyplot as plt
import MAP2_Analysis as mp2a
import os
import pandas as pd
import pickle
from pathlib import Path
from SNSPlottingWidget import SNSPlottingWidget
from Utility import *
# from scipy import stats

# The type of mRNA to analyse, possbile choices: "Gria1", "Gria2", "CNIH2", "All"

# number of cells for each mRNA 
num_cells = {"CNIH2":28,"Gria1":19,"Gria2":15}


# cells to exclude from mRNA analysis

Excluded_cells= {"CNIH2":[],"Gria1":[],"Gria2":[3,9,10]}

# channel_names


#  two scales_mRNA 
# scale1 = 2.408  #pixels/uM
# scale2 = 7.587 #pixels/uM
scales_mRNA = {"CNIH2":[2.408,7.587],"Gria1":[4.817,4.817],"Gria2":[4.817,7.587]}
# channels to include
channels = [1,3]

dend_data = {}
soma_data = {}
dend_data_meta = {}
soma_data_meta = {}

mRNA_COLOR_code_old = {"Gria1":'#606C38',"Gria2":'#283618',"CNIH2":'#DDA15E',"Camk2a":'#BC6C25'}

# based on the finalized color panel at: https://coolors.co/f94144-f3722c-f8961e-f9844a-f9c74f-90be6d-43aa8b-4d908e-577590-277da1
# colors number 1,3,5,7 selected from left from this panel
mRNA_COLOR_code_old1 = {"Gria1":'#F94144',"Gria2":'#F8961E',"CNIH2":'#F9C74F',"Camk2a":'#43AA8B'}

# based on the finalized color panel at: https://coolors.co/f72585-b5179e-7209b7-560bad-480ca8-3a0ca3-3f37c9-4361ee-4895ef-4cc9f0
# colors number 1,3,5,7 selected from left from this panel

mRNA_COLOR_code_final = {"Gria1":'#560BAD',"Gria2":'#07bebe',"CNIH2":'#4361EE',"Camk2a":'#7f3b8f'}

mRNA_COLOR_code = mRNA_COLOR_code_final


#  dendritic widths analysis
widths = [5.0,10,15.0]
Analysis_dend_width = 0

# Function used to create the mRNA data dictionary by reading individual files outputed from Puncta detection using SpyDen V1.
# not required anymore, directly read the files from pickle files

# def ReadFiles(mRNA_to_analyse,widths_to_analyse,soma_bins,dend_bins,bin_size,channel_3_mRNA,lengths):
#      soma_data,dend_data,dend_data_meta,soma_data_meta = {},{},{},{}
#      soma_total_count ,soma_cell_count ,soma_unit_count, soma_total_stat, soma_cell_stat, soma_unit_stat = {},{},{},{},{},{}
#      dend_total_count ,dend_cell_count ,dend_unit_count, dend_total_stat, dend_cell_stat, dend_unit_stat = {},{},{},{},{},{}
#
#      count1 = 0
#      # MAP2_mean = {}
#      # MAP2_std = {}
#      # MAP2_sem = {}
#      MAP2_Sum_norm = {}
#      for width in widths_to_analyse:
#          #  intializing dictionaries to store processed data length wise
#         soma_data[width],dend_data[width],dend_data_meta[width],soma_data_meta[width] = {},{},{},{}
#         soma_total_count[width] ,soma_cell_count[width] ,soma_unit_count[width], soma_total_stat[width], \
#             soma_cell_stat[width], soma_unit_stat[width] = {},{},{},{},{},{}
#         dend_total_count[width] ,dend_cell_count[width] ,dend_unit_count[width], dend_total_stat[width]\
#             , dend_cell_stat[width], dend_unit_stat[width] = {},{},{},{},{},{}
#
#
#
#         channel_3_count = 1
#         dend_data[width][channel_3_mRNA] = {}
#         soma_data[width][channel_3_mRNA] = {}
#
#         dend_data_meta[width][channel_3_mRNA] = {}
#         soma_data_meta[width][channel_3_mRNA] = {}
#         # MAP2_mean[width] = []
#         # MAP2_std[width] = []
#         # MAP2_sem[width] = []
#         MAP2_Sum_norm[width] = {}
#         for mRNA in mRNA_to_analyse:
#             print("call MAP2 analysis")
#             mp2a.plotAndSaveMap2(num_cells, Excluded_cells, dend_bins, bin_size, int(float(width)), mRNA,lengths)
#              #  intializing dictionaries to store processed data length wise
#             soma_data[width][mRNA],dend_data[width][mRNA],dend_data_meta[width][mRNA],soma_data_meta[width][mRNA] = {},{},{},{}
#             soma_total_count[width][mRNA] ,soma_cell_count[width][mRNA] ,soma_unit_count[width][mRNA], soma_total_stat[width][mRNA], \
#                 soma_cell_stat[width][mRNA], soma_unit_stat[width][mRNA] = {},{},{},{},{},{}
#             dend_total_count[width][mRNA] ,dend_cell_count[width][mRNA] ,dend_unit_count[width][mRNA], dend_total_stat[width][mRNA]\
#                 , dend_cell_stat[width][mRNA], dend_unit_stat[width][mRNA] = {},{},{},{},{},{}
#
#             cell_num = num_cells[mRNA]
#             cells_to_exclude = Excluded_cells[mRNA]
#             cells = range(1,cell_num)
#             total_d_count = 0;
#             #root folder, will be changed later if mRNA type is set to all
#
#             folder = os.path.abspath(os.getcwd())+"/../{}/cells/".format(mRNA)
#             Channel_names = {0:"Dapi",1:mRNA,2:"MAP2",3:"CAMKII"}
#             # for i in channels:
#
#             soma_data[width][mRNA] = {}
#             dend_data[width][mRNA] = {}
#
#             # MAP2_mean[width] = {}
#             # MAP2_std[width] = {}
#             # MAP2_sem[width] = {}
#             # MAP2_Sum_norm[width] = {}
#
#             for lx in lengths:
#                 Sum_norm_map2_folder = os.path.abspath(os.getcwd())+"/../{0}/MAP2-Figures/{1}/Sum_norm_MAP2_{1}_{2}.npy".format(mRNA,int(float(width)),lx)#folder+"../MAP2-Figures/{0}/Sum_norm_MAP2_{0}_{1}.npy".format(int(float(width)),lx)
#
#                 if not lx in  MAP2_Sum_norm[width].keys():
#                     MAP2_Sum_norm[width][lx] = np.load(Sum_norm_map2_folder)
#                 else:
#                     # breakpoint()
#                     try:
#                         MAP2_Sum_norm[width][lx] = np.concatenate((MAP2_Sum_norm[width][lx],np.load(Sum_norm_map2_folder)),axis=0)
#                     except:
#                         print("error in concatinating ",mRNA)
#                 # std_map2_folder = folder+"../MAP2-Figures/{0}/std_MAP2_{0}_{1}.npy".format(int(float(width)),lx)
#                 # MAP2_mean[width][mRNA][lx] = np.load(mean_map2_folder)
#                 # MAP2_std[width][mRNA][lx] = np.load(std_map2_folder)
#                 # MAP2_sem[width][mRNA][lx] = MAP2_std[width][mRNA][lx]/cell_num
#             # breakpoint()
#             for cell in cells:
#                 if cell not in cells_to_exclude:
#                     # print(cell)
#                      file =  folder+"cell_{0}/{1}/{2}/dend_puncta_channel{1}_{2}.json".format(cell,1,width)
#                      file2 = folder+"cell_{0}/{1}/{2}/dend_puncta_channel{1}_{2}.json".format(cell,2,width)
#                      file3 = folder+"cell_{0}/{1}/{2}/dend_puncta_channel{1}_{2}.json".format(cell,3,width)
#                      # print(file)
#                      # breakpoint()
#                      if Path(file).is_file():
#                         # print("file exists")
#                         dend_data_meta[width][mRNA][cell],soma_data_meta[width][mRNA][cell] = LoadDend(folder+"cell_{0}/".format(cell),scales_mRNA[mRNA][0])
#                         total_d_count += len(dend_data_meta[width][mRNA][cell].keys())
#                         f_d = open(file)
#                         dend_data[width][mRNA][cell] = json.load(f_d)
#                         f_d.close()
#                         f_s = open(folder+"cell_{0}/{1}/{2}/soma_puncta_channel{1}_{2}.json".format(cell,1,width))
#                         # print(f_d,f_s)
#                         soma_data[width][mRNA][cell] = json.load(f_s)
#                         f_s.close()
#
#                         dend_data_meta[width][channel_3_mRNA][channel_3_count],soma_data_meta[width][channel_3_mRNA][channel_3_count] = LoadDend(folder+"cell_{0}/".format(cell),scales_mRNA[mRNA][0])
#                         f_d = open(file3)
#                         dend_data[width][channel_3_mRNA][channel_3_count] = json.load(f_d)
#                         f_d.close()
#                         f_s = open(folder+"cell_{0}/{1}/{2}/soma_puncta_channel{1}_{2}.json".format(cell,3,width))
#                         # print(f_d,f_s)
#                         soma_data[width][channel_3_mRNA][channel_3_count] = json.load(f_s)
#                         f_s.close()
#                         channel_3_count += 1
#                      if Path(file2).is_file():
#                         dend_data_meta[width][mRNA][cell],soma_data_meta[width][mRNA][cell] = LoadDend(folder+"cell_"+str(cell)+"/",scales_mRNA[mRNA][1])
#                         total_d_count += len(dend_data_meta[width][mRNA][cell].keys())
#                         f_d = open(file2)
#                         dend_data[width][mRNA][cell] = json.load(f_d)
#                         f_d.close()
#                         f_s = open(folder+"cell_{0}/{1}/{2}/soma_puncta_channel{1}_{2}.json".format(cell,2,width))
#                         # print(f_d,f_s)
#                         soma_data[width][mRNA][cell] = json.load(f_s)
#                         f_s.close()
#             #  storing dictionaries with  processed data cell, unit wise in binned format
#
#             soma_total_count[width][mRNA] ,soma_cell_count[width][mRNA] ,soma_unit_count[width][mRNA], soma_total_stat[width][mRNA], \
#                 soma_cell_stat[width][mRNA], soma_unit_stat[width][mRNA] = GetPunctaDicts(soma_data[width][mRNA],bins=soma_bins)
#             dend_total_count[width][mRNA] ,dend_cell_count[width][mRNA] ,dend_unit_count[width][mRNA], dend_total_stat[width][mRNA], \
#                 dend_cell_stat[width][mRNA], dend_unit_stat[width][mRNA] = GetPunctaDicts(dend_data[width][mRNA],bins=dend_bins)
#         soma_total_count[width][channel_3_mRNA] ,soma_cell_count[width][channel_3_mRNA] ,soma_unit_count[width][channel_3_mRNA], soma_total_stat[width][channel_3_mRNA], \
#                 soma_cell_stat[width][channel_3_mRNA], soma_unit_stat[width][channel_3_mRNA] = GetPunctaDicts(soma_data[width][channel_3_mRNA],bins=soma_bins)
#         dend_total_count[width][channel_3_mRNA] ,dend_cell_count[width][channel_3_mRNA] ,dend_unit_count[width][channel_3_mRNA], dend_total_stat[width][channel_3_mRNA], \
#                 dend_cell_stat[width][channel_3_mRNA], dend_unit_stat[width][channel_3_mRNA] = GetPunctaDicts(dend_data[width][channel_3_mRNA],bins=dend_bins)
#      c_ls = GetLengthCounts(dend_data_meta)
#
#      MAP2_mean = {}
#      MAP2_std = {}
#      MAP2_sem = {}
#      for width in widths_to_analyse:
#          MAP2_mean[width] = {}
#          MAP2_std[width] = {}
#          MAP2_sem[width] = {}
#          for lx in lengths:
#               MAP2_mean[width][lx] = MAP2_Sum_norm[width][lx].mean(axis=0)
#               MAP2_std[width][lx] = MAP2_Sum_norm[width][lx].std(axis=0)
#               MAP2_sem[width][lx] = MAP2_std[width][lx]/MAP2_Sum_norm[width][lx].shape[0]
#      # breakpoint()
#      return soma_data,dend_data,dend_data_meta,\
#          soma_total_count ,soma_cell_count ,soma_unit_count, soma_total_stat, soma_cell_stat, soma_unit_stat,\
#              dend_total_count ,dend_cell_count ,dend_unit_count, dend_total_stat, dend_cell_stat, dend_unit_stat,c_ls,MAP2_mean,MAP2_std,MAP2_sem
     
def GetCellWiseRatio(soma,dend):
    
   
    """
    function to calculate ratio of RNA between dendrites and soma
     input:  soma puncta dict type
             dendrite puncta dict
     
     returns: dict with channels and cells as key containing dend/soma ratio for each channel  
    """
    
    cell_wise_ratio = {}
    
    # for channel in soma.keys():
    #     cell_wise_ratio[channel] = []
    #     count = 0
    #     for cell in soma[channel].keys():
    #         sum_dend  = GetMatSum(dend[channel][cell])
    #         # breakpoint()
    #         cell_wise_ratio[channel].append(np.divide(sum_dend,soma[channel][cell][0]))
    #     cell_wise_ratio[channel] = np.asarray(cell_wise_ratio[channel])
    # # breakpoint()
    for width in soma.keys():
        cell_wise_ratio[width] = {}
        for mrna in soma[width].keys():
            cell_wise_ratio[width][mrna] = []
            count = 0
            for cell in soma[width][mrna].keys():
                 sum_dend  = GetMatSum(dend[width][mrna][cell])
                 # breakpoint()
                 cell_wise_ratio[width][mrna].append(np.divide(sum_dend,soma[width][mrna][cell][0]))
            cell_wise_ratio[width][mrna] = np.asarray(cell_wise_ratio[width][mrna])
    return cell_wise_ratio

def GetSomaticDendriticFraction(soma,dend):
    cell_wise_fraction = {}
    cells = []
    total = {}
    for width in soma.keys():
        cell_wise_fraction[width] = {}
        total[width] = {}
        for mrna in soma[width].keys():
            cell_wise_fraction[width][mrna] = []
            total[width][mrna] = {}
            count = 0
            for cell in soma[width][mrna].keys():
                sum_soma = soma[width][mrna][cell][0]
                sum_dend  = GetMatSum(dend[width][mrna][cell])
                total[width][mrna][cell]  = np.add(sum_soma, sum_dend)
                soma_fraction = np.divide(sum_soma,total[width][mrna][cell])
                dend_fraction = np.divide(sum_dend,total[width][mrna][cell])
                # breakpoint()
                # cells.append(cell)
                cell_wise_fraction[width][mrna].append([soma_fraction,dend_fraction])
            cell_wise_fraction[width][mrna] = np.asarray(cell_wise_fraction[width][mrna])
            # print(cells)
    return cell_wise_fraction, total#, cells


def GetDendriticSpatialDistribution(soma_d,dend_d,dend_meta_data):
    dendritic_dist = {}
    for width in dend_d.keys():
        dendritic_dist[width] = {}
        for mrna in dend_d[width].keys():
            dendritic_dist[width][mrna] = {}
            for l in Lengths:
                dendritic_dist[width][mrna][l] = []
                    # dendritic_dist_cell_wise[l][channel] = {}
    for width in dend_d.keys():
        for mrna in dend_d[width].keys():
            for cell in dend_d[width][mrna].keys():
                for unit in dend_d[width][mrna][cell].keys():
                    for l in Lengths:
                        try:
                            if dend_meta_data[width][mrna][cell][unit] >= l:
                                    dendritic_dist[width][mrna][l].append(dend_d[width][mrna][cell][unit])    
                        except:
                            breakpoint()
    
    for width in dend_d.keys():
        for mrna in dend_d[width].keys():
            for cell in dend_d[width][mrna].keys():
                for l in Lengths:
                    dendritic_dist[width][mrna][l]= np.asarray(dendritic_dist[width][mrna][l])
    return dendritic_dist



def GetPunctaDicts(od,bins=None):
    total_count = 0
    cell_count = {}
    unit_count = {}
    total_stat = []
    cell_stat = {}
    unit_stat = {}
    
    
    # for channel in od.keys():
    #     total_count[channel] = 0
    #     cell_count[channel] = {}
    #     unit_count[channel] = {}
    #     total_stat[channel] = []
    #     cell_stat[channel] = {}
    #     unit_stat[channel] = {}
        
    for cell in od.keys():
        cell_count[cell] = 0
        unit_count[cell] = {}
        
        cell_stat[cell] = []
        unit_stat[cell] = {}
        
        unique_punctas = []
        for unit in od[cell].keys():
            all_punctas = od[cell][unit]
            total_count += 1
            cell_count[cell] += 1
            unit_count[cell][unit] = len(all_punctas)
            all_stats = GetAllPunctastat(all_punctas)
            if all_stats.shape[0] != 0:
                binned_sum = BinnedSum(all_stats,bins,-1,"{0}_{1}_{2}".format("no channel",cell,unit)) # binning based on lenght  
            else:
                # breakpoint()
                binned_sum = np.zeros((bins.shape[0],GetRNAPunctaStat(np.zeros((10,1))).shape[0]))
            unit_stat[cell][unit] = binned_sum
            
            if cell_stat[cell] == []:
                cell_stat[cell] = binned_sum
            else:
                # breakpoint()
                cell_stat[cell] = np.add(cell_stat[cell],binned_sum)
            if total_stat == []:
                total_stat = binned_sum
            else:
                total_stat= np.add(total_stat,binned_sum)
            # breakpoint()
            # all_unique_stats = GetAllPunctastat(unique_punctas)
            # cell_unique_stats[channel][cell] = BinnedSum(all_stats,bins,-1)
                # print(max_length)
    return total_count ,cell_count ,unit_count, total_stat,cell_stat,unit_stat
            




    



def GetAllPunctastat(all_punctas):
    all_stats = []
    for puncta in all_punctas:
        all_stats.append(GetRNAPunctaStat(puncta))
    return np.asarray(all_stats)





def GetLengthCounts(meta_dict):
    l_count=np.zeros(len(Lengths))
    for width in meta_dict.keys():
        for mRNA in meta_dict[width].keys():
            for cell in meta_dict[width][mRNA].keys():
                for unit in meta_dict[width][mRNA][cell]:
                    count = (meta_dict[width][mRNA][cell][unit] >= Lengths)
                    # if (meta_dict[width][mRNA][cell][unit] >= 100):
                    #     print(cell,unit,meta_dict[width][mRNA][cell][unit])
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

# def ChiSq(yd,y_fit,sigmas):
#     nzs = np.nonzero(sigmas)
#     print(nzs)
#     r_yd = np.take(yd,nzs)
#     r_yf = np.take(y_fit,nzs)
#     r_sgs = np.take(sigmas,nzs)
#     residuals = r_yd - r_yf
#     chi_squ = np.sum((residuals/r_sgs)**2)
#     return chi_squ


def DumpDict(datadict,fname):
    with open(fname,'wb') as outfile:
        pickle.dump(datadict,outfile, protocol=pickle.HIGHEST_PROTOCOL)
    print("{} saved!".format(fname))

def ReadDataDict(fname):
    with open(fname,'rb') as infile:
        d = pickle.load(infile)
    print("{} loaded!".format(fname))
    return d

def ReadData(mrna_name):
    # dend_cell_sum = {}
    curr_wd = os.getcwd()
    print(curr_wd)
    return [
        ReadDataDict(os.path.join(curr_wd,"Scripts/mRNA_data/{}/soma_data.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_data.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd,"Scripts/mRNA_data/{}/dend_data_meta.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_total_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_cell_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_unit_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_total_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_cell_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/soma_unit_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_total_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_cell_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_unit_count.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_total_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_cell_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/dend_unit_stat.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/c_ls.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/MAP2_mean.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/MAP2_std.pickle".format(mrna_name))),
        ReadDataDict(os.path.join(curr_wd, "Scripts/mRNA_data/{}/MAP2_sem.pickle".format(mrna_name)))
    ]

if __name__ == '__main__':
    
    """
        parsing argumenst on mRNA and widths to be analysed. Each length is analysed seperately. mRNAs can be analysed in combinations
    """
    parser = argparse.ArgumentParser(description='mRNA analysis.py file ')
    parser.add_argument('-m', "--mRNA", nargs="+", default = ["all"],
                        help='mRNA names to analyse works for Gria1, Gria2, CNIH2 or any combination, if all is given all four will be plotted togather. Comparison \
                            is always made against CMAKII mRNA')
    parser.add_argument('-w', "--width", nargs="+", default = ['5.0'],
                        help='dendritic width, allowed values: 5.0, 10, 15.0 or any combination')
    
    # reading the argumenst
    args = parser.parse_args()
    analyse_all = False
   
    if "ALL" in (m.upper() for m in args.mRNA):
         analyse_all = True
         mRNA_to_analyse = ["Gria2","Gria1","CNIH2"]
    else:
        mRNA_to_analyse = args.mRNA
    widths_to_analyse = args.width
    
    
    # setting up analysis paramters
    bin_size = 7.5
    soma_bins = np.arange(0,1,bin_size)
    dend_bins = np.arange(0,Lengths.max(),bin_size)
    lengths_to_analyse = Lengths[1:-3]
    channel_3_mRNA = "Camk2a"
    
    # reading all the files necessary for the analysis
    soma_data,dend_data,dend_data_meta,\
         soma_total_count ,soma_cell_count ,soma_unit_count, soma_total_stat, soma_cell_stat, soma_unit_stat,\
             dend_total_count ,dend_cell_count ,dend_unit_count, dend_total_stat, dend_cell_stat, dend_unit_stat,c_ls,\
                 MAP2_mean,MAP2_std,MAP2_sem  = ReadData(mRNA_to_analyse[0])

    # DumpDict(soma_data, "./mRNA_data/{}/soma_data.pickle".format(mRNA_to_analyse[0]))
    # breakpoint()

    #  pltotting class and setting
    pw = SNSPlottingWidget()
    
    in_set = 1  #plot in_set normalized plots ?
    save_it = 1 #save the plots or not (=1 for yes)
    ax_label = 1 #plot axis labels or not (=1 for yes)
    w_or_wo_ax_label = ["", "_with_ax_label"] #for chaning the filename suffix
    #  calculating fractions for total in soma and dendrites
    fractions, total_cell_counts = GetSomaticDendriticFraction(soma_cell_stat,dend_cell_stat)
    data_to_show = {}
    # breakpoint()
    #  statistics number to use as representation of mRNA copy-number
    
    stats_list = ["area","area X mu","area X med","mean","med","Puncta Count"]
    stat_no = 5; #5 is puncta count
    
    # plotting varaible initializing
    labels = {"Gria1":"Gria1","Gria2":"Gria2","CNIH2":"Cnih2","Camk2a":"Camk2a"}
    lab1 = []
    colors = []
    compartment = ["Soma","Dendrite"]
    x_lab,y_lab = ["mRNA species","mRNA fraction of total"]
    title = ""#"mRNA localization dendrite vs somata";
    
    # breakpoint()
    # setting up the data,labels, colors to plot in the correct format
    factions_data = pd.DataFrame
    for width in fractions.keys():
        data_to_show[width] = []
        for mrna in fractions[width].keys():
            # breakpoint()
            print(mrna,fractions[width][mrna][:,:,stat_no].shape)
            data_to_show[width].append(fractions[width][mrna][:,0,stat_no])
            data_to_show[width].append(fractions[width][mrna][:,1,stat_no])
            lab1.append(labels[mrna])
            colors.append(mRNA_COLOR_code[mrna])
        data_to_show[width] = np.asarray(data_to_show[width])
    # breakpoint()
    # loop to call the plotting function
    for width in widths_to_analyse:
        op_folder = os.path.abspath(os.getcwd())+"/Figures/{}/{}/".format(width,'_'.join(mRNA_to_analyse))
        pw.CreateFolderRecursive(op_folder)
        # breakpoint()
        pw.PlotCellFraction(data_to_show[width],
                            lab1,
                            compartment,
                            x_lab,
                            y_lab,
                            colors,
                            title,
                            op_folder+"soma_vs_dend_fractions_{0}_{1}_{2}".format(stats_list[stat_no],width,w_or_wo_ax_label[ax_label]),
                            [],
                            save_it = save_it,
                            set_axis_label=ax_label)

        

    # calculating the ratio of dendrite to soma mrna count 
    dend_soma_ratio = GetCellWiseRatio(soma_cell_stat,dend_cell_stat)
    
    # breakpoint()
    # plotting varaible initializing
    x_lab,y_lab = ["mRNA species",'mRNA count ratio (Dendrite/Soma)']
    title = ""#"Ratio between dendritic and somatic mRNA copy-number"
    file_prefix = "dend_vs_soma_ratio_"
    vp_data_to_show = {}
    vp_labs = []
    vp_colors = []
    
    # setting up the data,labels, colors to plot in the correct format
    for width in dend_soma_ratio.keys():
        vp_data_to_show[width] = []
        for mrna in dend_soma_ratio[width].keys():
            print(mrna,dend_soma_ratio[width][mrna][:,stat_no].shape)
            vp_data_to_show[width].append(dend_soma_ratio[width][mrna][:,stat_no])
            vp_labs.append(labels[mrna])
            vp_colors.append(mRNA_COLOR_code[mrna])
    
        vp_data_to_show[width] = np.asarray(vp_data_to_show[width])
    
    # loop to call the plotting function
    for width in widths_to_analyse:
        op_folder = os.path.abspath(os.getcwd())+"/Figures/{}/{}/".format(width,'_'.join(mRNA_to_analyse))
        pw.CreateFolderRecursive(op_folder)    
        pw.ViolinPlotStats(vp_data_to_show[width],vp_labs,x_lab,y_lab,vp_colors ,title,op_folder+file_prefix+"_{0}_{1}_{2}".format(stats_list[stat_no],width,w_or_wo_ax_label[ax_label])\
                           ,[dend_soma_ratio.keys()],save_it = save_it,set_axis_label=ax_label)
    
    
    # getting dendritic distribution of mRNA counts
    dend_sum_norm_distribution_dict = GetDendriticSpatialDistribution(None,dend_unit_stat,dend_data_meta) 
    # breakpoint()
    total_puncta_dist = {}
    
    for width in widths_to_analyse:
        total_puncta_dist[width] = {}
        for mrna in dend_sum_norm_distribution_dict[width].keys():
            total_puncta_dist[width][mrna] = {}
            for ldx,l in enumerate(lengths_to_analyse):
               l1 = l #str(l)
               if dend_sum_norm_distribution_dict[width][mrna][l1].shape[0] >0:
                   total_puncta_dist[width][mrna][l1] = dend_sum_norm_distribution_dict[width][mrna][l1][:,:,stat_no].sum(axis=0)
               else:
                   total_puncta_dist[width][mrna][l1] = np.zeros((100,1))
    
    # means = {}
    # stds = {}
    # sems = {}
    # breakpoint()
    MAP2_norm_means = {}
    MAP2_nrom_stds = {}
    MAP2_norm_sems = {}
    # breakpoint()
    count_fittings = [0,0,0,1,0]
    norm_fittings = [0,0,0,0,0 ]

    for width in  widths_to_analyse:
        for mrna in mRNA_to_analyse:
            labs = [labels[mrna],labels[channel_3_mRNA]]
            x_lab,y_lab,y_lab_norm = [r"Dendritic distance ($\mu m$)",'mRNA puncta density',"Normalized \n mRNA density"]
            title = ""#"Spatial distribution of mRNA copy-number"
            file_prefix = "Spatial_mRNA_distribution"
            plot_colors = [mRNA_COLOR_code[mrna],mRNA_COLOR_code[channel_3_mRNA]]
            for ldx,l1 in enumerate(lengths_to_analyse):
                dend_dist = {}
                xs = np.arange(0,l1,bin_size)
                means = np.zeros((2,xs.shape[0]))
                stds = np.zeros((2,xs.shape[0]))
                if in_set == 1:
                    # breakpoint()
                    MAP2_norm_means = np.zeros(means.shape)
                    MAP2_norm_stds = np.zeros(means.shape)
                    norm_density_mean = np.zeros(means.shape)
                    norm_density_std = np.zeros(means.shape)

                    means[0] = dend_sum_norm_distribution_dict[width][mrna][l1][:,:,stat_no].mean(axis=0)[1:xs.shape[0]+1]
                    stds[0] = dend_sum_norm_distribution_dict[width][mrna][l1][:,:,stat_no].std(axis=0)[1:xs.shape[0]+1]
                    
                    # breakpoint()
                    MAP2_norm_data = dend_sum_norm_distribution_dict[width][mrna][l1][:,:,stat_no][:,1:xs.shape[0]+1]/MAP2_mean[width][l1][1:xs.shape[0]+1]
                    MAP2_norm_data_camKII = dend_sum_norm_distribution_dict[width][channel_3_mRNA][l1][:,:,stat_no][:,1:xs.shape[0]+1]/MAP2_mean[width][l1][1:xs.shape[0]+1]
                    
                    MAP2_norm_data = (MAP2_norm_data.T/MAP2_norm_data[:,0]).T
                    MAP2_norm_data_camKII = (MAP2_norm_data_camKII.T/MAP2_norm_data_camKII[:,0]).T
                    
                    means[1] = dend_sum_norm_distribution_dict[width][channel_3_mRNA][l1][:,:,stat_no].mean(axis=0)[1:xs.shape[0]+1]
                    stds[1] = dend_sum_norm_distribution_dict[width][channel_3_mRNA][l1][:,:,stat_no].std(axis=0)[1:xs.shape[0]+1]
                    # breakpoint()
                    MAP2_norm_means = means/MAP2_mean[width][l1][1:xs.shape[0]+1]
                    MAP2_norm_stds = ((stds/means ) * MAP2_norm_means) #+ (MAP2_std[width][l1]/MAP2_mean[width][l1])[1:xs.shape[0]+1]
                    MAP2_norm_stds[0] = MAP2_norm_stds[0]/np.sqrt(dend_sum_norm_distribution_dict[width][mrna][l1][:,:,stat_no].shape[0])
                    MAP2_norm_stds[1] = MAP2_norm_stds[1]/np.sqrt(dend_sum_norm_distribution_dict[width][channel_3_mRNA][l1][:,:,stat_no].shape[0])
                    # MAP2_norm_means[cdx] = MAP2_norm_means[cdx]/MAP2_norm_means[cdx,0]
                    
                    norm_density_mean = (MAP2_norm_means.T/MAP2_norm_means[:,0]).T
                    # breakpoint()
                    norm_density_std[0] = ((MAP2_norm_stds[0]/MAP2_norm_means[0]) * norm_density_mean[0])# + MAP2_norm_stds[0,0]/MAP2_norm_means[0,0]/np.sqrt(dend_sum_norm_distribution_dict[width][mrna][l][:,:,stat_no].shape[0])
                    norm_density_std[1] = ((MAP2_norm_stds[1]/MAP2_norm_means[1]) * norm_density_mean[1])# + MAP2_norm_stds[1,0]/MAP2_norm_means[1,0]/np.sqrt(dend_sum_norm_distribution_dict[width][channel_3_mRNA][l][:,:,stat_no].shape[0])
                    # norm_density_std[0] = 2*MAP2_norm_stds[0]/MAP2_norm_means[0]
                    # norm_density_std[1:] = (MAP2_norm_stds[1:]/MAP2_norm_means[1:]+MAP2_norm_stds[0]/MAP2_norm_means[0])*norm_density_mean[1:]
                    # np.nan_to_num(norm_density_std
                
                MAP2_norm_stds = np.nan_to_num(MAP2_norm_stds) 
                norm_density_std = np.nan_to_num(norm_density_std)
                # removing all sample that has infinite or nan entries after nomalization 
                MAP2_norm_data_camKII = MAP2_norm_data_camKII[~np.isnan(MAP2_norm_data_camKII).any(axis=1)]
                MAP2_norm_data = MAP2_norm_data[~np.isnan(MAP2_norm_data).any(axis=1)]
                MAP2_norm_data_camKII = MAP2_norm_data_camKII[~np.isinf(MAP2_norm_data_camKII).any(axis=1)]
                MAP2_norm_data = MAP2_norm_data[~np.isinf(MAP2_norm_data).any(axis=1)]
                # breakpoint()
                # pw.PlotFittedCurves(xs, MAP2_norm_data, MAP2_norm_data_camKII, labs, x_lab, y_lab, y_lab_norm,plot_colors,title+"_curve_fit", op_folder+file_prefix+"_norm_{0}_{1}_{2}_len_{3}".format(stats_list[stat_no],mrna,width,l1)\
                #                   ,bin_size,save_it = save_it,set_axis_label=ax_label,exp_method="NormE")
                # #
                # breakpoint()
                pw.PlotBinnedStats(np.asarray([xs,xs]), MAP2_norm_means, MAP2_norm_stds,norm_density_mean, labs, x_lab, y_lab, y_lab_norm,plot_colors,title, op_folder+file_prefix+"_norm_{0}_{1}_{2}_len_{3}_{4}".format(stats_list[stat_no],mrna,width,l1,w_or_wo_ax_label[ax_label])\
                                  ,bin_size,save_it = save_it,fit_exp=count_fittings[ldx],in_set=in_set,set_axis_label=ax_label,exp_method="1E")
                # breakpoint()
                pw.PlotBinnedStats(np.asarray([xs,xs]), norm_density_mean, norm_density_std, norm_density_mean, labs, x_lab, y_lab, y_lab_norm, plot_colors,title, op_folder+file_prefix+"_norm__with_fit_{0}_{1}_{2}_len_{3}_{4}".format(stats_list[stat_no],mrna,width,l1,w_or_wo_ax_label[ax_label])\
                                    , bin_size,save_it = save_it,fit_exp=norm_fittings[ldx],in_set=0,set_axis_label=ax_label)
                # breakpoint()
                # op_folder = os.path.abspath(os.getcwd())+"/Figures/{}/{}/".format(width,'_'.join([channel_3_mRNA]))
                # pw.CreateFolderRecursive(op_folder)
                # pw.PlotBinnedStats1p(np.asarray([xs]), np.asarray([norm_density_mean[1]]), np.asarray([norm_density_std[1]]), norm_density_mean[1], [r'CaMKII$\alpha$'], x_lab, y_lab, y_lab_norm, plot_colors,title, op_folder+file_prefix+"_norm__with_fit_{0}_{1}_{2}_len_{3}".format(stats_list[stat_no],channel_3_mRNA,width,l1)\
                #                     , bin_size,save_it = 1,fit_exp=1,in_set=0,set_axis_label=0)
                        
        