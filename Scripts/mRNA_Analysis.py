#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 13:55:38 2021

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
mRNA = "CNIH2"
folder = os.path.abspath(os.getcwd())+"/{}/".format(mRNA)
num_cells = {"CNIH2":28,"Gria1":19,"Gria2":15}
Excluded_cells= {"CNIH2":[],"Gria1":[],"Gria2":[3]}
Channel_names = {0:"Dapi",1:"Cnih2",2:"MAP2",3:"CAMKII"}
scale1 = 2.41 #pixels/uM
scale2 = 7.587 #pixels/uM
channels = [1,3]
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
COLORS = [CB91_Blue,CB91_Green,CB91_Pink,CB91_Purple,CB91_Violet,CB91_Amber]
Lengths = np.asarray([25,50,100,150,200,250])
widths = [5.0,10,15.0]
Analysis_dend_width = 0
def ReadFiles(dir,cells_to_exclude,molecule,width,cell_num):
     cells = range(1,cell_num)
     # print(cells)
     total_d_count = 0;
     for i in channels:
         soma_data['channel_'+str(i)] = {}
         dend_data['channel_'+str(i)] = {}
         for cell in cells:
             # print(cell)
              file2 = folder+"cell_"+str(cell)+"/"+str(2)+"/"+str(width)+"/dend_puncta_channel{0}_{1}.json".format(2,width) 
              # print(file)
              if Path(file).is_file():
                 dend_data_meta['cell_'+str(cell)],soma_data_meta['cell_'+str(cell)] = LoadDend(folder+"cell_"+str(cell)+"/",scale1)
                 total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                 f_d = open(file)
                 dend_data['channel_'+str(i)]['cell_'+str(cell)] = json.load(f_d)
                 f_d.close()
                 f_s = open(folder+"cell_"+str(cell)+"/"+str(i)+"/"+str(width)+"/soma_puncta_channel{0}_{1}.json".format(i,width))
                 # print(f_d,f_s)
                 soma_data['channel_'+str(i)]['cell_'+str(cell)] = json.load(f_s)
                 f_s.close()
              if Path(file2).is_file():
                 dend_data_meta['cell_'+str(cell)],soma_data_meta['cell_'+str(cell)] = LoadDend(folder+"cell_"+str(cell)+"/",scale2)
                 total_d_count += len(dend_data_meta['cell_'+str(cell)].keys())
                 f_d = open(file2)
                 dend_data['channel_1']['cell_'+str(cell)] = json.load(f_d)
                 f_d.close()
                 f_s = open(folder+"cell_"+str(cell)+"/"+str(2)+"/"+str(width)+"/soma_puncta_channel{0}_{1}.json".format(2,width))
                 # print(f_d,f_s)
                 soma_data['channel_1']['cell_'+str(cell)] = json.load(f_s)
                 f_s.close()
                 # print(soma_data['channel_'+str(i)])
              # for idx,file in enumerate(os.listdir(folder+"cell_"+str(cell))):
              #    if file.endswith(".json"):
              #        print(folder+"cell_"+str(cell)+"/"+file)
         # print(soma_data['channel_1'].keys())
         # print("*"*20)
     # PlotSomaDistribution(soma_data)
     # print(soma_data_meta)
     # print("*"*20)
     # print(dend_data_meta)
     # print(soma_data['channel_1'].keys())
     # print(soma_data['channel_3'].keys())
     # print(dend_data['channel_1'].keys())
     # print(dend_data['channel_3'].keys())
     # sd_stats,sum_soma_punctas,sum_soma_cell_punctas,cell_soma_punctas = PlotSomaDistribution(soma_data)
     # dend_stats,sum_dend_punctas,sum_dend_cell_punctas,cell_dend_punctas = PlotSomaDistribution(dend_data)
     bin_size = 2
     soma_bins = np.arange(0,1,bin_size)
     dend_bins = np.arange(0,300,bin_size)
     breakpoint()
     soma_total_count ,soma_cell_count ,soma_unit_count, soma_total_stat, soma_cell_stat, soma_unit_stat, soma_cell_unique_stats = GetPunctaDicts(soma_data,bins=soma_bins)
     dend_total_count ,dend_cell_count ,dend_unit_count, dend_total_stat, dend_cell_stat, dend_unit_stat, dend_cell_unique_stats  = GetPunctaDicts(dend_data,bins=dend_bins)
     c_ls = GetLengthCounts(dend_data_meta)
     # breakpoint()
     
     # breakpoint()
     dend_cell_sum = {}
     # for c in dend_cell_stat.keys():
     #     dend_cell_sum[c] = {}
     #     for cell in dend_cell_stat[c].keys():
     #         dend_cell_sum[c][cell] = np.array([GetMatSum(dend_cell_stat[c][cell])])
     # breakpoint()
     
     # dend_distribution_dict = GetDendriticSpatialDistribution(soma_cell_unique_stats,dend_unit_stat,dend_data_meta) 
     # breakpoint()
     dend_sum_norm_distribution_dict = GetDendriticSpatialDistribution(None,dend_unit_stat,dend_data_meta) 
     
     mean_dist = {}
     std_dist = {}
     sem_dist = {}
     
     mean_sum_norm_dist = {}
     std_sum_norm_dist = {}
     sem_sum_norm_dist = {}
     for ldx,l in enumerate(Lengths):
        l1 = l #str(l)
        mean_dist[l1] = {}
        std_dist[l1] = {}
        sem_dist[l1] = {}
        mean_sum_norm_dist[l1] = {}
        std_sum_norm_dist[l1] = {}
        sem_sum_norm_dist[l1] = {}
        for cdx,c in enumerate(dend_sum_norm_distribution_dict[l].keys()):
            # mean_dist[l][c] = dend_distribution_dict[l][c].mean(axis=0)
            # std_dist[l][c] = dend_distribution_dict[l][c].std(axis=0)
            # sem_dist[l][c] = std_dist[l][c]/np.sqrt(dend_distribution_dict[l][c].shape[0])
            mean_sum_norm_dist[l1][c] = dend_sum_norm_distribution_dict[l][c].mean(axis=0)
            std_sum_norm_dist[l1][c] = dend_sum_norm_distribution_dict[l][c].std(axis=0)
            sem_sum_norm_dist[l1][c] = std_sum_norm_dist[l1][c]/np.sqrt(dend_sum_norm_distribution_dict[l][c].shape[0])
     # breakpoint()
     dist_data = {}
     dist_data['means'] = mean_sum_norm_dist
     dist_data['std'] = std_sum_norm_dist
     dist_data['sem'] = sem_sum_norm_dist
     # breakpoint()
     # with open(folder+"dend_dist_data_bin_size_"+str(bin_size)+".json", 'w') as fp: json.dump(dist_data,fp,indent=4)
     pw = PlottingWidgetmRNA()
     
     op_folder = folder+str(width)+"/"
     pw.CreateFolderRecursive(op_folder)

    
     save_it = 1
     fractions,cell_total_stats = GetSomaticDendriticFraction(soma_cell_unique_stats,dend_cell_unique_stats)
     lab1= ['{0}-Soma'.format(Channel_names[channels[0]]),'{0}-Dendrite'.format(Channel_names[channels[0]]),'{0}-Soma'.format(Channel_names[channels[1]]),'{0}-Dendrite'.format(Channel_names[channels[1]])]
     x_lab,y_lab = ["mRNA",'fraction']
     title = "Fraction of mRNA in Soma Vs. dendrites";
     stat_no = 2;
     data_to_show = np.asarray([fractions['channel_1'][:,0,0],fractions['channel_1'][:,1,0],fractions['channel_3'][:,0,0],fractions['channel_3'][:,1,0]])
     # breakpoint()
     pw.PlotCellFraction(data_to_show,lab1,x_lab,y_lab,title,op_folder+"soma_vs_dend_fractions_area_x_mean",save_it = save_it)
     # x_lab,y_lab = ["mRNA",'fraction']
     # pw.PlotCellFraction(data_to_show[2],data_to_show[3],lab1[2],lab1[3],x_lab,y_lab,title,folder+"soma_vs_dend_fractions_CaMKII",save_it = save_it)
     # breakpoint()
     breakpoint()
     soma_dend_ratio = GetCellWiseRatio(soma_cell_unique_stats,dend_cell_unique_stats)
     lab1,lab2 = [Channel_names[channels[0]],Channel_names[channels[1]]]
     x_lab,y_lab = ["mRNAs",'Ratio (Dendrite/Soma)']
     title = "Dendritic vs. Somatic mRNA ratio"
     file_prefix = "dend_vs_soma_ratio_"
     pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,0],soma_dend_ratio['channel_3'][:,0],lab1,lab2,x_lab,y_lab ,title+r" ( copy number = $A*\mu$ )",op_folder+file_prefix+"area_x_mean_10_um",save_it = save_it)
     pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,1],soma_dend_ratio['channel_3'][:,1],lab1,lab2,x_lab,y_lab ,title+r" ( copy number = $A*\tilde{x}$ )",op_folder+file_prefix+"area_x_med_10_um",save_it = save_it)
     pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,2],soma_dend_ratio['channel_3'][:,2],lab1,lab2,x_lab,y_lab ,title+r" ( copy number = $\mu$ )",op_folder+file_prefix+"mean_10_um",save_it = save_it)
     pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,3],soma_dend_ratio['channel_3'][:,3],lab1,lab2,x_lab,y_lab ,title+r" ( copy number = $\tilde{x}$ )",op_folder+file_prefix+"med_10_um",save_it = save_it)
     # pw.ViolinPlotStats(soma_dend_ratio['channel_1'][:,4],soma_dend_ratio['channel_3'][:,4],lab1,lab2,x_lab,y_lab ,title+r" ( copy number = $A$ )",folder+file_prefix+"dend_vs_soma_ratio_area_10_um",save_it = save_it)
     
     
     
     
     # pw.PlottwoStats(sd_stats[1][:,2],sd_stats[1][:,3],'mean vs median',"mean","median","coorealtion between mean and median","./sasd",save_it = 0)
    # pw.PlottwoStats(sd_stats[:,2],sd_stats[:,4],'mean vs area',"mean","area","coorealtion between mean and area","./sasd",save_it = 0)
    # pw.PlottwoStats(sd_stats[:,3],sd_stats[:,4],'median vs area',"median","area","coorealtion between median and area","./sasd",save_it = 0)
     
     # binned_dend_punctas,binned_dend_punctas_cellwise,count_bins = GetDendriticSpatialDistribution(sum_soma_cell_punctas,sum_dend_punctas,bins )
     # print(binned_dend_punctas['channel_1'][:,2])
     # print(binned_dend_punctas,binned_dend_punctas_cellwise,count_bins)
     # pw.PlotBinnedStats(bins,binned_dend_punctas['channel_1'][:][:,3].mean(),binned_dend_punctas['channel_3'][:,3].mean(),'GriA1','CaMKII','Dendritic distance','mRNA number',r"dendritic mRNA( copy number = A*$\mu$ )","./sasd",save_it = 0)
     # print(dend_data)


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
   
    for channel in dend_d.keys():
        for cell in dend_d[channel].keys():
            for unit in dend_d[channel][cell].keys():
                # m_unit = "Dendrite"+unit.split("_")[1]
                for l in Lengths:
                    # breakpoint()
                    print("meta data for = {0}_{1} = {2}".format(cell,unit,dend_meta_data[cell][unit]))
                    if dend_meta_data[cell][unit] >= l:
                        div = GetMatSum(dend_d[channel][cell][unit])[0:5]
                        if (div==np.zeros(div.shape)).any():
                            print("zero sum found. No puncta found in the dendrite {} of dendrite {} of lenght = {}".format(unit,cell,dend_meta_data[cell][unit]))
                        if soma_d != None:
                            print("som is provided")
                            div = soma_d[channel][cell][0][0:5]
                        soma_norm = np.divide(dend_d[channel][cell][unit][:,0:5],div)
                        # print(soma_norm)
                        if not np.isnan(soma_norm).any():
                            dendritic_dist[l][channel].append(soma_norm)
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
    cell_unique_stats = {}
    
    for channel in od.keys():
        total_count[channel] = 0
        cell_count[channel] = {}
        unit_count[channel] = {}
        total_stat[channel] = []
        cell_stat[channel] = {}
        unit_stat[channel] = {}
        cell_unique_stats[channel] = {}
        for cell in od[channel].keys():
            cell_count[channel][cell] = 0
            unit_count[channel][cell] = {}
            
            cell_stat[channel][cell] = []
            unit_stat[channel][cell] = {}
            
            unique_punctas = []
            for unit in od[channel][cell].keys():
                all_punctas = od[channel][cell][unit]
                # print(len(all_punctas),end=" ")
                if unique_punctas == []:
                    unique_punctas = all_punctas
                else:
                    unique_punctas = GetUniqueRows(np.concatenate((all_punctas, unique_punctas)))
                total_count[channel] += 1
                cell_count[channel][cell] += 1
                unit_count[channel][cell][unit] = len(all_punctas)
                all_stats = GetAllPunctastat(all_punctas)
                binned_sum = BinnedSum(all_stats,bins,-1,"{0}_{1}_{2}".format(channel,cell,unit)) # binning based on lenght  
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
            all_unique_stats = GetAllPunctastat(unique_punctas)
            cell_unique_stats[channel][cell] = BinnedSum(all_stats,bins,-1)
                # print(max_length)
    return total_count ,cell_count ,unit_count, total_stat,cell_stat,unit_stat, cell_unique_stats 
            
def BinnedSum(arr,bins,num_col=-1,name= None):
    # print("binned_sum for =",name)
    if len(arr.shape) == 2:
        rr,cc = arr.shape 
        binned_sum = np.zeros((len(bins),cc))
        arr = SortPunctas(arr,num_col)
        # max_lenght = arr[-1,-1]
        digitized = bins.searchsorted(arr[:,num_col])
        # if len(digitized) > 1:
        #     digitized[0] = digitized[1]
        # if name=="channel_1_cell_2_dendrite_0":
        #     breakpoint()
        for c in range(0,cc):
            binned_sum[:,c] = np.bincount(digitized, weights=arr[:,c], minlength=len(bins))
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
    return [stat1,stat2,stat3,stat4,stat5,dfs]

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
        if('AllDendrite' in x):
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
    
    def PlotBinnedStats(self,x,mean_1,mean_2,std_1,std_2,lab1,lab2,xlab,ylab,title_string,file_name,bin_size,width=10,height=8,fsize=16,save_it = 1,fit_exp =1,auto_close=1):
        fig, ax = plt.subplots(figsize=(width, height))
        plt.rc('font', **{'family':'serif','serif':['Palatino']})
        plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        ax.plot(x,np.zeros(x.shape),'k--',label='=0')
        ax.errorbar(x+bin_size/2,mean_1,std_1,label=lab1,color=COLORS[0],marker='o')
        ax.errorbar(x+bin_size/2,mean_2,std_2,label=lab2,color=COLORS[1],marker='o')
        ax.set_xlabel(xlab,fontsize=fsize)
        ax.set_ylabel(ylab,fontsize=fsize)
        plt.title(title_string,fontsize=fsize)
        
        folder = "."
        
        if fit_exp == 1:
            y1_fit, r1_squared = ExpFit(x,mean_1,1,-1)
            y2_fit, r2_squared = ExpFit(x,mean_2,1,-1)
            ax.plot(x+bin_size/2,y1_fit,'^--',c=COLORS[0],label=lab1+"-fit, $r^2$ =%0.2f" %(r1_squared))
            ax.plot(x+bin_size/2,y2_fit,'^--',c=COLORS[1],label=lab2+"-fit, $r^2$ =%0.2f" %(r2_squared))
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
    ReadFiles(folder,Excluded_cells[mRNA],mRNA,widths[Analysis_dend_width_cells[mRNA])
    