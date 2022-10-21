#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:22:33 2022

@author: surbhitwagle
"""

import pandas as pd
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import lmfit
from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
from operator import add
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal,ttest_ind
import seaborn as sns
import scikit_posthocs as sp

"""
Parent class for All plotting stuff
"""
class SNSPlottingWidget():
    def __init__(self,fsize=16,tsize=25,fam='Source Code Pro',pixelden = 100,lw=3.0,width=10,height=8):
        rc('font',
            family=fam,
            size=fsize)
        # rc('text', 
        #     # linewidth=lw,
        #     fsize=10,
        #     )
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
        """
            function creates folders upto  a path recursively
            arguments : folderpath
        """
        Path(folder).mkdir(parents=True, exist_ok=True)
    
    def SaveFigures(self,filename,ext_list = [".png",".svg",".pdf"],dpi=300):
        """
            function to save figures
            required arguments:
                filename
        """
        for ext in ext_list:
            plt.savefig(filename+ext,dpi=dpi)
        
    def PlotBinnedStats(self,xs,means,stds,MAP2_norm,labs,xlab,ylab,ylab_norm,color,title_string,file_name,bin_size,save_it = 1,fit_exp =1,auto_close=1,in_set=1,lw=3.0,set_axis_label=1,exp_method = "NormE"):
        """
            function to plot binned statistics
            
        """
        fig, ax = plt.subplots()
        # ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0' ,markersize=4*lw)
        if in_set ==1:
            left, bottom, width, height = [0.35, 0.7, 0.2, 0.2]
            ax2 = fig.add_axes([left, bottom, width, height])
            ax2.set_ylabel(ylab_norm)
        if xs.shape[0] == means.shape[0]:
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i],means[i],stds[i],label=labs[i],color=color[i],marker='d',linestyle='None' )
               
               if in_set==1:
                   ax2.plot(xs[i],MAP2_norm[i],color=color[i],marker='o',markersize=4,linestyle='dashed')
            if set_axis_label == 1:
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            
            folder = "."
            
            if fit_exp == 1:
                for i in range(xs.shape[0]):
                    # breakpoint()
                    yi_fit, ri_squared,chi_squ = ExpFit(exp_method,xs[i],means[i],stds[i,:],0,+1,labs[i])
                    ax.plot(xs[i],yi_fit,marker='None',c=color[i],label=labs[i]+"-fit")
            
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
            
    
    def PlotFittedCurves(self,xs,data1,data2,labs,xlab,ylab,ylab_norm,color,title_string,file_name,bin_size,save_it = 1,lw=3.0,set_axis_label=1,exp_method = "NormE"):
        """
            function to plot binned statistics
            
        """
        fig, ax = plt.subplots()
        # ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0' ,markersize=4*lw)
       
        if xs.shape[0] == data1.shape[1] and xs.shape[0] == data2.shape[1]:
               
            if set_axis_label == 1:
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            
            folder = "."
            
            for i in range(0,data1.shape[0]):
                # breakpoint()
                yi_fit, ri_squared,chi_squ = ExpFit(exp_method,xs,data1[i],np.zeros(data1[i].shape),0,+1,labs[0])
                ax.plot(xs,yi_fit,marker='None',c=color[0],alpha=0.1)
            
            # plt.show()
            # breakpoint()
           
            
            for i in range(0,data2.shape[0]):
                # breakpoint()
                yi_fit, ri_squared,chi_squ = ExpFit(exp_method,xs,data2[i],np.zeros(data2[i].shape),0,+1,labs[1])
                ax.plot(xs,yi_fit,marker='None',c=color[1],alpha=0.1)
            
            # plt.show()
            means1 = data1.mean(axis = 0)
            yi_fit, ri_squared,chi_squ = ExpFit(exp_method,xs,means1,np.zeros(means1.shape),0,+1,labs[0])
            ax.plot(xs,yi_fit,marker='None',c=color[0],alpha=1)
            means2 = data2.mean(axis = 0)
            yi_fit, ri_squared,chi_squ = ExpFit(exp_method,xs,means2,np.zeros(means2.shape),0,+1,labs[1])
            ax.plot(xs,yi_fit,marker='None',c=color[1],alpha=1)
            
            fig.tight_layout()
            ax.legend()
            if save_it == 1: 
                self.SaveFigures(file_name)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
          
        else:
            print("Not same length of xs and means",xs.shape,data.shape)
    def PlotCellFraction(self,fractions,lab,compartment,xlab,ylab,color,title_string,file_name,molecules,groups=2,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        # breakpoint()
        num_plots = int(fractions.shape[0]/groups)
        pos = np.linspace(1,2,groups)
        # breakpoint()
        x_points = []
        pairs = [] 
        x_tics = []

        for i in range(num_plots):
            # breakpoint()
            print(i)
            bp1 = ax.boxplot(fractions[groups*i],widths = 0.5,positions=[i*(groups+1)+pos[0]],showfliers=False)
            bp2 = ax.boxplot(fractions[groups*i+groups-1],widths = 0.5,positions=[i*(groups+1)+pos[1]],showfliers=False)
            x_points.append(i*(groups+1)+pos)
            x_points.append([i*(groups+1)+pos[1],(num_plots-1)*(groups+1)+pos[1]])
            x_tics.append(i*(groups+1)+pos.mean())
            pairs.append([i*groups+1,i*groups+2])
            pairs.append([i*groups+2,(num_plots-1)*groups+2])
            # means = []
            # stds = []
            # for k in range(groups):
            #     means.append(fractions[groups*i+k-1].mean())
            #     stds.append(fractions[groups*i+k-1].std())
            # # sp = sns.swarmplot(y=fractions[2*i:2*i+2],x=i*2+pos)
            # for j, line in enumerate(bp['means']):
            #     x, y = line.get_xydata()[1]
            #     text = r'${:.2f}(\pm {:.2f})$'.format(means[j], stds[j])
            #     ax.annotate(text, xy=(x, y))
            self.setBoxColors(bp1,color[1])
            self.setBoxColors(bp2,color[2])
        # breakpoint()
        plt.plot([], c=color[1], label=compartment[0])
        plt.plot([], c=color[2], label=compartment[1])
        plt.xticks(x_tics, lab)
        if fractions.shape[0]>2:
            p_values = sp.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        else:
            # breakpoint()
            p_val = ttest_ind(fractions[0],fractions[1]).pvalue
            max_ind = np.asarray(pairs).max()
            p_values = np.ones((max_ind+1,max_ind+1))*p_val
        
        x_points = np.asarray(x_points).flatten()
        # breakpoint()
        pairs = np.array(pairs)
        plt.legend(loc="upper right")
        y_max = 1
        for idx,pair in enumerate(pairs):
            txt = ''
            print(p_values[pair[0]][pair[1]],pair)
            if p_values[pair[0]][pair[1]] <= 0.05:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.01:
                txt += '*'
            if p_values[pair[0]][pair[1]] <= 0.001:
                txt += '*'
            y_max += 0.02
            self.AnnotateText(ax,x_points[idx*2],x_points[idx*2+1],y_max,0.01,txt,'k')
        
        if set_axis_label == 1:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        # plt.ylim([0,1.3])
        plt.title(title_string)
        # plt.tick_params(
        #     axis='x',          # changes apply to the x-axis
        #     which='both',      # both major and minor ticks are affected
        #     bottom=False,      # ticks along the bottom edge are off
        #     top=False,         # ticks along the top edge are off
        #     labelbottom=False) # labels along the bottom edge are off
        # means = []
        # stds = []
        # num_means =fractions.shape[0]
        # for i in range(num_means):
        #     means.append(fractions[i].mean())
        #     stds.append(fractions[i].std())
        # # breakpoint()
        # for j, line in enumerate(bp['means']):
        #     x, y = line.get_xydata()[1]
        #     text = r'${:.2f}(\pm {:.2f})$'.format(means[j], stds[j])
        #     ax.annotate(text, xy=(x, y))
        # # for j, line in enumerate(bp2['means']):
        # #     x, y = line.get_xydata()[1]
        # #     text = r' ${:.2f}(\pm {:.2f})$'.format(means[2+j], stds[2+j])
        # #     ax.annotate(text, xy=(x, y))
        fig.tight_layout()
        folder = "."
        if save_it == 1:
            self.SaveFigures(file_name)
            print("saved figures to: {%s/%s}" %(folder, file_name))
        else:
            print("Plots not saved")
        plt.show()
    
    def ViolinPlotStats(self,data,labs,xlab,ylab,colors,title_string,file_name,molecules,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        
        # ax.set_axis_style([lab1,lab2])
        # labels=[lab1+"\n N=%d cells"%(data1.shape[0]),lab2+"\n N=%d cells"%(data2.shape[0])]
        num_plots = data.shape[0]
        x = np.arange(len(labs))
        ax.set_xticks(x)
        ax.set_xticklabels(labs)
        xs = np.zeros(data.shape)
        df = pd.DataFrame()
        for i in range(0,num_plots):
            # breakpoint()
            violin_parts = ax.violinplot(data[i],positions=[x[i]],showextrema=False,showmeans = False,showmedians = False,points=len(data[i]))
            
            
            quartile1, median, quartile3 = np.percentile(data[i], [25, 50, 75], axis=0)
            whiskers_min,whiskers_max = data[i].min(),data[i].max()
            ax.scatter(x[i],median,marker='s',color='white', s=10, zorder=100)
            ax.vlines(x[i], quartile1, quartile3, color='k', linestyle='-', lw=6, zorder=99)
            # ax.vlines(x[i], whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
            # violin_parts['cmeans'].set_edgecolor('red')
            # violin_parts['cmedians'].set_edgecolor('green')
            for pc in violin_parts['bodies']:
                pc.set_facecolor(colors[i])
                pc.set_edgecolor(colors[i])
                pc.set_alpha(0.3)
            df1 = pd.DataFrame({labs[i]:data[i]})
            df = pd.concat([df,df1],axis=1)
            # violin_parts['cbars'].set_edgecolor(None)
            # violin_parts['cbars'].set_alpha(0)
        # for i in range(num_plots):
        #     ax.scatter(xs,data,color=colors[i],alpha=0.6)
        
        # breakpoint()
        ax = sns.swarmplot(data=df,palette=colors)
        # print(violin_parts['bodies'])
        
        
        
        if set_axis_label == 1:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
        # stat,p_val = kruskal(data1, data2)
        plt.title(title_string)#+"\n stat = %0.2e, p-value = %0.2e")%(stat,p_val),fontsize=fsize)
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

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def normExponential(x, b):
    return np.exp(b*x)

def oneExponential(x,a,b):
    return a*np.exp(b*x)
def twoExponential(x,a,b,c,d):
    return oneExponential(x,a,b) + oneExponential(x,c,d)

def ExpFit(ftype,xdata,ydata,sigmas,Fidx,Lidx,molecule):
    
    """
        Fit a function to a given distribution
        
    """
    if ftype == "NormE":
        param_bounds=([-np.inf],[np.inf])
        popt, pcov = curve_fit(normExponential, xdata[Fidx:], ydata[Fidx:],bounds =param_bounds, maxfev=5000)
        y_fit = normExponential(xdata, *popt)
    elif ftype == "1E":
        param_bounds=([-np.inf,-np.inf],[+np.inf,0])
        popt, pcov = curve_fit(oneExponential, xdata[Fidx:], ydata[Fidx:],bounds =param_bounds)
        y_fit = oneExponential(xdata, *popt)
    elif ftype == "2E":
        param_bounds=([-np.inf,-np.inf,-np.inf,-np.inf],[+np.inf,0,+np.inf,0])
        popt, pcov = curve_fit(twoExponential, xdata[Fidx:], ydata[Fidx:],bounds =param_bounds)
        y_fit = twoExponential(xdata, *popt)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))
        
    print("fitted "+ ftype, popt)
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