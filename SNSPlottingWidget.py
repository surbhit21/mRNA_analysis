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
from Fonkeu_et_al_2019 import GetmRNADist
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
            # breakpoint()
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
                    yi_fit,chi_squ = ExpFitWithMinimize(exp_method,xs[i],means[i],stds[i,:],0,+1,labs[i])
                    # breakpoint()
                    ax.plot(xs[i],yi_fit,marker='None',c=color[i],label=labs[i]+r"-fit,$\chi^2$ = {:.2f}".format(chi_squ))
            
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
    
    def PlotBinnedStats1p(self,xs,means,stds,MAP2_norm,labs,xlab,ylab,ylab_norm,color,title_string,file_name,bin_size,save_it = 1,fit_exp =1,auto_close=1,in_set=1,lw=3.0,set_axis_label=1,exp_method = "NormE"):
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
            # breakpoint()
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i],means[i],stds[i],label=labs[i],color=color[i],marker='d',linestyle='dashed' )
               
               if in_set==1:
                   ax2.plot(xs[i],MAP2_norm[i],color=color[i],marker='o',markersize=4,linestyle='dashed')
            if set_axis_label == 1:
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)
            
            folder = "."
            
            if fit_exp == 1:
                for i in range(xs.shape[0]):
                    # breakpoint()
                    yi_fit, ri_squared,chi_squ = ExpFitWithMinimize(exp_method,xs[i],means[i],stds[i,:],0,+1,labs[i])
                    ax.plot(xs[i],yi_fit,marker='None',c=color[i],label=labs[i]+"-fit")
                    fonkey_fit_y = GetmRNADist(xs[i],0.0018,0.0018)
                    fonkey_fit_y1 = GetmRNADist(xs[i],0.0008,0.0018)
                    fonkey_fit_y2 = GetmRNADist(xs[i],0.0028,0.0018)
                    # 0.0008,0.0028
                    ax.plot(xs[i],fonkey_fit_y,marker='None',c='k',label="Fonkeu et al. (2019)")
                    ax.plot(xs[i],fonkey_fit_y1,'k--')
                    ax.plot(xs[i],fonkey_fit_y2,'k--')
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
            
            
    def PlotMAP2(self,xs,means,stds,labs,xlab,ylab,title_string,file_name,bin_size,width=10,height=8,fsize=16,save_it = 1,fit_exp =1,auto_close=1):
        fig, ax = plt.subplots(figsize=(width, height))
        # plt.rc('font', **{'family':'serif','serif':['Palatino']})
        # plt.rc('text', usetex=True)
        # ax.scatter(x_data1,y_data1,label=lab1)
        # ax.scatter(x_data2,y_data2,label=lab2)
        ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0')
        if xs.shape[0] == means.shape[0]:
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i]+bin_size/2,means[i],stds[i],label=labs[i],color='k',marker='o')
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
                plt.savefig(file_name+".pdf",dpi=150)
                print("saved figures to: {%s/%s}" %(folder, file_name))
            else:
                print("Plots not saved")
            plt.show()
            if auto_close == 1:
                plt.close()
        else:
            print("Not same number of xs and means")

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
        num_plots = int(fractions.shape[0]/groups)
        pos = np.linspace(1,2,groups)
        # breakpoint()
        x_points = []
        pairs = [] 
        x_tics = []

        for i in range(num_plots):
            bp1 = ax.boxplot(fractions[groups*i],widths = 0.5,positions=[i*(groups+1)+pos[0]],showfliers=False,labels=[compartment[0]])
            bp2 = ax.boxplot(fractions[groups*i+groups-1],widths = 0.5,positions=[i*(groups+1)+pos[1]],showfliers=False,labels=[compartment[1]],patch_artist=True, )
            x_points.append(i*(groups+1)+pos)
            x_points.append([i*(groups+1)+pos[1],(num_plots-1)*(groups+1)+pos[1]])
            x_tics.append(i*(groups+1)+pos.mean())
            pairs.append([i*groups+1,i*groups+2])
            pairs.append([i*groups+2,(num_plots-1)*groups+2])
            self.setBoxColors(bp1,color[i],1)
            self.setBoxColors(bp2,color[i],1,True)
        
        # plt.plot([], c=color[0], label=compartment[0])
        # plt.plot([], c=color[1], label=compartment[1])
        plt.xticks(x_tics, lab)
        if fractions.shape[0]>2:
            p_values = sp.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        else:
            p_val = ttest_ind(fractions[0],fractions[1]).pvalue
            max_ind = np.asarray(pairs).max()
            p_values = np.ones((max_ind+1,max_ind+1))*p_val
        
        x_points = np.asarray(x_points).flatten()
        # breakpoint()
        pairs = np.array(pairs)
        # plt.legend(loc="upper right")
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
            fac = np.abs(x2-x1)*0.04
            plt.plot([x1,x1, x2,x2], [y+fac, y+h+fac, y+h+fac, y+fac], lw=1.5, c=color)
            plt.text((x1+ x2)*0.5,y+h+fac,txt, ha=ha, va=va, color=color)
        
    def setBoxColors(self,bp,c,a=1,flipped=False):
        setp(bp['boxes'], color=c,alpha=a)
        setp(bp['caps'], color=c,alpha=a )
        setp(bp['caps'], color=c ,alpha=a)
        setp(bp['whiskers'], color=c ,alpha=a)
        setp(bp['whiskers'], color=c ,alpha=a)
        setp(bp['fliers'], color=c ,alpha=a)
        setp(bp['fliers'], color=c ,alpha=a)
        if flipped == True:
            setp(bp['medians'], color='w',alpha=a )
        else:
            setp(bp['medians'], color=c,alpha=a )
        # setp(bp['mean'], color=c ,alpha=0.8)
        # if fill == True:
        #     setp(bp['boxes'], facecolor=c)

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def normExponential(x, params):
    b = params['b'].value
    return np.exp(b*x)

def oneExponential(x,params):
    a = params['a'].value
    b = params['b'].value
    return a*np.exp(b*x)

def twoExponential(x,params):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    return a*np.exp(b*x) + c*np.exp(d*x)

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
    print("chi-squared = ",chi_squ)
    return y_fit,r_squared,chi_squ

def ExpFitWithMinimize(ftype,xdata,ydata,sigmas,Fidx,Lidx,molecule):
    
    """
        Fit a function to a given distribution using lmfit minimize method
        
    """
    fit_paramas = Parameters()
    # np.random.seed(2022)
    exp_min = -2.0
    exp_max = 0
    pref_min = 0
    pref_max = 200

    
    
    
    if ftype == "NormE":
        b_init = np.random.uniform(exp_min,exp_max)
        fit_paramas.add('b',b_init,min=exp_min,max=exp_max)
        
        residuals = Residual(fit_paramas,normExponential,xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual,params=fit_paramas,method='leastsq',args=(normExponential,xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = normExponential(xdata[Fidx:],out2.params)
        
        # breakpoint()
        
    elif ftype == "1E":
        
        a_init = np.random.uniform(pref_min,pref_max)
        b_init = np.random.uniform(exp_min,exp_max)
        fit_paramas.add('a',a_init,min=pref_min,max=pref_max)
        fit_paramas.add('b',b_init,min=exp_min,max=exp_max)
        residuals = Residual(fit_paramas,oneExponential,xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual,params=fit_paramas,method='leastsq',args=(oneExponential,xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = oneExponential(xdata[Fidx:],out2.params)
       
    elif ftype == "2E":
        a_init = np.random.uniform(pref_min,pref_max)
        b_init = np.random.uniform(exp_min,exp_max)
        c_init = np.random.uniform(pref_min,pref_max)
        d_init = np.random.uniform(exp_min,exp_max)
        fit_paramas.add('a',a_init,min=pref_min,max=pref_max)
        fit_paramas.add('b',b_init,min=exp_min,max=exp_max)
        fit_paramas.add('c',c_init,min=pref_min,max=pref_max)
        fit_paramas.add('d',d_init,min=exp_min,max=exp_max)
        residuals = Residual(fit_paramas,twoExponential,xdata[Fidx:], ydata[Fidx:])
        out2 = minimize(Residual,params=fit_paramas,method='leastsq',args=(twoExponential,xdata[Fidx:], ydata[Fidx:]))
        report_fit(out2.params)
        y_fit = twoExponential(xdata[Fidx:],out2.params)
    else:
        raise NotImplementedError("ftype: {} not implemented, contact author or define it yourself".format(ftype))
        
    # print("fitted "+ ftype, popt)
    
    chi_squ = out2.chisqr
    
    return y_fit,chi_squ

def Residual(paras,fun,x,data):
    expected_vals = fun(x,paras)
    res = expected_vals - data
    return res

def ChiSq(yd,y_fit,sigmas):
    nzs = np.nonzero(sigmas)
    print(nzs)
    r_yd = np.take(yd,nzs)
    r_yf = np.take(y_fit,nzs)
    r_sgs = np.take(sigmas,nzs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals/r_sgs)**2)
    return chi_squ