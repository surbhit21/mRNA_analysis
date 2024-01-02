#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:22:33 2022

@author: surbhitwagle
"""
import matplotlib.pyplot as plt
from CNIH2_protein_model_fitting import FitModelProtein
from Fonkeu_et_al_2019 import GetmRNADist
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import rc
from matplotlib.pyplot import setp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mRNA_model_fitting import *
import numpy as np
from pathlib import Path
from scipy.stats import ks_2samp, kruskal,ttest_ind,spearmanr,pearsonr
import seaborn as sns
import scikit_posthocs as sp


"""
Parent class for All plotting stuff
"""
class SNSPlottingWidget():
    def __init__(self,fsize=18,tsize=25,fam='serif',pixelden = 100,lw=3.0,width=8,height=6):
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
            titlesize=fsize,
            titleweight='regular',
            labelsize=fsize,
            )
    
        rc('legend',
           loc='best',
           frameon=False,
            fontsize=fsize)
        rc('legend.frameon')
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
        self.fsize = fsize
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
        plt.tight_layout()
        for ext in ext_list:
            plt.savefig(filename+ext,dpi=dpi)

    def getclosestindex(self,arr, per):
        total = np.sum(arr)
        perc = total*per/100

        for i in range(0,len(arr)):
            perc -= arr[i]
            if perc <= 0:
                return i
        return len(arr)-1
    def PlotBinnedStats(self,xs,means,stds,MAP2_norm,labs,xlab,ylab,ylab_norm,color,title_string,file_name,bin_size,save_it = 1,fit_exp =1,auto_close=1,in_set=1,lw=3.0,set_axis_label=1,exp_method = "2E"):
        """
            function to plot binned statistics
            
        """
        fig, ax = plt.subplots()
        ax.tick_params(axis='both', which='major', labelsize=self.fsize)
        # ax.set_ylim([0,1.5])
        ax.spines[['right', 'top']].set_visible(False)
        # ax.plot(xs[0],np.zeros(xs[0].shape),'k--',label='=0' ,markersize=4*lw)
        # if in_set ==1:
        #     left, bottom, width, height = [0.25, 0.75, 0.2, 0.2]
        #     # ax2 = fig.add_axes([left, bottom, width, height])
        #     # ax2.set_ylabel(ylab_norm)
        #     ax2 = ax.inset_axes(bounds = [left,bottom,width,height],zorder=4)#, width=3, height=2, loc="upper center")
        #     ax2.set_ylabel(ylab_norm)
        #     ax2.yaxis.tick_right()
        if xs.shape[0] == means.shape[0]:
            # breakpoint()

            # if in_set == 1:
            #     for i in range(xs.shape[0]):
            #         ax2.plot(xs[i], MAP2_norm[i] , color=color[i], marker='o', markersize=4, linestyle='dashed')
            for i in range(xs.shape[0]):
               ax.errorbar(xs[i],means[i],stds[i],label=labs[i],color=color[i],marker='d',linestyle='None' ,zorder=6)
               # per = 75
               # i_close = self.getclosestindex(means[i],per)
               # ax.vlines(x=xs[i,i_close],ymin=0,ymax=(means+stds).max(),
               #           color=color[i],linestyle="-")
               # print(i_close)
            if set_axis_label == 1:
                ax.set_xlabel(xlab)
                ax.set_ylabel(ylab)

            folder = "."
            # breakpoint()
            if fit_exp == 1:
                left, bottom, width, height = [0.2,0.75,0.2,0.2]
                tics = np.arange(0,1.2,0.5)
                ax3 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)
                ax3.set_title("Normalized \n fits")
                # ax3.yaxis.tick_right()
                ax3.set_ylim([0,1.1])
                ax3.set_yticks(tics)
                for i in range(xs.shape[0]):
                    # breakpoint()
                    yi_fit,chi_squ = ExpFitWithMinimize(exp_method,xs[i],means[i],stds[i,:],0,+1,labs[i])
                    # breakpoint()
                    ax.plot(xs[i],yi_fit,marker='None',c=color[i],label=labs[i]+r"-fit,$\chi^2_\nu$ = {:.2f}".format(chi_squ))
                    ax3.plot(xs[i],yi_fit/yi_fit[0],c=color[i])
            elif fit_exp == 2:
                j_r_factors = [1,10]
                left, bottom, width, height = [0.3, 0.75, 0.2, 0.2]
                ax3 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)
                ax3.set_ylabel("Normalized \n fits")
                ax3.yaxis.tick_right()
                ax3.set_ylim([0, 1.1])
                tics = np.arange(0, 1.2, 0.5)
                ax3.set_yticks(tics)
                for i in range(xs.shape[0]):
                    # breakpoint()
                    print("fitting Fonkeu model for {}".format(labs[i]))
                    x1,yi_fit,chi_squ,paras,mini,out2 = FitModel( xs[i], means[i],stds[i,:],j_r_fcator=j_r_factors[i])
                    breakpoint()
                    # chi_squ = ChiSq(means[i],yi_fit,stds[i])
                    ax.plot(x1, yi_fit, marker='None', c=color[i],
                            label=labs[i] + r"-fit,$\chi^2$ = {:.2f}".format(chi_squ))
                    ax3.plot(xs[i],yi_fit/yi_fit[0],c=color[i])
            elif fit_exp == 3:
                # left, bottom, width, height = [0.3, 0.75, 0.2, 0.2]
                # ax3 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)
                # ax3.set_ylabel("Normalized \n fits")
                # ax3.yaxis.tick_right()
                # ax3.set_ylim([0, 1.1])
                # tics = np.arange(0, 1.2, 0.5)
                # ax3.set_yticks(tics)
                for i in range(xs.shape[0]):
                    # breakpoint()
                    print("fitting Fonkeu model for {}".format(labs[i]))
                    x1,yi_fit,chi_squ,paras,mini,out2 = FitModelProtein( xs[i], means[i],stds[i,:])
                    # breakpoint()
                    # chi_squ = ChiSq(means[i],yi_fit,stds[i])
                    ax.plot(x1, yi_fit, marker='None', c=color[i],
                            label=labs[i] + r"-fit,$\chi^2$ = {:.2f}".format(chi_squ))
                    # ax3.plot(xs[i],yi_fit/yi_fit[0],c=color[i])

            if in_set==1:
                ax.set_ylim([-2,means.max()+10])
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

    def PlotBinnedStats1p(self,xs,means,stds,MAP2_norm,labs,xlab,ylab,ylab_norm,color,title_string,file_name,bin_size,save_it = 1,fit_exp =1,auto_close=1,in_set=1,lw=3.0,set_axis_label=1,exp_method = "2E"):
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
            print("Not same length of xs and means",xs.shape,data1.shape)
    def PlotCellFraction(self,fractions,lab,compartment,xlab,ylab,color,title_string,file_name,molecules,groups=2,save_it = 1,set_axis_label=1):
        fig, ax = plt.subplots()
        ax.spines[['right', 'top']].set_visible(False)
        # y_range = [0.0,0.2,0.4,0.6,0.8,1]
        # ax.set_yticklabels(y_range)
        # ax.set_ylim([0, 1])
        # ax.spines.left.set_bounds((-0.02,1))
        num_plots = int(fractions.shape[0]/groups)
        pos = np.linspace(1,2,groups)
        # breakpoint()
        x_points = []
        pairs = [] 
        x_tics = []
        ax.tick_params(axis='both', which='major', labelsize=self.fsize)
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
        print(pairs)
        # plt.plot([], c=color[0], label=compartment[0])
        # plt.plot([], c=color[1], label=compartment[1])
        plt.xticks(x_tics, lab)
        if fractions.shape[0]>2:
            p_values = sp.posthoc_dunn(fractions, p_adjust = 'bonferroni')
        else:
            p_val = kruskal(fractions[0],fractions[1]).pvalue
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
        ax.spines[['right', 'top']].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=self.fsize)
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
        # fig.tight_layout()
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
            trans = ax.get_xaxis_transform()
            plt.hlines(xmin=x1,xmax=x2, y = y+fac, lw=1.5, color=color)#, transform=trans)
            plt.text((x1+ x2)*0.5,y+h+fac,txt, ha=ha, va=va, color=color)#,transform=trans )
        
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

    def CorrelationCalculationAndPlotting(self,data_to_show, localization, compartment, x_param, y_param, y_bg_param,op_file,
                                          f_scale=1,save_it = 0):
        fig, ax = plt.subplots(figsize=(16, 6), ncols=2, nrows=1)
        sns.set(font_scale=f_scale)
        sns.scatterplot(x=x_param, y=y_param, data=data_to_show, ax=ax[0], hue='compartment').set(
            title="Area vs  {} intensity in {} rois".format(y_param, localization))
        sns.scatterplot(x=x_param, y=y_bg_param, data=data_to_show, ax=ax[1], hue='compartment').set(
            title="Area vs  {}  intensity in {} background rois".format(y_param, localization))
        sns.regplot(x=x_param, y=y_param, data=data_to_show[data_to_show["compartment"] == compartment[0]], ax=ax[0])
        sns.regplot(x=x_param, y=y_param, data=data_to_show[data_to_show["compartment"] == compartment[1]], ax=ax[0])
        sns.regplot(x=x_param, y=y_bg_param, data=data_to_show[data_to_show["compartment"] == compartment[0]], ax=ax[1])
        sns.regplot(x=x_param, y=y_bg_param, data=data_to_show[data_to_show["compartment"] == compartment[1]], ax=ax[1])
        fig.tight_layout()
        corr1 = pearsonr(data_to_show[data_to_show['compartment'] == compartment[0]][x_param],
                         data_to_show[data_to_show['compartment'] == compartment[0]][y_param])
        corr2 = pearsonr(data_to_show[data_to_show['compartment'] == compartment[0]][x_param],
                         data_to_show[data_to_show['compartment'] == compartment[0]][y_bg_param])
        corr3 = pearsonr(data_to_show[data_to_show['compartment'] == compartment[1]][x_param],
                         data_to_show[data_to_show['compartment'] == compartment[1]][y_param])
        corr4 = pearsonr(data_to_show[data_to_show['compartment'] == compartment[1]][x_param],
                         data_to_show[data_to_show['compartment'] == compartment[1]][y_bg_param])
        print("Correlation between {} , {} in {} {} = ".format(x_param, y_param, compartment[0], localization), corr1)
        print("Correlation between {} , {} in {} {} = ".format(x_param, y_bg_param, compartment[0], localization),
              corr2)
        print("Correlation between {} , {} in {} {} = ".format(x_param, y_param, compartment[1], localization), corr3)
        print("Correlation between {} , {} in {} {} = ".format(x_param, y_bg_param, compartment[1], localization),
              corr4)
        if save_it == 1:
            self.SaveFigures(op_file)
        plt.show()
    def Histogramplotting(self,df,num_bins,stats,hist_stat,alphas,colors,xlab,ylab,op_file= "",titles = [],legends =True,n_rows=1,n_cols=2,save_it = 0):
        fig,ax = plt.subplots(figsize=(6*n_cols,8*n_rows),nrows=n_rows,ncols=n_cols )
        try:
            ax = ax.flatten()
        except:
            ax = [ax]
        for i in range(n_cols):
            for j in range(len(stats)):
                sns.histplot(df[stats[j]], alpha=alphas[j],bins=num_bins,stat=hist_stat, legend=legends, color= colors[i], ax=ax[i], edgecolor='w').\
                    set( xlabel=xlab,ylabel=ylab)

            if not titles == []:
                ax[i].set_title(titles[i])
        fig.tight_layout()
        if save_it == 1:
            self.SaveFigures(op_file)
        plt.show()
