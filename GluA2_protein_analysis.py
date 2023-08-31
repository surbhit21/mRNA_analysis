import csv
# from functools import reduce
# import json
# from lmfit import conf_interval, minimize,Minimizer, Parameters, Parameter, report_fit, printfuncs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
# from matplotlib import rc
# from multiprocessing import Pool
import numpy as np
import os

# import PlotBinned
# from PlottingWidgetAMPA import *
import pandas as pd
# from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
#
# import scikit_posthocs as sp
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, kruskal,spearmanr,pearsonr
# from scipy import signal
from SNSPlottingWidget import SNSPlottingWidget
import seaborn as sns
from statannotations.Annotator import Annotator
# from Utility import *


scale_GluA2_protein = 0.24034
Lengths = np.array([25,50,75,100,150,250])
bin_size = 3
bins = np.arange(0, Lengths.max(), bin_size)
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
COLORS_dict = {"spine":"#005f73","shaft":'#CA6702',"spine_s":"#005f73","spine_i":'#CA6702',"shaft_s":"#005f73","shaft_i":'#CA6702'}
class GluA2DataAnalysis():
    def __init__(self, DataFolder):
        """
            class that loads data from folder
            DataFolder: Root folder that should have cell folders named as Neuron1, Neuron2 ...
        """
        self.df = DataFolder
        self.molecules = ["GFP", "surf-GluA2", "int-GluA2"];  # order of channels
        self.GFP_channel = 0
        self.int_channel = 2
        self.surf_channel = 1
        self.soma_stats = ["area", "mean", "stddev", "min", "max", "intden", "median", "rawintden"]
        self.dend_stats = ["area", "mean", "stddev", "min", "max", "intden", "median", "rawintden", "length"]
        self.conditions = ["", "-bg"]

    def LoadData(self, bins, num_cells, exclude_cells=[]):
        """
        Assumes a folder strucutre. Give the outermost folder which contains folder for each image
       df ==> cell_i ==> Rois ==> int-GLuA2/suef-GluA2/GFP
        """
        # fig,(ax0,ax1) = plt.subplots(1, 2,figsize=(20, 12), sharey=True)
        files = os.listdir(self.df)
        #         set of dictonaries to hold data
        int_glua2_data = {}
        surf_glua2_data = {}
        int_glua2_density = {}
        surf_glua2_density = {}
        ratio_int_surf = {}
        total_ratio_int_surf = {}
        GFP_data = {}
        soma_int_data = {}
        soma_surf_data = {}
        soma_gfp_data = {}
        raw_int_data = {}
        raw_surf_data = {}
        int_glua2_bg = {}
        surf_glua2_bg = {}

        for l in Lengths:
            surf_glua2_data[l] = []
            int_glua2_data[l] = []
            surf_glua2_density[l] = []
            int_glua2_density[l] = []
            ratio_int_surf[l] = []
            total_ratio_int_surf[l] = []
            GFP_data[l] = []
            raw_int_data[l] = []
            raw_surf_data[l] = []
            int_glua2_bg[l] = []
            surf_glua2_bg[l] = []
        for fdx in range(1, num_cells + 1):
            # print(file)
            file = "Neuron{}".format(fdx)
            if os.path.isdir(os.path.join(self.df, file)) and fdx not in exclude_cells:
                sub_folder = "/Rois/SUM/"
                actual_F = "/data"
                background_F = "/background"
                soma_measures = "/measures/all_soma_data.csv"
                dend_measures = "/measures/all_dend_data.csv"
                all_measures = "/measures/all_data.csv"
                file1 = file + sub_folder

                int_fname = os.path.join(self.df, file1 + self.molecules[self.int_channel] + actual_F);
                surf_fname = os.path.join(self.df, file1 + self.molecules[self.surf_channel] + actual_F)
                GFP_fname = os.path.join(self.df, file1 + self.molecules[self.GFP_channel] + actual_F)

                int_glua2 = os.listdir(int_fname)
                surf_glua2 = os.listdir(surf_fname)
                gfp_files = os.listdir(GFP_fname)
                # print(int_glua2==surf_glua2,surf_glua2==gfp_files)

                int_soma_file = os.path.join(self.df, file1 + self.molecules[self.int_channel] + all_measures);
                surf_soma_file = os.path.join(self.df, file1 + self.molecules[self.surf_channel] + all_measures);
                GFP_soma_file = os.path.join(self.df, file1 + self.molecules[self.GFP_channel] + all_measures);

                int_bg_fname = os.path.join(self.df, file1 + self.molecules[self.int_channel] + background_F);
                surf_bg_fname = os.path.join(self.df, file1 + self.molecules[self.surf_channel] + background_F)
                gfp_bg_fname = os.path.join(self.df, file1 + self.molecules[self.GFP_channel] + background_F)
                k = 0
                for (idx, sdx) in zip(int_glua2, surf_glua2):
                    # print(idx)
                    #  reading the data from csv files and removing the 1st, title line
                    int_data = self.ReadCSVFull(int_fname + "/" + idx)
                    surf_data = self.ReadCSVFull(surf_fname + "/" + sdx)
                    gfp_data = self.ReadCSVFull(GFP_fname + "/" + idx)
                    # breakpoint()

                    #                     reading the background data
                    int_bg_data = self.ReadCSVFull(int_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")
                    surf_bg_data = self.ReadCSVFull(surf_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")
                    gfp_bg_data = self.ReadCSVFull(gfp_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")

                    #                   background reduced data
                    """
                        The mean of the background roi is reduced from the roi data
                    """
                    corrected_int_data = int_data
                    corrected_int_data[:, 1] -= int_bg_data[:, 1].mean()  # - int_bg_data[:,-1]
                    corrected_surf_data = surf_data
                    corrected_surf_data[:, 1] -= surf_bg_data[:, 1].mean()  # - surf_bg_data[:,-1]
                    corrected_gfp_data = gfp_data
                    corrected_gfp_data[:, 1] -= gfp_bg_data[:, 1].mean()  # - gfp_bg_data[:,-1]

                    # binning the data for all three channels
                    corrected_binned_surf_data = self.BinnedSum(corrected_surf_data, bins, 0, idx)
                    corrected_binned_int_data = self.BinnedSum(corrected_int_data, bins, 0, idx)
                    corrected_binned_gfp_data = self.BinnedSum(corrected_gfp_data, bins, 0, idx)
                    # print(corrected_int_data,corrected_binned_surf_data)

                    #                   GFP normalized distribution
                    GFP_normed_surf_data = corrected_binned_surf_data[:, 1] / np.sqrt(corrected_binned_gfp_data[:, 1])
                    GFP_normed_int_data = corrected_binned_int_data[:, 1] / corrected_binned_gfp_data[:, 1]

                    # calculating the ratio for dendrites and for soma
                    binned_bg_surf_data = self.BinnedSum(surf_bg_data, bins, 0, idx)
                    binned_bg_int_data = self.BinnedSum(int_bg_data, bins, 0, idx)
                    binned_sum_ratio = corrected_binned_surf_data[:, 1] / corrected_binned_int_data[:, 1]
                    # print(binned_sum_ratio)
                    total_sum_ratio = corrected_binned_surf_data[:, 1].sum() / corrected_binned_int_data[:, 1].sum()
                    # print(total_sum_ratio)
                    # if (total_sum_ratio > 10).any():
                    #     print(int_fname)
                    # print(int_fname,idx,np.nanmax(binned_sum_ratio))
                    # calculating the binned sum

                    # binned_sum_int = self.BinnedSum(int_data,bins,0,idx)
                    # binned_sum_surf = self.BinnedSum(surf_data,bins,0,sdx)
                    # binned_gfp_data = self.BinnedSum(gfp_data,bins,0,sdx)
                    for l in Lengths:
                        if l <= int_data[-1, 0]:
                            # if(l==100):
                            #     print(int_data.shape,idx,fdx)
                            # appending the binned ratios for appropriate lengths
                            surf_glua2_data[l].append(corrected_binned_surf_data[:, 1])
                            int_glua2_data[l].append(corrected_binned_int_data[:, 1])
                            surf_glua2_density[l].append(GFP_normed_surf_data)
                            int_glua2_density[l].append(GFP_normed_int_data)
                            GFP_data[l].append(corrected_binned_gfp_data[:, 1])
                            # appending the ratios for appropriate lengths
                            ratio_int_surf[l].append(binned_sum_ratio)
                            total_ratio_int_surf[l].append(total_sum_ratio)
                            int_glua2_bg[l].append(binned_bg_surf_data)
                            surf_glua2_bg[l].append(binned_bg_int_data)

                            # getting raw values only upto length l

                            x_n = int(np.ceil(l / scale_GluA2_protein))
                            raw_int_data[l].append(int_data[:, 1])
                            raw_surf_data[l].append(surf_data[:, 1])

                self.ReadCSVFull(int_soma_file, soma_int_data, self.dend_stats, fdx, 1)
                self.ReadCSVFull(surf_soma_file, soma_surf_data, self.dend_stats, fdx, 1)
                self.ReadCSVFull(GFP_soma_file, soma_gfp_data, self.dend_stats, fdx, 1)

        for l in Lengths:
            surf_glua2_data[l] = np.asarray(surf_glua2_data[l])
            int_glua2_data[l] = np.asarray(int_glua2_data[l])
            surf_glua2_density[l] = np.asarray(surf_glua2_density[l])
            int_glua2_density[l] = np.asarray(int_glua2_density[l])
            ratio_int_surf[l] = np.asarray(ratio_int_surf[l])
            total_ratio_int_surf[l] = np.asarray(total_ratio_int_surf[l])
            raw_int_data[l] = np.asarray(raw_int_data[l])
            raw_surf_data[l] = np.asarray(raw_surf_data[l])
            GFP_data[l] = np.asanyarray(GFP_data[l])
            int_glua2_bg[l] = np.asanyarray(int_glua2_bg[l])
            surf_glua2_bg[l] = np.asanyarray(surf_glua2_bg[l])
        col_names = [item_two + item_one for item_one in self.conditions for item_two in self.dend_stats] + [
            "compartment"]
        int_df = pd.DataFrame.from_dict(soma_int_data,
                                        orient='index', columns=col_names)
        surf_df = pd.DataFrame.from_dict(soma_surf_data,
                                         orient='index', columns=col_names)
        gfp_df = pd.DataFrame.from_dict(soma_gfp_data,
                                        orient='index', columns=col_names)
        return int_glua2_data, surf_glua2_data, int_glua2_density, surf_glua2_density, \
               ratio_int_surf, total_ratio_int_surf, \
               int_df, surf_df, raw_int_data, raw_surf_data, GFP_data, soma_gfp_data, int_glua2_bg, surf_glua2_bg

    def BinnedSum(self, arr, bins, num_col=-1, name=None):
        # print(name)
        if len(arr.shape) == 2:
            rr, cc = arr.shape
            binned_sum = np.zeros((len(bins), cc))
            digitized = bins.searchsorted(arr[:, num_col])
            #             print(digitized,bins,arr[:,0])
            # breakpoint()
            digitized[0] = digitized[1]
            for c in range(0, cc):
                binned_sum[:, c] = np.bincount(digitized, weights=arr[:, c], minlength=len(bins))
            binned_sum[:, num_col] = bins
            return binned_sum[1:]
        else:
            print("quite not the shape", arr.shape)
            return np.zeros((len(bins), arr.shape[1]))

    def DendWiseDict(self, num_samples, data, data_dict, cols, cell_id):
        # rr,cc = data.shape
        # num_samples = rr//4
        # print(data)
        # data_dict = {}
        # print(data_dict)
        keys = data.keys()  # list(data_dict.keys())
        original_keys = [key for key in keys if "-bg" not in key]
        original_keys.remove('Label')
        for dict_key in original_keys:
            # print(dict_key)
            data_dict["Neuron_{}_{}".format(cell_id, dict_key)] = [float(ida) for ida in data[dict_key]] + [float(idb)
                                                                                                            for idb in
                                                                                                            data[
                                                                                                                dict_key + '-bg']]
            if dict_key.lower() == 'soma':
                data_dict["Neuron_{}_{}".format(cell_id, dict_key)].append("Soma")
            else:
                data_dict["Neuron_{}_{}".format(cell_id, dict_key)].append("Dendrite")
        return data_dict

    def ReadCSVFull(self, filename, data_dict=None, cols=None, cell_id=None, p_or_m=0):
        """
        function to read CSV data files:
            argumentrs : filename : path of the .csv file
                data_dict: data dict to store the read data, used only in case p_or_m == 1
                cols: number of columns in the file, used only in case p_or_m == 1
                cell_id: id of the cell
                p_or_m: =1 if it is a imageJ measure file =0 if it is a imagJ plot profile file

        """

        with open(filename) as csvDataFile:
            # read file as csv file
            csvReader = csv.reader(csvDataFile)
            # loop over rows
            num_row = -1
            if p_or_m == 1:
                csv_data = {}
                for row in csvReader:
                    # add cell [0] to list of dates
                    # print(row[1].split(':')[-1])
                    num_row += 1
                    key = row[1].split(':')[-1]
                    csv_data[key] = row[-1 * len(cols):]
                # print(csv_data)
                return self.DendWiseDict(num_row // 2, csv_data, data_dict, cols, cell_id)
            elif p_or_m == 0:
                csv_data = []
                for row in csvReader:
                    csv_data.append(row)
                data = np.asarray(csv_data[1:]).astype("float")
                # print(data[0,-2:])
                # print(data[:,1:])
                return data[:, 1:]
            else:
                raise ("some other type file sent, not proccessible")

    def GetSumNormDistribution(self, data):
        sums = data.sum(axis=1)
        norm_data = np.transpose(data) / sums
        return np.transpose(norm_data)

    def GetSomaNormDistribution(self, data, index=1):
        d_vec = data[:, index]
        norm_data = data / d_vec[:, None]
        return norm_data

    def GetNormalizedDendriticDistribution(self, int_data, surf_data, soma_norm=1):

        norm_int = {}
        norm_surf = {}
        mean_int = {}
        mean_surf = {}
        std_int = {}
        std_surf = {}
        sem_int = {}
        sem_surf = {}
        ratio_mean = {}
        ratio_std = {}
        ratio_sem = {}
        off_set = 0
        for l in Lengths[:-1]:
            x = np.arange(0, l, bin_size)
            if soma_norm == 1:
                norm_int[l] = self.GetSomaNormDistribution(int_data[l], index=off_set)
                norm_surf[l] = self.GetSomaNormDistribution(surf_data[l], index=off_set)
            else:
                norm_int[l] = int_data[l]
                norm_surf[l] = surf_data[l]
            print(int_data[l].shape)

            mean_int[l] = norm_int[l].mean(axis=0)[off_set:x.shape[0]+off_set]
            std_int[l] = norm_int[l].std(axis=0)[off_set:x.shape[0]+off_set]  # (int_std[1]/int_mean[1] + int_std/int_mean)*mean_int[l]
            sem_int[l] = std_int[l] / np.sqrt(norm_int[l].shape[0])



            mean_surf[l] = norm_surf[l].mean(axis=0)[off_set:x.shape[0]+off_set]
            std_surf[l] = norm_surf[l].std(axis=0)[off_set:x.shape[0]+off_set]  # (surf_std[1]/surf_mean[1] + surf_std/surf_mean)*mean_surf[l]
            sem_surf[l] = std_surf[l] / np.sqrt(norm_surf[l].shape[0])
            # print((surf_std[1]/surf_mean[1] + surf_std/surf_mean)*mean_surf[l],norm_surf[l].std(axis=0))
            # ratio_mean[l] = np.nanmean(ratio_int_surf[l],axis=0)
            # ratio_std[l] = np.nanstd(ratio_int_surf[l],axis=0)
            # ratio_sem[l] = ratio_std[l]/np.sqrt(ratio_std[l].shape[0])


        return mean_surf, mean_int, sem_surf, sem_int, std_surf, std_int,  norm_int, norm_surf


G2DA_c = GluA2DataAnalysis("/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/Control/whole_neuron")

int_glua2_data, surf_glua2_data,int_glua2_density, surf_glua2_density,ratio_int_surf,\
    total_ratio_int_surf,int_df,\
        surf_df, raw_int_data,raw_surf_data,\
GFP_data,GFP_soma,int_glua2_bg,surf_glua2_bg = G2DA_c.LoadData(bins,15,exclude_cells=[])
mean_surf,mean_int,sem_surf,sem_int,std_surf,std_int,norm_int_c,norm_surf_c = G2DA_c.GetNormalizedDendriticDistribution(int_glua2_data,surf_glua2_data)
mean_surf_density,mean_int_density,sem_surf_density,sem_int_density,std_surf_density,std_int_density,norm_int_c_density,norm_surf_c_density = G2DA_c.GetNormalizedDendriticDistribution(int_glua2_density,surf_glua2_density)
mean_GFP,mean_GFP,sem_GFP,sem_GFP,std_GFP,std_GFP,norm_GFP,norm_GFP = G2DA_c.GetNormalizedDendriticDistribution(GFP_data,GFP_data)
# breakpoint()
plt_widget = SNSPlottingWidget()
roi_area_ratio = pd.DataFrame()
roi_area_ratio['area'] = surf_df['area']/surf_df['area-bg']
dendrites_surf = surf_df[surf_df.compartment == 'Dendrite']
somas_surf = surf_df[surf_df.compartment == 'Soma']
dendrites_int = int_df[surf_df.compartment == 'Dendrite']
somas_int = int_df[surf_df.compartment == 'Soma']
roi_area_ratio['length'] = dendrites_int["length"]/dendrites_int["length-bg"]
fig,ax = plt.subplots(figsize=(12,8),ncols=1,nrows=1)
# ratio_of_ratio = ratio_s_d_ind.eval("int / surf").rename("Ratio of Int/surf ratio")
x1 = np.ones(roi_area_ratio['length'].size)
sns.scatterplot(data=roi_area_ratio,ax=ax).set(xticklabels=[])
y1 = x1
sns.lineplot(x1,y1,color='r',label="y=1")
plt.show()

## histogram of foreground and background Fluorescent

num_bins = 10
stats = np.array([["mean","mean-bg"]])
alphas = np.array([0.8,0.2])
colors = np.array([COLORS_dict["spine_i"],COLORS_dict["shaft_i"]])
op_folder = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Max-Kracht/GluA2/Figures/"
plt_widget.CreateFolderRecursive(op_folder)
x_lab = "Mean fluorescent intesity [a.b.u]"
y_lab = "Count"
save_it = 1
# plt_widget.Histogramplotting(dendrites_int,num_bins,stats,"count",alphas,colors,x_lab,y_lab,titles = ["Dendrite Internal"], \
# op_file=os.path.join(op_folder,"Mean_intensity_dend_internal_glua2_histogram"),n_cols=1,save_it=save_it)
# plt_widget.Histogramplotting(dendrites_surf,num_bins,stats,"count",alphas,colors,x_lab,y_lab,titles = ["Dendrite Surface"], \
# op_file=os.path.join(op_folder,"Mean_intensity_dend_surface_glua2_histogram"),n_cols=1,save_it=save_it)
# plt_widget.Histogramplotting(somas_int,num_bins,stats,"count",alphas,colors,x_lab,y_lab,titles = ["Soma Internal"], \
# op_file=os.path.join(op_folder,"Mean_intensity_soma_internal_glua2_histogram"),n_cols=1,save_it=save_it)
# plt_widget.Histogramplotting(somas_surf,num_bins,stats,"count",alphas,colors,x_lab,y_lab,titles = ["Soma Surface"], \
# op_file=os.path.join(op_folder,"Mean_intensity_soma_surface_glua2_histogram"),n_cols=1,save_it=save_it)

## Coorelation plot between ROI mean intensity and area
stat1 = "mean"
stat2 = "area"
compartments = ["Soma","Dendrite"]
y_param = stat1
x_param = stat2
y_bg_param = stat1+"-bg"
# dts = concat_data
# plt_widget.CorrelationCalculationAndPlotting(surf_df,"surface",compartments,x_param,y_param,y_bg_param,\
# op_file=os.path.join(op_folder,"Coor_plot_bw_{}_{}_Surface_GluA2".format(stat1,stat2)),f_scale=1,save_it=save_it)
# plt_widget.CorrelationCalculationAndPlotting(int_df,"internal",compartments,x_param,y_param,y_bg_param,\
# op_file=os.path.join(op_folder,"Coor_plot_bw_{}_{}_internal_GluA2".format(stat1,stat2)),f_scale=1,save_it=save_it)


"""
### Surface to internal ratio in Soma vs Dendritic Compartment

stat = "mean"
stat_bg = stat+"-bg"
ratio_s_i_ind = pd.DataFrame()
ratio_s_i_ind['ratio']  = (surf_df[stat]-surf_df[stat_bg])/(int_df[stat] - int_df[stat_bg])
ratio_s_i_ind["compartment"] = surf_df.compartment
ratios_soma = ratio_s_i_ind[ratio_s_i_ind['compartment'] == "Soma"]['ratio']
ratios_dendrite = ratio_s_i_ind[ratio_s_i_ind['compartment'] == "Dendrite"]['ratio']
fig,ax = plt.subplots(figsize=(8,6),nrows=1,ncols=1)
sns.set(font_scale = 1.5)
sns.boxplot(data=ratio_s_i_ind,x="compartment",y="ratio",notch=True,flierprops={"marker": "8","linestyle":''},ax=ax,palette = [COLORS_dict["shaft_s"],COLORS_dict["shaft_i"]]).set(
    xlabel="Compartment",ylabel=r"Fluorescent intesity ratio $\frac{surface}{internal}$",title="surface to internal {} intesity ratio".format(stat))
pairs = [("Soma","Dendrite")]
annotator = Annotator(ax, data=ratio_s_i_ind,x="compartment",y="ratio",pairs=pairs)
annotator.configure(test="Kruskal")
annotator.apply_and_annotate()
ax.set_ylim([0.2,0.6])
print("Soma: Mean {:.2f}, std ({:.3f})".format(ratios_soma.mean(),ratios_soma.std()))
print("Dendrite: Mean {:.2f}, std ({:.3f})".format(ratios_dendrite.mean(),ratios_dendrite.std()),ratios_dendrite.shape,ratios_soma.shape)
# print(ratio_s_i_ind.head(),ratio_s_i_ind.columns.values)
plt.show()

"""
## dendritic distribution of the ratio of surface to internal GluA2

def line_fit(x, m,c):
    return m*x + c

def exp_fit(x,a,b):
    return a*np.exp(-b*x)

def R_seq(y_fit,y_orig):
    ss_res = ((y_orig-y_fit)**2).sum()
    ss_tot = ((y_orig-y_orig.mean())**2).sum()
    # print("in R_Seq =",ss_tot,ss_res)
    return 1 - (ss_res/ss_tot)

def ChiSq(yd,y_fit,sigmas):
    nzs = np.nonzero(sigmas)
    # print(nzs)
    r_yd = np.take(yd,nzs)
    r_yf = np.take(y_fit,nzs)
    r_sgs = np.take(sigmas,nzs)
    residuals = r_yd - r_yf
    chi_squ = np.sum((residuals/ r_sgs)**2)
    return chi_squ
to_analyse = Lengths[:-2]
fig,ax1 = plt.subplots(figsize=(8,6*len(to_analyse)),ncols=1,nrows=len(to_analyse))

off_set = 0
for ldx,l1 in enumerate(to_analyse):
    # print(ldx,l1)
    x = np.arange(0, l1, bin_size)
    mean_ratio = ratio_int_surf[l1].mean(axis=0)[off_set:x.shape[0]+off_set]
    # print("For Length {}, params = ".format(l1),ratio_int_surf[l1].shape,mean_ratio)
    std_ratio = ratio_int_surf[l1].std(axis=0)[off_set:x.shape[0]+off_set]
    sem_ratio = (std_ratio/np.sqrt(ratio_int_surf[l1].shape[0]))
    popt, pcov = curve_fit(line_fit, x, mean_ratio)
    popt_e, pcov_e = curve_fit(exp_fit, x, mean_ratio)
    # print("For Length {}, params = ".format(l1),popt_e)
    # print("ratios = ",ratio_int_surf[l1][:,13])
    mean_fit = line_fit(x,*popt)
    mean_fit_exp = exp_fit(x,*popt_e)
    r2 = R_seq(mean_fit,mean_ratio)
    r2_exp = ChiSq(mean_ratio,mean_fit_exp,sem_ratio)
    # print(r2,r2_exp)
    # print("popt = ",popt,"pcov = ",pcov)
    print("popt_e = ",popt_e,"pcov_e = ",pcov_e)
     # ax1[ldx].plot(x,mean_fit,'r-.',label=r"line-fit,$R^2 = %.2f$" % r2)
    ax1[ldx].errorbar(x+bin_size/2,mean_ratio,sem_ratio,color='#cc8963ff',label='data',marker='o',linestyle="None")
    ax1[ldx].plot(x,mean_fit_exp,'k--',label=r"exp-fit,$\chi^2 = %.2f$" % r2_exp)
    # sns.regplot(x,mean_ratio,ax=ax1[ldx],line_kws={"alpha": 0.2})
    pea_ratio = spearmanr(x,mean_ratio)
    f_size = 12
    loc = "upper right"
    # ax1[ldx].text(.05, .8, 'r={:.2f}, p={:.2g}'.format(pea_ratio[0],pea_ratio[1]),transform=ax1[ldx].transAxes)
    ax1[ldx].set_xlabel(r"Dendritic distance (in $\mu m$), N ={}".format( ratio_int_surf[l1].shape[0]),size=f_size)
    ax1[ldx].set_ylabel(r"Fluorescent intensity ratio $\frac{surface}{internal}$",size=f_size)
    ax1[ldx].set_title("Dendritic distribution of surface to internal fluorescent intesity ratio",size=f_size)
    ax1[ldx].legend(fontsize=f_size,loc=loc)
    ax1[ldx].set_xlim([-1,l1+1])

    # ax1[ldx].set_ylim([0.3,0.6])
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Surface_to_internal_glua2_ratio_exp_fitted"))
plt.show()

l1 = 100
x_range = np.arange(0,l1,bin_size)
ratio_mean = ratio_int_surf[l1].mean(axis=0)[off_set:x.shape[0]+off_set]
ratio_std = ratio_int_surf[l1].std(axis=0)[off_set:x.shape[0]+off_set]
ratio_sem = ratio_std/np.sqrt(ratio_int_surf[l1].shape[0])
plt_widget.PlotBinnedStats(np.asarray([x_range]), np.asarray([ratio_mean/ratio_mean[0]]), np.asarray([ratio_sem]),
                           np.asarray([ratio_mean]), ["cnih2"], x_lab, y_lab, ["Norm CNIH2 density"], ["#DDA15E"],"",
                           op_folder+"_norm__with_CNIH2_ss_mode_fit",
                                    bin_size,save_it = 0,fit_exp=3,in_set=0,set_axis_label=0)


fig,ax2 = plt.subplots(figsize=(8,6*len(to_analyse)),ncols=1,nrows=len(to_analyse))
# breakpoint()
off_set = 0
for ldx,l1 in enumerate(to_analyse):
    x=  np.arange(0, l1, bin_size)
    ax2[ldx].errorbar(x+bin_size/2,mean_surf[l1],sem_surf[l1],color=COLORS_dict["shaft_s"],label="surface",marker='s',linestyle="-")
    ax2[ldx].errorbar(x+bin_size/2,mean_int[l1],sem_int[l1],color=COLORS_dict["shaft_i"],label="internal",marker='o',linestyle="--")
    # ax2[ldx].errorbar(x+bin_size/2,mean_GFP[l1],sem_GFP[l1],label = "GFP",color="green",marker='8',linestyle="-." )
    ax2[ldx].set_ylabel("Normalized fluorescence [a.b.u]",fontsize=14)

    ax2[ldx].set_xlabel(r"Dendritic distance in $\mu m$, N= {}".format(int_glua2_data[l1].shape[0]),fontsize=14)
    ax2[ldx].set_title("Dendritic distribution of normalized fluorescence intensity",fontsize=14)
    ax2[ldx].legend()
    # ax2[1].legend()
    ax2[ldx].set_xlim([-1,l1+2])

    # ax2[1].set_xlim([-4,100])
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Glua2_normalized_fluorescent_intensity"))
plt.show()


"""
Fitting the AMPA GLuA2 model to the normalized fluorescent intensity 
"""

length_to_fit = 100 #fitting the 100 micron long data