import csv
import json
import os
import pandas as pd
import pickle
from protein_model_fitting import *
import Total_AMPAR
import seaborn as sns
from SNSPlottingWidget import SNSPlottingWidget
from Utility import *

scale_CNIH2_protein = 0.12
COLORS = ["#005f73","#9b2226","#CA6702","#337357"]
COLORS_dict = {"spine":"#005f73","shaft":'#CA6702',"spine_s":"#005f73","spine_i":'#CA6702',"shaft_s":"#005f73","shaft_i":'#CA6702'}


# A class that reads, bins, normalizes data
class CNHI2DataAnalysis():
    def __init__(self, DataFolder):
        """
            class that loads data from folder
            DataFolder: Root folder that should have cell folders named as Neuron1, Neuron2 ...
        """
        self.df = DataFolder
        self.molecules = ["DAPI", "MAP2","CNIH2", "GluA2"];  # order of channels
        self.DAPI_channel = 0
        self.MAP2_channel = 1
        self.CNIH2_channel = 2
        self.Glua2_channel = 3
        self.soma_stats = ["area", "mean", "stddev", "min", "max", "intden", "median", "rawintden"]
        self.dend_stats = ["area", "mean", "stddev", "min", "max", "intden", "median", "rawintden", "length"]
        self.conditions = ["", "-bg"]

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

    def GetNormalizedDendriticDistribution(self, Glua2_data, CNIH2_data, soma_norm=1):

        norm_CNIH2 = {}
        norm_GluA2 = {}
        mean_CNIH2 = {}
        mean_GluA2 = {}
        std_CNIH2 = {}
        std_GluA2 = {}
        sem_CNIH2 = {}
        sem_GluA2 = {}
        ratio_mean = {}
        ratio_std = {}
        ratio_sem = {}
        off_set = 0
        for l in Lengths[:-2]:
            x = np.arange(0, l, bin_size)
            x1 = np.arange(0, l, scale_CNIH2_protein)
            if l in CNIH2_data.keys():
                if soma_norm == 1 :
                    norm_CNIH2[l] = self.GetSomaNormDistribution(CNIH2_data[l], index=off_set)
                    norm_GluA2[l] = self.GetSomaNormDistribution(Glua2_data[l], index=off_set)
                else:
                    norm_CNIH2[l] = CNIH2_data[l]
                    norm_GluA2[l] = Glua2_data[l]
                print(CNIH2_data[l].shape)
                CNIH2_mean = CNIH2_data[l].mean(axis=0)
                CNIH2_std = CNIH2_data[l].std(axis=0)

                mean_CNIH2[l] = norm_CNIH2[l].mean(axis=0)[off_set:x.shape[0] + off_set]
                std_CNIH2[l] = norm_CNIH2[l].std(axis=0)[off_set:x.shape[
                                                                     0] + off_set]  # (CNIH2_std[1]/CNIH2_mean[1] + CNIH2_std/CNIH2_mean)*mean_CNIH2[l]
                sem_CNIH2[l] = std_CNIH2[l] / np.sqrt(norm_CNIH2[l].shape[0])

                Glua2_mean = Glua2_data[l].mean(axis=0)
                Glua2_std = Glua2_data[l].std(axis=0)

                mean_GluA2[l] = norm_GluA2[l].mean(axis=0)[off_set:x.shape[0] + off_set]
                std_GluA2[l] = norm_GluA2[l].std(axis=0)[off_set:x.shape[
                                                                     0] + off_set]  # (Glua2_std[1]/Glua2_mean[1] + Glua2_std/Glua2_mean)*mean_GluA2[l]
                sem_GluA2[l] = std_GluA2[l] / np.sqrt(norm_GluA2[l].shape[0])

                norm_CNIH2[l] = norm_CNIH2[l][off_set:x.shape[0] + off_set]
                norm_GluA2[l] = norm_GluA2[l][off_set:x.shape[0] + off_set]

        return mean_GluA2, sem_GluA2, std_GluA2, norm_GluA2, \
               mean_CNIH2, sem_CNIH2, std_CNIH2, norm_CNIH2

    def LoadData(self, bins, num_cells, exclude_cells=[]):
        """
        Assumes a folder strucutre. Give the outermost folder which contains folder for each image
       df ==> cell_i ==> Rois ==> int-GLuA2/suef-GluA2/MAP2
        """
        # fig,(ax0,ax1) = plt.subplots(1, 2,figsize=(20, 12), sharey=True)
        files = os.listdir(self.df)
        #         set of dictonaries to hold data
        CNIH2_length_data = {}
        Glua2_length_data = {}
        CNIH2_length_density = {}
        Glua2_length_density = {}
        ratio_CNIH2_GluA2 = {}
        total_ratio_CNIH2_GluA2 = {}
        MAP2_length_data = {}
        soma_CNIH2_data = {}
        soma_Glua2_data = {}
        soma_MAP2_data = {}
        raw_CNIH2_data = {}
        raw_Glua2_data = {}
        raw_MAP2_data = {}
        CNIH2_bg = {}
        Glua2_bg = {}

        for l in Lengths:
            Glua2_length_data[l] = []
            CNIH2_length_data[l] = []
            Glua2_length_density[l] = []
            CNIH2_length_density[l] = []
            ratio_CNIH2_GluA2[l] = []
            total_ratio_CNIH2_GluA2[l] = []
            MAP2_length_data[l] = []
            raw_CNIH2_data[l] = []
            raw_Glua2_data[l] = []
            raw_MAP2_data[l] = []
            CNIH2_bg[l] = []
            Glua2_bg[l] = []
        for fdx in range(1, num_cells + 1):
            file = "cell_{}".format(fdx)
            # print(file)
            if os.path.isdir(os.path.join(self.df, file)) and fdx not in exclude_cells:
                sub_folder = "/Rois/SUM/"
                actual_F = "/data"
                background_F = "/background"
                soma_measures = "/measures/all_soma_data.csv"
                dend_measures = "/measures/all_dend_data.csv"
                all_measures = "/measures/all_data.csv"
                file1 = file + sub_folder

                CNIH2_fname = os.path.join(self.df, file1 + self.molecules[self.CNIH2_channel] + actual_F);
                Glua2_fname = os.path.join(self.df, file1 + self.molecules[self.Glua2_channel] + actual_F)
                MAP2_fname = os.path.join(self.df, file1 + self.molecules[self.MAP2_channel] + actual_F)

                CNIH2_files = os.listdir(CNIH2_fname)
                Glua2_files = os.listdir(Glua2_fname)
                MAP2_files = os.listdir(MAP2_fname)
                # print(CNIH2_files,Glua2_files,MAP2_files)

                CNIH2_soma_file = os.path.join(self.df, file1 + self.molecules[self.CNIH2_channel] + all_measures);
                Glua2_soma_file = os.path.join(self.df, file1 + self.molecules[self.Glua2_channel] + all_measures);
                MAP2_soma_file = os.path.join(self.df, file1 + self.molecules[self.MAP2_channel] + all_measures);

                CNIH2_bg_fname = os.path.join(self.df, file1 + self.molecules[self.CNIH2_channel] + background_F);
                Glua2_bg_fname = os.path.join(self.df, file1 + self.molecules[self.Glua2_channel] + background_F)
                MAP2_bg_fname = os.path.join(self.df, file1 + self.molecules[self.MAP2_channel] + background_F)
                k = 0
                for (idx, sdx) in zip(CNIH2_files, Glua2_files):
                    # print(idx)
                    #  reading the data from csv files and removing the 1st, title line
                    CNIH2_data = self.ReadCSVFull(CNIH2_fname + "/" + idx)
                    Glua2_data = self.ReadCSVFull(Glua2_fname + "/" + idx)
                    MAP2_data = self.ReadCSVFull(MAP2_fname + "/" + idx)
                    # breakpoint()

                    #                    reading the background data
                    CNIH2_bg_data = self.ReadCSVFull(CNIH2_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")
                    Glua2_bg_data = self.ReadCSVFull(Glua2_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")
                    MAP2_bg_data = self.ReadCSVFull(MAP2_bg_fname + "/" + idx.split(".")[0] + "-bg.csv")

                    #                   background reduced data
                    """
                        The mean of the background roi is reduced from the roi data
                    """
                    corrected_CNIH2_data = CNIH2_data
                    # corrected_CNIH2_data[:,1] -= CNIH2_bg_data[:,1].mean() #- CNIH2_bg_data[:,-1]
                    corrected_Glua2_data = Glua2_data
                    # corrected_Glua2_data[:,1] -= Glua2_bg_data[:,1].mean()# - Glua2_bg_data[:,-1]
                    corrected_MAP2_data = MAP2_data
                    # corrected_MAP2_data[:,1] -= MAP2_bg_data[:,1].mean() #- MAP2_bg_data[:,-1]
                    if (np.max(Glua2_data[:,1]) > 975 ):
                        print("cell = ",fdx,idx)
                    # binning the data for all three channels
                    corrected_binned_Glua2_data = self.BinnedSum(corrected_Glua2_data, bins, 0, idx)
                    corrected_binned_CNIH2_data = self.BinnedSum(corrected_CNIH2_data, bins, 0, idx)
                    corrected_binned_MAP2_data = self.BinnedSum(corrected_MAP2_data, bins, 0, idx)
                    # print(corrected_CNIH2_data,corrected_binned_Glua2_data)

                    #                   MAP2 normalized distribution
                    MAP2_normed_Glua2_data = corrected_binned_Glua2_data[:, 1] / corrected_binned_MAP2_data[:, 1]
                    MAP2_normed_CNIH2_data = corrected_binned_CNIH2_data[:, 1] / corrected_binned_MAP2_data[:, 1]

                    # calculating the ratio for dendrites and for soma
                    binned_bg_Glua2_data = self.BinnedSum(Glua2_bg_data, bins, 0, idx)
                    binned_bg_CNIH2_data = self.BinnedSum(CNIH2_bg_data, bins, 0, idx)
                    binned_sum_ratio = corrected_binned_Glua2_data[:, 1] / corrected_binned_CNIH2_data[:, 1]
                    # print(binned_sum_ratio)
                    total_sum_ratio = corrected_binned_Glua2_data[:, 1].sum() / corrected_binned_CNIH2_data[:, 1].sum()
                    # print(total_sum_ratio)

                    # calculating the binned sum

                    # binned_sum_CNIH2 = self.BinnedSum(CNIH2_data,bins,0,idx)
                    # binned_sum_GluA2 = self.BinnedSum(Glua2_data,bins,0,sdx)
                    # binned_MAP2_data = self.BinnedSum(MAP2_data,bins,0,sdx)
                    for l in Lengths:
                        if l <= CNIH2_data[-1, 0]:
                            # if(l==100):
                            #     print(CNIH2_data.shape,idx,fdx)
                            # appending the binned ratios for appropriate lengths
                            Glua2_length_data[l].append(corrected_binned_Glua2_data[:, 1])
                            CNIH2_length_data[l].append(corrected_binned_CNIH2_data[:, 1])
                            Glua2_length_density[l].append(MAP2_normed_Glua2_data)
                            CNIH2_length_density[l].append(MAP2_normed_CNIH2_data)
                            MAP2_length_data[l].append(corrected_binned_MAP2_data[:, 1])
                            # appending the ratios for appropriate lengths
                            ratio_CNIH2_GluA2[l].append(binned_sum_ratio)
                            total_ratio_CNIH2_GluA2[l].append(total_sum_ratio)
                            CNIH2_bg[l].append(binned_bg_Glua2_data)
                            Glua2_bg[l].append(binned_bg_CNIH2_data)

                            # getting raw values only upto length l

                            x_n = int(np.ceil(l / scale_CNIH2_protein))
                            raw_CNIH2_data[l].append(CNIH2_data[0:x_n, 1])
                            raw_Glua2_data[l].append(Glua2_data[0:x_n, 1])
                            raw_MAP2_data[l].append(MAP2_data[0:x_n,1])
                self.ReadCSVFull(CNIH2_soma_file, soma_CNIH2_data, self.dend_stats, fdx, 1)
                self.ReadCSVFull(Glua2_soma_file, soma_Glua2_data, self.dend_stats, fdx, 1)
                self.ReadCSVFull(MAP2_soma_file, soma_MAP2_data, self.dend_stats, fdx, 1)

        for l in Lengths:
            Glua2_length_data[l] = np.asarray(Glua2_length_data[l])
            CNIH2_length_data[l] = np.asarray(CNIH2_length_data[l])
            Glua2_length_density[l] = np.asarray(Glua2_length_density[l])
            CNIH2_length_density[l] = np.asarray(CNIH2_length_density[l])
            ratio_CNIH2_GluA2[l] = np.asarray(ratio_CNIH2_GluA2[l])
            total_ratio_CNIH2_GluA2[l] = np.asarray(total_ratio_CNIH2_GluA2[l])
            raw_CNIH2_data[l] = np.asarray(raw_CNIH2_data[l])
            raw_Glua2_data[l] = np.asarray(raw_Glua2_data[l])
            raw_MAP2_data[l] = np.asarray(raw_MAP2_data[l])
            MAP2_length_data[l] = np.asanyarray(MAP2_length_data[l])
            CNIH2_bg[l] = np.asanyarray(CNIH2_bg[l])
            Glua2_bg[l] = np.asanyarray(Glua2_bg[l])
        col_names = [item_two + item_one for item_one in self.conditions for item_two in self.dend_stats] + [
            "compartment"]
        CNIH2_df = pd.DataFrame.from_dict(soma_CNIH2_data,
                                          orient='index', columns=col_names)
        Glua2_df = pd.DataFrame.from_dict(soma_Glua2_data,
                                          orient='index', columns=col_names)
        MAP2_df = pd.DataFrame.from_dict(soma_MAP2_data,
                                         orient='index', columns=col_names)
        # breakpoint()
        return CNIH2_length_data, Glua2_length_data, CNIH2_length_density, Glua2_length_density, \
               ratio_CNIH2_GluA2, total_ratio_CNIH2_GluA2, \
               CNIH2_df, Glua2_df, raw_CNIH2_data, raw_Glua2_data,raw_MAP2_data, MAP2_length_data, soma_MAP2_data, CNIH2_bg, Glua2_bg


def DumpDict(datadict,fname):
    with open(fname,'wb') as outfile:
        pickle.dump(datadict,outfile, protocol=pickle.HIGHEST_PROTOCOL)
    print("{} saved!".format(fname))

def ReadDataDict(fname):
    with open(fname,'rb') as infile:
        d = pickle.load(infile)
    print("{} loaded!".format(fname))
    return d

def ReadData():
    # dend_cell_sum = {}
    curr_wd = os.getcwd()
    print(curr_wd)
    return [
        ReadDataDict(os.path.join(curr_wd,"Protein_data/CNIH2_GluA2/CNIH2_length_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/Glua2_length_data.pickle")),
        ReadDataDict(os.path.join(curr_wd,"Protein_data/CNIH2_GluA2/CNIH2_length_density.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/Glua2_length_density.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/ratio_CNIH2_GluA2.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/total_ratio_CNIH2_GluA2.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/CNIH2_df.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/Glua2_df.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/raw_CNIH2_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/raw_GluA2_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/raw_MAP2_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/MAP2_length_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/soma_MAP2_data.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/CNIH2_bg.pickle")),
        ReadDataDict(os.path.join(curr_wd, "Protein_data/CNIH2_GluA2/Glua2_bg.pickle")),
    ]

G2DA_c = CNHI2DataAnalysis("")

CNIH2_length_data, \
Glua2_length_data,\
CNIH2_length_density, \
Glua2_length_density,\
ratio_CNIH2_GluA2,\
total_ratio_CNIH2_GluA2, \
CNIH2_df,\
Glua2_df,\
raw_CNIH2_data,\
raw_Glua2_data,\
raw_MAP2_data,\
MAP2_length_data,\
soma_MAP2_data,\
CNIH2_bg,\
Glua2_bg = ReadData()
#

# mean_Glua2,sem_GluA2,std_Glua2,norm_GluA2_c,\
# mean_CNIH2,sem_CNIH2,std_CNIH2,norm_CNIH2_c = G2DA_c.GetNormalizedDendriticDistribution(Glua2_length_data,CNIH2_length_data)
mean_Glua2_density,sem_Glua2_density,std_Glua2_density,norm_GluA2_c_density,\
mean_CNIH2_density,sem_CNIH2_density,std_CNIH2_density,norm_CNIH2_c_density= G2DA_c.GetNormalizedDendriticDistribution(Glua2_length_density,CNIH2_length_density)
# print(MAP2_length_data)
mean_MAP2_density,sem_MAP2_density,std_MAP2_density,norm_MAP2_c_density,\
mean_MAP2_density,sem_MAP2_density,std_MAP2_density,norm_MAP2_c_density= G2DA_c.GetNormalizedDendriticDistribution(MAP2_length_data,MAP2_length_data)

breakpoint()

# loading the total glua2 data from M. Kracht
total_glua2_mk = {}
with open("./total_glua2_data.json") as fp:
    total_glua2_mk = json.load(fp)
total_glua2_mk = {int(k):np.array(v) for k,v in total_glua2_mk.items()}
mean_total_Glua2_mk,sem_total_Glua2_mk,std_total_Glua2_mk,norm_total_Glua2_mk,\
mean_total_Glua2_mk,sem_total_Glua2_mk,std_total_Glua2_mk,norm_total_Glua2_mk = G2DA_c.GetNormalizedDendriticDistribution(total_glua2_mk,total_glua2_mk)
# breakpoint()

# breakpoint()
plt_widget = SNSPlottingWidget()
# ## sanity check to see if the Rois and Background has same area
#
# roi_area_ratio = pd.DataFrame()
# roi_area_ratio['area'] = Glua2_df['area']/Glua2_df['area-bg']
# dendrites_Glua2 = Glua2_df[Glua2_df.compartment == 'Dendrite']
# somas_Glua2 = Glua2_df[Glua2_df.compartment == 'Soma']
# dendrites_CNIH2 = CNIH2_df[Glua2_df.compartment == 'Dendrite']
# somas_CNIH2 = CNIH2_df[CNIH2_df.compartment == 'Soma']
# roi_area_ratio['length'] = dendrites_CNIH2["length"]/dendrites_CNIH2["length-bg"]
# fig,ax = plt.subplots(figsize=(12,8),ncols=1,nrows=1)
# # ratio_of_ratio = ratio_s_d_ind.eval("int / surf").rename("Ratio of Int/surf ratio")
# x1 = np.ones(roi_area_ratio['length'].size)
# sns.scatterplot(data=roi_area_ratio,ax=ax).set(xticklabels=[])
# # plt.rc('font', family='Arial')
# y1 = x1
# sns.lineplot(x1,y1,color='r',label="y=1")
# plt.show()




## histogram of foreground and background Fluorescent

num_bins = 10
stats = np.array([["mean","mean-bg"]])
alphas = np.array([0.8,0.2])
colors = np.array([COLORS_dict["spine_i"],COLORS_dict["shaft_i"]])

# plotting function signature
# plt_widget.Histogramplotting(df,num_bins,stats,hist_stat,alphas,colors,xlab,ylab,op_file= "",titles = [],legends =True,n_rows=1,n_cols=2,save_fig = 0)
op_folder = "./Figures/Protein/"
plt_widget.CreateFolderRecursive(op_folder)
x_lab = "Mean fluorescent intesity [a.b.u]"
y_lab = "Count"
save_it = 1
# plt_widget.Histogramplotting(dendrites_CNIH2,10,stats,"count",alphas,colors,x_lab,y_lab,op_file=os.path.join(op_folder,"Mean_intensity_dend_CNIH2_histogram"),n_cols=1,save_it=save_it)
# plt_widget.Histogramplotting(dendrites_Glua2,10,stats,"count",alphas,colors,x_lab,y_lab,op_file=os.path.join(op_folder,"Mean_intensity_dend_GluA2_histogram"),n_cols=1,save_it=save_it)
#
# plt_widget.Histogramplotting(somas_CNIH2,5,stats,"count",alphas,colors,x_lab,y_lab,op_file=os.path.join(op_folder,"Mean_intensity_soma_GluA2_histogram"),n_cols=1,save_it=save_it)
# plt_widget.Histogramplotting(somas_Glua2,5,stats,"count",alphas,colors,x_lab,y_lab,op_file=os.path.join(op_folder,"Mean_intensity_soma_GluA2_histogram"),n_cols=1,save_it=save_it)

stat1 = "mean"
stat2 = "area"
compartments = ["Soma","Dendrite"]
y_param = stat1
x_param = stat2
y_bg_param = stat1+"-bg"
# dts = concat_data
# plt_widget.CorrelationCalculationAndPlotting(Glua2_df,"GluA2",compartments,x_param,y_param,y_bg_param,op_file=os.path.join(op_folder,"Coor_plot_bw_{}_{}_GluA2".format(stat1,stat2)),f_scale_CNIH2_protein=1,save_it=save_it)
# plt_widget.CorrelationCalculationAndPlotting(CNIH2_df,"CNIH2",compartments,x_param,y_param,y_bg_param,op_file=os.path.join(op_folder,"Coor_plot_bw_{}_{}_CNIH2".format(stat1,stat2)),f_scale_CNIH2_protein=1,save_it=save_it)

fsize = plt_widget.fsize
to_analyse = Lengths[-5:-4]
win_len = 10
for l1 in to_analyse:
    num_bins = np.arange(0,int(l1/scale_CNIH2_protein),int(bin_size/scale_CNIH2_protein))
    raw_MAP_norm_GluA2 = np.divide(raw_Glua2_data[l1],raw_MAP2_data[l1])
    binned_sum = np.add.reduceat(raw_MAP_norm_GluA2,num_bins,axis=1)
    # breakpoint()
    mean_sw = np.median(binned_sum,axis=0)[:-1]
    std_sw = binned_sum.std(axis=0)[:-1]
    x = np.linspace(0,l1,mean_sw.shape[0])
    x_range = bins[0:num_bins.shape[0]]
    # print(spearmanr(x,mean_sw))
    # breakpoint()
    print(binned_sum.mean(axis=0))
    fig,ax = plt.subplots(figsize=(8,6),ncols=1,nrows=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylabel("MAP2 normalized \n GluA2 [a.b.u]", fontsize=fsize)
    ax.set_xlabel("Length of dendrite ($\mu$m) ", fontsize=fsize)
    # for i in range(binned_sum.shape[-1]):
    flyprops = {'markersize': 0.01}
    colorprops = None
    meanprops = {'marker':"s",'markersize':plt_widget.msize,'color':"#2ca02c"}
    medianprops =meanprops# {'linewidth': 0}
    whiskerprops = {'color': '#c0c0c0ff','linewidth':7}
    boxprops = dict(facecolor='k', color='k')
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(2)
    ax.boxplot(binned_sum[:,:-1],positions = bins[0:num_bins.shape[0]-1]+bin_size/2,
               showmeans=False,
               showfliers=False, showcaps=False,
               whiskerprops=whiskerprops, boxprops=boxprops, medianprops = medianprops,
               widths=1,patch_artist=True)
    yi_fit,rseq,param = ExpFitWithMinimize("2E",x,mean_sw,std_sw,0,+1,"Glua2")
    print("drop =",1-yi_fit[-1]/yi_fit[0])
    # breakpoint()
    x_tics = np.arange(0,l1+1,25)
    no_labels = x_tics.shape[0]#int(np.floor(l1 / 25))
    label = [f'{i * l1 / no_labels:.0f}' for i in range(no_labels + 1)]
    ax.set_xticks(range(0, l1 + 1, 20))
    ax.set_xticklabels(label)
    print(x_tics)
    plt.xticks(x_tics)
    lw=3
    # y_tics = np.arange(0,2,0.5)
    # plt.yticks(y_tics)
    # ax.plot(x,yi_fit,color = COLORS_dict["shaft_s"],linestyle="--",linewidth=lw,
    #         label="2E-fit")
    # for i in range(raw_MAP_norm_GluA2.shape[0]):
    #     if np.max(raw_MAP_norm_GluA2[i]> 3):
    #         print("i_val = ",i)
    #     ax.plot(x, raw_MAP_norm_GluA2[i], color='gray', linewidth=lw / 3, alpha=0.1)
    # ax.plot(x,mean_sw,color=COLORS_dict["shaft_i"],label="Glua2",linewidth=lw)
    ax.set_xlim([-1, l1+1])
plt.legend()
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder, "raw_glua2_density_profile_{}".format(l1)))
plt.show()

for l1 in to_analyse:
    num_bins = np.arange(0,int(l1/scale_CNIH2_protein),int(bin_size/scale_CNIH2_protein))
    raw_MAP_norm_CNIH2 = np.divide(raw_CNIH2_data[l1],raw_MAP2_data[l1])
    binned_sum = np.add.reduceat(raw_MAP_norm_CNIH2,num_bins,axis=1)
    # breakpoint()
    mean_sw = np.median(binned_sum,axis=0)[:-1]
    std_sw = binned_sum.std(axis=0)[:-1]
    x = np.linspace(0,l1,mean_sw.shape[0])
    x_range = bins[0:num_bins.shape[0]]
    # print(spearmanr(x,mean_sw))
    # breakpoint()
    print(binned_sum.mean(axis=0))
    fig,ax = plt.subplots(figsize=(8,6),ncols=1,nrows=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylabel("MAP2 normalized \n CNIH2 [a.b.u]", fontsize=fsize)
    ax.set_xlabel("Dendritic distance ($\mu$m) ", fontsize=fsize)
    # for i in range(binned_sum.shape[-1]):
    flyprops = {'markersize': 0.01}
    colorprops = None
    meanprops = {'marker':"s",'markersize':plt_widget.msize,'color':"#2ca02c"}
    medianprops =meanprops# {'linewidth': 0}
    whiskerprops = {'color': '#c0c0c0ff','linewidth':7}
    boxprops = dict(facecolor='k', color='k')
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(2)
    ax.boxplot(binned_sum[:,:-1],positions = bins[0:num_bins.shape[0]-1]+bin_size/2,
               showmeans=False,
               showfliers=False, showcaps=False,
               whiskerprops=whiskerprops, boxprops=boxprops, medianprops = medianprops,
               widths=1,patch_artist=True)
    yi_fit,rseq,param = ExpFitWithMinimize("2E",x,mean_sw,std_sw,0,+1,"CNIH2")
    print("drop =",1-yi_fit[-1]/yi_fit[0])
    # breakpoint()
    x_tics = np.arange(0,l1+1,25)
    no_labels = x_tics.shape[0]#int(np.floor(l1 / 25))
    label = [f'{i * l1 / no_labels:.0f}' for i in range(no_labels + 1)]
    ax.set_xticks(range(0, l1 + 1, 20))
    ax.set_xticklabels(label)
    print(x_tics)
    plt.xticks(x_tics)
    lw=3
    # y_tics = np.arange(0,2,0.5)
    # plt.yticks(y_tics)
    # ax.plot(x,yi_fit,color = COLORS_dict["shaft_s"],linestyle="--",linewidth=lw,
    #         label="2E-fit")
    # for i in range(raw_MAP_norm_GluA2.shape[0]):
    #     if np.max(raw_MAP_norm_GluA2[i]> 3):
    #         print("i_val = ",i)
    #     ax.plot(x, raw_MAP_norm_GluA2[i], color='gray', linewidth=lw / 3, alpha=0.1)
    # ax.plot(x,mean_sw,color=COLORS_dict["shaft_i"],label="Glua2",linewidth=lw)
    ax.set_xlim([-1, l1+1])
plt.legend()
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder, "raw_cnih2_density_profile_{}".format(l1)))
plt.show()

# breakpoint()
fig,ax2 = plt.subplots(figsize=(10,8*len(to_analyse)),ncols=1,nrows=len(to_analyse))
# lta = 100
ax2 = [ax2]
fitted_vals = []
for ldx,l1 in enumerate(to_analyse):
    x_range =  np.arange(0, l1, bin_size)
    # ax2[ldx].plot(x_range,mean_CNIH2_density[l1],color=COLORS_dict["shaft_s"],label="CNIH2",marker='d',linestyle="--")
    # ax2[ldx].plot(x_range,mean_MAP2_density[l1],color="r",label="MAP2",marker='s',linestyle="-")
    # ax2[ldx].spines[["right", "top"]].set_visible(False)
    ax2[ldx].errorbar(x_range,mean_Glua2_density[l1],sem_Glua2_density[l1],color=COLORS_dict["shaft_i"],
                      label="Experimental data",marker='o',markersize=7,linestyle="--",alpha=0.8)
    # ax2[ldx].plot(x_range, mean_MAP2_density[l1], color=COLORS_dict["shaft_s"], label="MAP2", marker='o',
    #               markersize=7, linestyle="--", alpha=0.8)

    # breakpoint()

    x_lit, yi_lit = Total_AMPAR.RunSSProtein(D_P=1., v_P=1.46, x_range=[0, l1])
    ax2[ldx].plot(x_lit, yi_lit/yi_lit[0],
                  'g-',
                  label=r"$D_p = {:.2f}$ $\mu m^2/s$, $v_P = {:.2f}$ $\mu m/s$".
                  format(1., 1.46))

    x_lit, yi_lit = Total_AMPAR.RunSSProtein(D_P=1.5e-2, v_P=0, x_range=[0, l1])
    # breakpoint()
    ax2[ldx].plot(x_lit, yi_lit/yi_lit[0],
                  'y-',
                  label=r"$D_p = {:.3f}$ $\mu m^2/s$, $v_P = {:.2f}$ $\mu m/s$".
                  format(1.5e-2, 0))
    # x_min = x_min[int(0)]
    x, yi_fit, chi_squ, paras, mini, out2 = FitModelProtein(x_range, mean_Glua2_density[l1] / mean_Glua2_density[l1][0],
                                                            sem_Glua2_density[l1],
                                                            molecule="GluA2")  # curve_fit(exp_fit, x2, d2)
    print("chi squ = ",chi_squ)
    # breakpoint()
    print(yi_fit[0]-yi_fit[-1])
    fitted_vals = [10 ** paras['D_P'].value,10 ** paras['v_P'].value]
    x_lit, yi_fit = Total_AMPAR.RunSSProtein(D_P=fitted_vals[0], v_P=fitted_vals[1])
    # breakpoint()
    left, bottom, width, height = [0.4, 0.45, 0.3, 0.3]
    tics = np.arange(0, 1.2, 0.5)
    ax1 = ax2[ldx].inset_axes(bounds=[left, bottom, width, height], zorder=4)
    # ax1.set_title("Normalized \n density")
    ax1.set_ylabel("Normalized \n protein density")
    # ax1.set_yticks(tics)
    ax.set_xlim(x_lit[0], x_lit[-1] + 1)
    x_tics = np.arange(x_lit[0], x_lit[-1] + 1, (x_lit[-1] - x_lit[0] + 1) // 5)
    # ax.set_xticks(x_tics)
    ax1.plot(x_lit, yi_fit / yi_fit[0], color=COLORS_dict["shaft_i"])
    ax2[ldx].plot(x_lit,yi_fit/yi_fit[0],
                  color=COLORS_dict["shaft_i"],
                  linestyle="-",
                  label=r"Theory fit: $D_p = {:.2f}$ $\mu m^2/s$, $v_P = {:.4f}$ $\mu m/s$".
                  format(10 ** paras['D_P'].value,10 ** paras['v_P'].value))

    # breakpoint()
    ax2[ldx].set_ylabel("Normalized density",fontsize=fsize)
    ax2[ldx].set_xlabel("Length of dendrite ($\mu$m) ",fontsize=fsize)
    ax2[ldx].set_title("GluA2 protein",fontsize=fsize)
    ax2[ldx].legend(loc='lower left',fontsize=fsize)
    ax2[ldx].set_xlim([-1,l1+1])
    ax2[ldx].set_ylim([0.0,1.1])


    # ax3[ldx].set_xlim([-1, l1])
    # ax3[ldx].set_ylim([0.0, 1.1])
#
# plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Normlized_intesity_Dend_dist_GluA2_fitted"))
plt.show()


# fitting to the total GluA2 from MK data
fig,ax2 = plt.subplots(figsize=(10,8*len(to_analyse)),ncols=1,nrows=len(to_analyse))
ax2 = [ax2]
for l1 in to_analyse:
    x_range = np.arange(0, l1, bin_size)
    ax2[ldx].errorbar(x_range, mean_total_Glua2_mk[l1], sem_total_Glua2_mk[l1], color=COLORS_dict["shaft_i"],
                      label="MK data", marker='o', markersize=7, linestyle="--", alpha=0.8)

    x, yi_fit, chi_squ, paras, mini, out2 = FitModelProtein(x_range, mean_total_Glua2_mk[l1] / mean_total_Glua2_mk[l1][0],
                                                            sem_total_Glua2_mk[l1],
                                                            molecule="GluA2")  # curve_fit(exp_fit, x2, d2)
    print("chi_seq = ", chi_squ)
    # breakpoint()
    print("Theory fit: $D_p = {}$ $\mu m^2/s$, $v_P = {}$ $\mu m/s$".
                  format(10 ** paras['D_P'].value, 10 ** paras['v_P'].value))
    print(yi_fit[0] - yi_fit[-1])
    # fitted_vals = [10 ** paras['D_P'].value, 10 ** paras['v_P'].value]
    x_lit, yi_fit = Total_AMPAR.RunSSProtein(D_P=10 ** paras['D_P'], v_P=10 ** paras['v_P'].value)
    # breakpoint()
    left, bottom, width, height = [0.4, 0.45, 0.3, 0.3]
    tics = np.arange(0, 1.2, 0.5)
    ax1 = ax2[ldx].inset_axes(bounds=[left, bottom, width, height], zorder=4)
    # ax1.set_title("Normalized \n density")
    ax1.set_ylabel("Normalized \n protein density")
    # ax1.set_yticks(tics)
    ax.set_xlim(x_lit[0], x_lit[-1] + 1)
    x_tics = np.arange(x_lit[0], x_lit[-1] + 1, (x_lit[-1] - x_lit[0] + 1) // 5)
    # ax.set_xticks(x_tics)
    ax1.plot(x_lit, yi_fit / yi_fit[0], color=COLORS_dict["shaft_i"])
    ax2[ldx].plot(x_lit, yi_fit / yi_fit[0],
                  color=COLORS_dict["shaft_i"],
                  linestyle="-",
                  label=r"Theory fit: $D_p = {:.2f}$ $\mu m^2/s$, $v_P = {:.4f}$ $\mu m/s$".
                  format(10 ** paras['D_P'].value, 10 ** paras['v_P'].value))

    # breakpoint()
    ax2[ldx].set_ylabel("Normalized density", fontsize=fsize)
    ax2[ldx].set_xlabel("Length of dendrite ($\mu$m) ", fontsize=fsize)
    ax2[ldx].set_title("GluA2 protein", fontsize=fsize)
    ax2[ldx].legend(loc='lower left', fontsize=fsize)
    ax2[ldx].set_xlim([-1, l1 + 1])
    ax2[ldx].set_ylim([0.0, 2])
plt_widget.SaveFigures(os.path.join(op_folder,"Normlized_intesity_Dend_MK_dist_GluA2_fitted"))
plt.show()
    # breakpoint()
def plot_velocity_impact(dp,vp,color,fname,prefix="",log=False,y_lim=False,loc="upper right"):
    fig,ax = plt.subplots(figsize=(10,8),ncols=1,nrows=1)
    x_lit, yi_lit = Total_AMPAR.RunSSProtein(dp,vp)
    print((1-yi_lit[-1]/yi_lit[0])*100)
    left, bottom, width, height = [0.4,0.45,0.3,0.3]
    tics = np.arange(0,1.2,0.5)
    ax1 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)
    # ax1.set_title("Normalized \n density")
    ax1.set_ylabel("Normalized \n protein density")
    # ax1.set_yticks(tics)
    ax.set_xlim(x_lit[0],x_lit[-1]+1)
    x_tics = np.arange(x_lit[0],x_lit[-1]+1,(x_lit[-1]-x_lit[0]+1)//5)
    # ax.set_xticks(x_tics)
    ax1.plot(x_lit,yi_lit/yi_lit[0],color=color)
    ax.plot(x_lit, yi_lit,
            color=color,
            linestyle='-',
            label=r"{}$D_p = {:.2f}$ $\mu m^2/s$, $v_P = {:.4f}$ $\mu m/s$".format(prefix,dp,vp))

    ax.set_xlabel("Length of dendrite ($\mu$m) ", fontsize=fsize)
    ax.set_title("GluA2 protein", fontsize=fsize)
    ax.set_ylabel("Protein density", fontsize=fsize)
    ax.legend(loc=loc, fontsize=fsize)

    if log:
        plt.yscale("log")
        ax1.set_yscale("log")
        # ax.set_ylim([0,int(1.1*np.max(np.log(yi_lit))) ])
    # else:
    ax.set_ylim([0, 1.1 * np.max(yi_lit)])
    if y_lim:
        ax1.set_ylim([0,1.1])
    plt_widget.SaveFigures(os.path.join(op_folder,"{}_{}_{}".format(fname,dp,vp)))
    plt.show()

plot_velocity_impact(1.0,1.46,'g',"Normlized_intesity_Dend_dist_GluA2_full_",log=True,loc="upper left")
plot_velocity_impact(1e-2,0,'y',"Normlized_intesity_Dend_dist_GluA2_full_",y_lim=True)
plot_velocity_impact(fitted_vals[0], fitted_vals[1],COLORS_dict["shaft_i"],
                     "Normlized_intesity_Dend_dist_GluA2_full_",prefix="Theory fit: \n",y_lim=False)


"""

fig, ax1 = plt.subplots(figsize=(8, 30), ncols=1, nrows=len(to_analyse))
# x_100 =  np.arange(0, 100, bin_size)
for ldx, l1 in enumerate(to_analyse):
    x_range = np.arange(0, l1, bin_size)
    ax1[ldx].plot(x_range, mean_CNIH2_density[l1], color=COLORS_dict["shaft_s"],
                      label="CNIH2", marker="o",linestyle="-",alpha=0.5)
    # ax1[ldx].errorbar(x_range+bin_size/2,mean_Glua2_density[l1],sem_Glua2_density[l1],color=COLORS_dict["shaft_i"],label="Glua2")

    # popt, pcov = curve_fit(line_fit, x_range+bin_size/2, mean_CNIH2_density[l1])
    popt_e, pcov_e = curve_fit(exp_fit, x_range, mean_CNIH2_density[l1])

    bp = int(12 / bin_size)  # untill 10 microns
    d1, d2 = np.split(mean_CNIH2_density[l1], [bp])
    x1, x2 = np.split(x_range, [bp])
    print(x1,x2)
    sem1, sem2 = np.split(sem_CNIH2_density[l1], [bp])
    popt_e1, pcov_e1 = curve_fit(exp_fit, x1, d1)
    x3,yi_fit,chi_squ,paras,mini,out2 = FitModelProtein(x2,d2/d2[0],sem2)#curve_fit(exp_fit, x2, d2)
    # print(popt_e1, popt_e2, popt_e)
    mean_fit_e1 = exp_fit(x1, *popt_e1)
    # mean_fit_e2 = exp_fit(x2, *popt_e2)
    mean_fit_exp = exp_fit(x_range, *popt_e)
    # r2 = R_seq(mean_fit,mean_CNIH2_density[l1])
    r2_exp = ChiSq(mean_CNIH2_density[l1], mean_fit_exp, sem_CNIH2_density[l1])
    r2_e1 = ChiSq(d1, mean_fit_e1, sem1)
    # r2_e2 = ChiSq(d2, mean_fit_e2, sem2)
    ax1[ldx].plot(x_range, mean_fit_exp, 'r--', label=r"exp-fit,$\chi^2 = %.2f$" % r2_exp)
    ax1[ldx].plot(x1, mean_fit_e1, 'g-', label=r"exp1-fit,$\chi^2 = %.2f$" % r2_e1)
    ax1[ldx].plot(x3+x2[0], d2[0]*yi_fit, 'k-', label=r"exp2-fit,$\chi^2 = %.2f$" % chi_squ)
    # ax1[ldx].plot(x_range+bin_size/2,mean_fit_exp,'g--',label=r"line-fit,$R^2 = %.2f$" % r2)
    # f_size = 12
    loc = "upper right"
    # print(pearsonr(mean_Glua2_density[l1],mean_CNIH2_density[l1]))
    # ax3.errorbar(x_100,mean_int_density,sem_int_density,color=COLORS_dict["shaft_i"],label="Glua2")
    ax1[ldx].set_ylabel("Normalized CNIH2 fluorescent", size=fsize)
    ax1[ldx].set_xlabel(r"Dendritic distance in $\mu m$, N={}".format(CNIH2_length_data[l1].shape[0]), size=fsize)
    ax1[ldx].set_title("Dendritic distribution of GFP normalized CNIH2 density", size=fsize)
    ax1[ldx].legend(fontsize=fsize, loc=loc)
    ax1[ldx].set_xlim([-1, l1])
    # ax1[ldx].yaxis.tick_right()
plt.tight_layout()
plt_widget.SaveFigures(os.path.join(op_folder,"Normlized_density_Dend_dist_CNIH2_exp_fitted"))
plt.show()

l1 = 100
x_range = np.arange(0,l1,bin_size)
cnih2_mean = mean_CNIH2[l1]
cnih2_sem = sem_CNIH2_density[l1]
cnih2_m_den = mean_CNIH2_density[l1]

plt_widget.PlotBinnedStats(np.asarray([x_range]), np.asarray([cnih2_m_den]), np.asarray([sem_CNIH2_density[l1]]),
                           np.asarray([cnih2_m_den]), ["cnih2"], x_lab, y_lab, ["Norm CNIH2 density"], ["#DDA15E"],"",
                           op_folder+"_norm__with_CNIH2_ss_mode_fit",
                                    bin_size,save_it = 0,fit_exp=3,in_set=0,set_axis_label=0)
                #
plt_widget.PlotBinnedStats(np.asarray([x_range]), np.asarray([cnih2_m_den]), np.asarray([sem_CNIH2_density[l1]]),
                           np.asarray([cnih2_m_den]), ["cnih2"], x_lab, y_lab, ["Norm CNIH2 density"], ["#DDA15E"],"",
                           op_folder+"_norm__with_CNIH2_ss_mode_fit",
                                    bin_size,save_it = 0,fit_exp=1,in_set=0,set_axis_label=0,exp_method="2E")

"""
