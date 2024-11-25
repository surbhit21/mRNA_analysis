#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14

@author: surbhitwagle
"""
import CNIH2_TemporalIntegration
import matplotlib.pyplot as plt
import os.path
from Utility import *


soma_tra = "11_11_2024_11_56_23"
brief_dend_tra = "11_11_2024_10_24_43"
prolonged_dend_tra = "11_08_2024_16_13_52"
per = "median"
labels = ["CNIH2"]
soma_tra_median = np.load(os.path.join(os.getcwd(), "Scripts/CNIH2-translation/{0}/CNIH2_{0}_median.npy".format(soma_tra)))
brief_dend_tra_median = np.load(os.path.join(os.getcwd() , "Scripts/CNIH2-translation/{0}/CNIH2_{0}_median.npy".format(brief_dend_tra)))
prolonged_dend_tra_median = np.load(os.path.join(os.getcwd() , "Scripts/CNIH2-translation/{0}/CNIH2_{0}_median.npy".format(prolonged_dend_tra)))

soma_tra_timepoints = np.load(os.path.join(os.getcwd() , "Scripts/CNIH2-translation/{0}/timepoints_{0}_100_percent.npy".format(soma_tra)))
brief_dend_tra_timepoints = np.load(os.path.join(os.getcwd() , "Scripts/CNIH2-translation/{0}/timepoints_{0}_100_percent.npy".format(brief_dend_tra)))
prolonged_dend_timepoints = np.load(os.path.join(os.getcwd() , "Scripts/CNIH2-translation/{0}/timepoints_{0}_100_percent.npy".format(prolonged_dend_tra)))
op_folder = os.path.join(os.getcwd(),"Figures")

plt.rc('font', family='Arial')
f_size= 30
dpi = 300
def SaveFigures(filename, ext_list=[".png", ".svg", ".pdf"], dpi=300):
    """

        function to save figures
        required arguments:
            filename
    """
    for ext in ext_list:
        plt.savefig(filename + ext, dpi=dpi)


fig, ax = plt.subplots(figsize=(8, 6), ncols=1)
ax.plot(soma_tra_timepoints/60,soma_tra_median/soma_tra_median[0],"#ADB5BD",label = "Soma syn. ")
ax.plot(brief_dend_tra_timepoints/60,brief_dend_tra_median/brief_dend_tra_median[0],"#6C757D",label = "local syn. (Brief)")
ax.plot(prolonged_dend_timepoints/60,prolonged_dend_tra_median/prolonged_dend_tra_median[0],"#2ca02c",label = "local syn. (Prolonged)")
ax.vlines(x=60,ymin=1,ymax=6,color = "r",linestyle="--")
ax.set_xlabel("Time (Mins)",fontsize=f_size)
ax.set_ylabel("Normalized median change",fontsize=f_size)
plt.legend(frameon=False,fontsize=f_size,loc="upper right",labelcolor='linecolor')
SaveFigures("{}/Median_CNIH2_comparison".format(op_folder))
plt.show()

