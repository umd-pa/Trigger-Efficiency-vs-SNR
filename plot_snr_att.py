import os, sys, shutil, glob
import argparse
import math
import numpy as np
import ROOT
import pickle
import scipy
import matplotlib.pyplot as plt
from ROOT import gStyle, gPad, kRed

# load the RNO-G library
ROOT.gSystem.Load(os.environ.get('RNO_G_INSTALL_DIR')+"/lib/libmattak.so")


#path to files with snr arrays for each attenuation 

indir = "/data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/snr_npy_23_3_switch"


files = sorted(glob.glob(os.path.join(indir, "*")))


#calculating mean snr value in each attenuation bin 

snr_means = []
snr_std = []
atts = []
atts_12 = []
snr_means_12 = []
for f in files:
    att_one = f.split("/")[-1].split("_")[-1].split(".")[0]
    att_two = f.split("/")[-1].split("_")[-1].split(".")[1]
    att = int(att_one) + 0.1*int(att_two)
    
    snr_arr = np.load(f)

    if (len(snr_arr) != 0):
        atts.append(att)
        snr_means.append(np.average(snr_arr))
        snr_std.append(np.std(snr_arr))
        if (att >= 6 and att <= 15):
            atts_12.append(att)
            snr_means_12.append(np.average(snr_arr))
        

#exponential fit for attenuation between 6-15 dB

fit = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  atts_12, snr_means_12,  p0=(17.5, -0.2))

x_vals = np.arange(0,32, 0.5)
y_vals = []

for x in x_vals:
    y_vals.append(fit[0][0]*np.exp(fit[0][1]*x))


#plot snr vs attenuation 

plt.scatter(atts,snr_means, s = 15)
plt.plot(x_vals, y_vals)
plt.errorbar(atts,snr_means, yerr = snr_std, ls = "None")
plt.title("Event SNR vs Attenuation")
plt.xlabel("Attenuation (dB)")
plt.ylabel("Event SNR")
plt.legend()
plt.savefig("snr_vs_att_23_3_switch.png")



    
