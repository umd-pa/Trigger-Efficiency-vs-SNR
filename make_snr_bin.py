import os
import argparse
import math
import numpy as np
import ROOT
import pickle
import matplotlib.pyplot as plt
from ROOT import gStyle, gPad, kRed
import scipy
import scipy.optimize as opt
from scipy import signal
from scipy.fft import fft, fftfreq, rfft, irfft
from array import array
from ROOT import gStyle, gPad, kRed, TMath
from NuRadioReco.utilities import bandpass_filter
from NuRadioReco.utilities import fft as fft_reco
from NuRadioReco.detector.RNO_G import analog_components

# load the RNO-G library
ROOT.gSystem.Load(os.environ.get('RNO_G_INSTALL_DIR')+"/lib/libmattak.so")

# make sure we have enough arguments to proceed
parser = argparse.ArgumentParser(description='daqstatus example')
parser.add_argument('--file', dest='file', required=True)
parser.add_argument('--att', type=float, required=True)
parser.add_argument('--ch', type=int, required=True)
args = parser.parse_args()
filename = args.file
att_lim = args.att
ch = args.ch

#voltage calibration coeffs

cal_path = "/data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/volCalConsts_pol9_s23_1697181551-1697183024.root"
fIn = ROOT.TFile.Open(filename)
combinedTree = fIn.Get("combined")


volCalib = ROOT.mattak.VoltageCalibration()
volCalib.readFitCoeffsFromFile(cal_path)

d = ROOT.mattak.DAQStatus()
wf = ROOT.mattak.Waveforms()
hdr = ROOT.mattak.Header()

combinedTree.SetBranchAddress("daqstatus", ROOT.AddressOf(d))
combinedTree.SetBranchAddress("waveforms", ROOT.AddressOf(wf))
combinedTree.SetBranchAddress("header", ROOT.AddressOf(hdr))

num_events = combinedTree.GetEntries()

rms_all = []


#real low pass filter 

def get_filter(frequencies):
    temp = 298.15
    amp_response = analog_components.load_amp_response(amp_type = 'phased_array', temp = 298.15)
    amp_response = amp_response['gain'](frequencies, temp) * amp_response['phase'](frequencies)

    return amp_response

#high pass filter (butterworth) to eliminate noise 

def apply_butterworth(spectrum, frequencies, passband, order=8):
    f = np.zeros_like(frequencies, dtype=complex)
    mask = frequencies > 0
    b, a = scipy.signal.butter(order, passband, 'highpass', analog=True)
    w, h = scipy.signal.freqs(b, a, frequencies[mask])
    f[mask] = h
    filtered_spectrum = f * spectrum

    return filtered_spectrum

#approximation of the real low pass filter 

def apply_butterworth_low(spectrum, frequencies, passband, order=8):
    f = np.zeros_like(frequencies, dtype=complex)
    mask = frequencies > 0
    b, a = scipy.signal.butter(order, passband, 'lowpass', analog=True)
    w, h = scipy.signal.freqs(b, a, frequencies[mask])
    f[mask] = h
    filtered_spectrum = f * spectrum

    return filtered_spectrum


for event in range(num_events):
    combinedTree.GetEntry(event)

    sysclk = hdr.sysclk
    sysclk_last_pps = hdr.sysclk_last_pps
    sys_diff = sysclk - sysclk_last_pps   
    
    #calculating the noise rms 

    if (hdr.trigger_info.force_trigger == True):
        c = ROOT.mattak.CalibratedWaveforms(wf, hdr, volCalib, False)


        g = c.makeGraph(ch)
    
        voltage = g.GetY()
        time = g.GetX()
        
        
        voltage_0 = []

        for v in voltage:
            voltage_0.append(v - voltage[0])
        
        n_samples = len(time)
        sampling_frequency = 1/(time[1] - time[0])
        pb = 0.23 #for the butterworth low pass filter
        pb2 = 0.05  #for the high pass filter

        spectrum = fft_reco.time2freq(voltage_0, sampling_frequency)
        frequencies = np.fft.rfftfreq(n_samples, 1 / sampling_frequency)
        
        filtered_hi = apply_butterworth(spectrum, frequencies, pb2, 8)
        filtered_spectrum = filtered_hi*get_filter(frequencies)
        filtered_trace = fft_reco.freq2time(filtered_spectrum, sampling_frequency)

        
        rms = np.sqrt(np.mean(filtered_trace**(2)))
        rms_all.append(rms)
        
        del g
        del c

rms_avg = np.average(rms_all)


snr_arr = []


for event in range(num_events):
    combinedTree.GetEntry(event)

    sysclk = hdr.sysclk
    sysclk_last_pps = hdr.sysclk_last_pps
    sys_diff = sysclk - sysclk_last_pps

    att = d.calinfo.attenuation
    
    #calculating snr at a particular attenuation (att_lim)

    if att == att_lim:
        if (sys_diff <= 200*10**(3)):
            c = ROOT.mattak.CalibratedWaveforms(wf, hdr, volCalib, False)
            g = c.makeGraph(ch)
            
            voltage = g.GetY()
            time = g.GetX()
        
            voltage_0 = []

            for v in voltage:
                voltage_0.append(v - voltage[0])

            n_samples = len(time)
            sampling_frequency = 1/(time[1] - time[0])
            pb = 0.23
            pb2 = 0.05

            spectrum = fft_reco.time2freq(voltage_0, sampling_frequency)
            frequencies = np.fft.rfftfreq(n_samples, 1 / sampling_frequency)

            #applying high pass filter followed by real low pass filter 

            filtered_hi = apply_butterworth(spectrum, frequencies, pb2, 8)
            filtered_spectrum = filtered_hi*get_filter(frequencies)
            filtered_trace = fft_reco.freq2time(filtered_spectrum, sampling_frequency)
        
            #snr calculation
            
            high = np.max(filtered_trace)
            low = np.min(filtered_trace)
            snr = (high-low)/2./rms_avg
                
            del g 

            snr_arr.append(snr)





#save an array of snr values for each attenuation 

np.save(f'/data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/snr_npy_23_{ch}_switch/snr_{att_lim}.npy', snr_arr)
