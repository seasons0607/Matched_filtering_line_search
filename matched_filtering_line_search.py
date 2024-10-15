# -*- coding: utf-8 -*-
# Last modified 2023-07-06
from astropy.io import fits
from matplotlib import gridspec
from astropy.stats import bayesian_blocks
import os
import sys
import numpy as np
import pandas as pd 
from astropy.time import Time
import matplotlib.pyplot as plt 
import yaml
import statistics
import math
from xspec import *
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib import ticker
from tqdm import tqdm
import subprocess

def gaussian(E, E_i, sigma):
    return (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((E - E_i) / sigma)**2)

def resolution_sigma(E, a0, b0, c0, d0):
    return (a0*(E**3.0)+b0*(E**2.0)+c0*E+d0)/(2.35*1000.0) #keV

def matched_filter(energy, energy_err, rate, a1, b1, c1, d1):
    N = len(energy)  # The number of energy bins
    filtered_rate = np.zeros(N)  
    
    for i in range(N):
        E_i = energy[i]
        sum_value = 0
        for j in range(N):
            E_j = energy[j]
            sigma = resolution_sigma(E_i, a1, b1, c1, d1)
            if abs(E_j - E_i) <= 3 * sigma:
                sum_value += rate[j] * gaussian(E_j, E_i, sigma) * (2.0 * energy_err[j])  
        filtered_rate[i] = sum_value
    
    return filtered_rate

def extract_integers(start, end):
    start = int(start) if start == int(start) else int(start) + 1
    end = int(end)
    return list(range(start, end + 1))

def extract_integers_as_strings(start, end):
    start = int(start) if start == int(start) else int(start) + 1
    end = int(end)
    integer_list = list(range(start, end + 1))
    string_list = [str(i) for i in integer_list]
    return string_list

def calc_fwhm(rmf_path):
    #rmf_path = "/Volumes/crx_raid0/inoue_22/projects/nicer_systematic/RS_CVn/Sigma_Gem/1200040104/analysis/spec/block0/fake_spectra/1200040104_block0.rmf"
    rmf_name = os.path.splitext(os.path.basename(rmf_path))[0]

    fitsFile_rmf = fits.open(rmf_path)
    fitsFile_rmf_header = fitsFile_rmf[2].header

    fitsFile_rmf_data = fitsFile_rmf[1].data
    energy_binmin = fitsFile_rmf_data['E_MIN']
    energy_binmax = fitsFile_rmf_data['E_MAX']
    energy_bin = []
    for i in range(len(energy_binmin)):
        energy_bin.append((energy_binmin[i]+energy_binmax[i])/2.0)

    fitsFile_rmf_data = fitsFile_rmf[2].data
    energy_lo = fitsFile_rmf_data['ENERG_LO']
    energy_hi = fitsFile_rmf_data['ENERG_HI']
    rmfmatrix = fitsFile_rmf_data['MATRIX']

    energy_all = [0.0111*x for x in range(200, 1100, 2)]
    fwhm_all = []
    for i in tqdm(range(len(energy_all)), desc="Calculating FWHM of "+str(rmf_name)+".rmf"):
        energy = energy_all[i]
        index = [i for i, (lo, hi) in enumerate(zip(energy_lo, energy_hi)) if lo < energy < hi][0]
        fwhm_energy_list =[]
        for j in range(len(energy_binmin)):
            if rmfmatrix[index][j] > 0.5 * (max(rmfmatrix[index])):
                fwhm_energy_list.append([energy_binmin[j], energy_binmax[j]])
        
        def MATRIX_fiting_function(x, amplitude, mean, sigma):
            return amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))
        
        initial_guess = [max(rmfmatrix[index]), energy, 0.1]
        popt, pcov = curve_fit(MATRIX_fiting_function, energy_bin, rmfmatrix[index], p0=initial_guess)
        amplitude_fit, mean_fit, sigma_fit = popt
        
        fwhm = 2.35 * sigma_fit * 1000.0 #eV
        fwhm_all.append(fwhm) #eV

    #Plotting 
    fig=plt.figure(figsize=(11, 7))
    gs = gridspec.GridSpec(1, 1)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    axs = []
    for i in range(1):
        axs.append(fig.add_subplot(gs[i, 0])) 
    label_fontsize = 20.5
    legend_size = 17
    markersize = 0.0
    linewidth_data = 2.0
    linewidth_model = 2.7
    str_fontsize = 23
    
    # axs[0].plot(energy_all, fwhm_all, marker="s")
    axs[0].tick_params(direction = "in", right = True, labelright = False, top = True, labeltop = False, labelleft=True, labelbottom = True, labelsize=label_fontsize)
    axs[0].set_ylabel(r"FWHM (eV)", fontsize = label_fontsize+2.5)
    axs[0].set_xlabel(r"Energy (keV)", fontsize = label_fontsize+2.5)
    axs[0].set_xlim(min(energy_all), max(energy_all))
    axs[0].set_ylim(min(fwhm_all), max(fwhm_all))
    axs[0].minorticks_on()
    axs[0].grid(False, which="minor")
    axs[0].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='x')
    axs[0].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='y')

    def FWHM_fiting_function(x, a, b, c, d):
        return a*(x**3.0)+b*(x**2.0)+c*x+d
    popt, pcov = curve_fit(FWHM_fiting_function, energy_all, fwhm_all, p0=[1.0, 1.0, 1.0, 1.0])
    xd = np.arange(1.0, 13.0, 0.001)
    estimated_curve = FWHM_fiting_function(xd, *popt)

    axs[0].plot(xd, estimated_curve, color="navy", linewidth=linewidth_model)
    axs[0].text(0.025, 0.92, "NICER", fontsize=str_fontsize+5, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[0].transAxes, fontweight = "bold", ha="left")
    axs[0].text(0.025, 0.85, rmf_name, fontsize=str_fontsize+5, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[0].transAxes, fontweight = "bold", ha="left")


    plt.subplots_adjust(hspace=0, left=0.105, right = 1-0.105, top=0.985, bottom=0.1)
    plt.savefig(str(rmf_name)+"_FWHM.pdf", format="pdf", dpi=1000)

    return popt

def plot_spectrum_with_MC(rmf_path, arf_path, spectrum_txt_path, continuum_txt_path, continuum_xcm_path, exposure, trial_number, continuum_energy_low, continuum_energy_upp, a, b, c, d):
    #### Reading data ####
    spectrum_data = []
    with open(spectrum_txt_path,"r")as fin:
        for line in fin.readlines():
            row = []
            toks = line.split(' ')
            for tok in toks:
                try:
                    num = float(tok)
                except ValueError:
                    #print(e, file=sys.stderr)
                    continue
                row.append(num)
    
            spectrum_data.append(row)
    spectrum_energy = []
    spectrum_energy_err = []
    spectrum_rates = []
    spectrum_rates_err = []

    
    for i in range(len(spectrum_data)):
        spectrum_energy.append(spectrum_data[i][0])
        spectrum_energy_err.append(spectrum_data[i][1])
        spectrum_rates.append(spectrum_data[i][2])
        spectrum_rates_err.append(spectrum_data[i][3])

    continuum_data = []
    with open(continuum_txt_path,"r")as fin:
        for line in fin.readlines():
            row = []
            toks = line.split(' ')
            for tok in toks:
                try:
                    num = float(tok)
                except ValueError:
                    #print(e, file=sys.stderr)
                    continue
                row.append(num)
    
            continuum_data.append(row)
    continuum_model_all = []
    
    for i in range(len(continuum_data)):
        continuum_model_all.append(continuum_data[i][4])
    
    #### Generating fake spectra ####
    xcm_dir = os.path.dirname(continuum_xcm_path)
    xcm_name = os.path.basename(continuum_xcm_path)
    rmf_name = os.path.basename(rmf_path)
    arf_name = os.path.basename(arf_path)

    if os.path.isdir(xcm_dir+"/fake_spectra") == False:
        os.mkdir(xcm_dir+"/fake_spectra")
    if os.path.isfile(xcm_dir+"/fake_spectra/"+rmf_name) == False:
        os.system("cp "+rmf_path+" "+xcm_dir+"/fake_spectra/"+rmf_name)
    if os.path.isfile(xcm_dir+"/fake_spectra/"+arf_name) == False:
        os.system("cp "+arf_path+" "+xcm_dir+"/fake_spectra/"+arf_name)

    for i in tqdm(range(trial_number), desc="Generating fake spectra"):
        cmd = "cd "+xcm_dir+"\n"
        if os.path.isfile(xcm_dir+"/fake_spectra/simulated_spectrum_no"+str(i)+".fak") == False:
            cmd += "xspec <<EOF\n"
            cmd += "@"+str(xcm_name)+"\n"
            cmd += "fakeit\n"
            cmd += "y\n"
            cmd += "\n"
            cmd += "simulated_spectrum_no"+str(i)+".fak\n"
            cmd += str(exposure)+", 1.0, "+str(exposure)+"\n"
            cmd += "quit\n"
            cmd += "yes\n"
            cmd += "EOF\n"
            cmd += "mv "+str(xcm_dir)+"/simulated_spectrum_no"+str(i)+".fak "+str(xcm_dir)+"/fake_spectra/\n"
            cmd += "mv "+str(xcm_dir)+"/simulated_spectrum_no"+str(i)+"_bkg.fak "+str(xcm_dir)+"/fake_spectra/\n"

        if os.path.isfile(xcm_dir+"/fake_spectra/simulated_spectrum_no"+str(i)+".qdp") == False:
            cmd += "cd "+xcm_dir+"/fake_spectra/\n"
            cmd += "xspec <<EOF\n"
            cmd += "data simulated_spectrum_no"+str(i)+".fak\n"
            cmd += "back simulated_spectrum_no"+str(i)+"_bkg.fak\n"
            cmd += "ignore **-"+str(continuum_energy_low)+" "+str(continuum_energy_upp)+"-**\n"
            cmd += "iplot\n"
            cmd += "wd simulated_spectrum_no"+str(i)+".qdp\n"
            cmd += "quit\n"
            cmd += "quit\n"
            cmd += "yes\n"
            cmd += "EOF\n"
    
        tqdm.write(cmd)
        cmd_result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if cmd_result.stdout:
            tqdm.write(cmd_result.stdout)
        
        if cmd_result.stderr:
            tqdm.write(cmd_result.stderr)
    
    #### Matched-filtering spectra ####
    names = ['e','de','rate','rate_err']
    df = pd.read_table(xcm_dir+"/fake_spectra/simulated_spectrum_no0.qdp",skiprows=3,names=names, delimiter=' ')
    simulated_energy = df.e
    simulated_rate_all = [[] for i in range(len(simulated_energy))]
    simulated_rate_all_unfiltered = [[] for i in range(len(simulated_energy))]
    simulation_rate_median = []
    simulation_rate_sigma_upp = []
    simulation_rate_sigma_low = []
    simulation_rate_3sigma_upp = []
    simulation_rate_3sigma_low = []

    for i in tqdm(range(trial_number), desc="Matched-filtering fake spectra"):  
        names = ['e','de','rate','rate_err']
        df = pd.read_table(xcm_dir+"/fake_spectra/simulated_spectrum_no"+str(i)+".qdp",skiprows=3,names=names, delimiter=' ')
        simulated_energy = df.e
        simulated_energy_err = df.de
        simulated_rate = df.rate
        matched_filtered_curve = matched_filter(simulated_energy, simulated_energy_err, simulated_rate, a, b, c, d)

        for j in range(len(simulated_energy)):
            simulated_rate_all[j].append(matched_filtered_curve[j])
            simulated_rate_all_unfiltered[j].append(simulated_rate[j])
    

    for i in range(len(simulated_energy)):
        simulation_rate_median.append(statistics.median(simulated_rate_all[i]))
        simulation_rate_sigma_upp.append(statistics.median(simulated_rate_all[i])+statistics.stdev(simulated_rate_all[i]))
        simulation_rate_sigma_low.append(statistics.median(simulated_rate_all[i])-statistics.stdev(simulated_rate_all[i]))
        simulation_rate_3sigma_upp.append(statistics.median(simulated_rate_all[i])+3.0*statistics.stdev(simulated_rate_all[i]))
        simulation_rate_3sigma_low.append(statistics.median(simulated_rate_all[i])-3.0*statistics.stdev(simulated_rate_all[i]))

    #Plotting 
    fig=plt.figure(figsize=(12.5, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=(3.0, 1.2))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    axs = []
    for i in range(2):
        axs.append(fig.add_subplot(gs[i, 0])) 
    label_fontsize = 20.5
    legend_size = 17
    markersize = 6
    linewidth_data = 2.0
    linewidth_model = 2.5
    str_fontsize = 23
    color2 = (32/255, 51/255, 136/255)
    for ax in axs:
        ax.tick_params(direction = "in", right = True, labelright = False, top = True, labeltop = False, labelleft=True, labelbottom = True, labelsize=label_fontsize)
        ax.set_xlim(simulated_energy.iloc[5], simulated_energy.iloc[-5])
        ax.set_xscale('log')
        ax.set_xticks(extract_integers(simulated_energy.iloc[5], simulated_energy.iloc[-5]))
        ax.set_xticklabels(extract_integers_as_strings(simulated_energy.iloc[5], simulated_energy.iloc[-5]))
        ax.yaxis.set_label_coords(-0.080, 0.50)
    
    axs[0].set_yscale('log')
    axs[0].text(0.015, 0.93, "(a)", fontsize=str_fontsize+5, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[0].transAxes, fontweight = "bold", ha="left")
    axs[0].tick_params(direction = "in", right = True, labelright = False, top = True, labeltop = False, labelleft=True, labelbottom = False, labelsize=label_fontsize)
    axs[0].set_ylim(min(spectrum_rates), max(spectrum_rates))
    # axs[0].set_yticks([1, 1.5, 2])
    # axs[0].set_yticklabels(["1.0", "1.5", "2.0"])
    axs[0].plot(spectrum_energy, matched_filter(spectrum_energy, spectrum_energy_err, spectrum_rates, a, b, c, d), color="black", linewidth=linewidth_model+0.9, marker="o", markersize=markersize, markerfacecolor="white", markeredgewidth=linewidth_model-0.2, zorder=3, label="Matched filtered spectrum")
    axs[0].plot(spectrum_energy, continuum_model_all, color="dimgrey", linewidth=linewidth_model-0.8, linestyle="dashdot", zorder=2, label="Best-fitting continuum model")
    axs[0].set_ylabel(r"Counts ($\mathrm{s^{-1}}$ $\mathrm{keV^{-1}}$)", fontsize = label_fontsize+2.5)
    axs[0].plot(simulated_energy, simulation_rate_sigma_low, color="darkgreen", linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2, label=r"$\pm 1 \sigma$ limit")
    axs[0].plot(simulated_energy, simulation_rate_sigma_upp, color="darkgreen", linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2)
    axs[0].plot(simulated_energy, simulation_rate_3sigma_low, color=color2, linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2, label=r"$\pm 3 \sigma$ limit")
    axs[0].plot(simulated_energy, simulation_rate_3sigma_upp, color=color2, linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2)


    linear_interp_sigma_low = interp1d(simulated_energy, simulation_rate_sigma_low, kind='linear')
    linear_interp_sigma_upp = interp1d(simulated_energy, simulation_rate_sigma_upp, kind='linear')
    linear_interp_3sigma_low = interp1d(simulated_energy, simulation_rate_3sigma_low, kind='linear')
    linear_interp_3sigma_upp = interp1d(simulated_energy, simulation_rate_3sigma_upp, kind='linear')

    xd = np.linspace(simulated_energy.iloc[2], simulated_energy.iloc[-2], 3000)
    axs[0].fill_between(xd, linear_interp_sigma_low(xd), linear_interp_sigma_upp(xd), fc="darkgreen", alpha=0.10)
    axs[0].fill_between(xd, linear_interp_sigma_upp(xd), linear_interp_3sigma_upp(xd), fc=color2, alpha=0.10)
    axs[0].fill_between(xd, linear_interp_3sigma_low(xd), linear_interp_sigma_low(xd), fc=color2, alpha=0.10)
    axs[0].minorticks_on()
    axs[0].grid(False, which="minor")
    axs[0].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='x', zorder=1)
    axs[0].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='y', zorder=1)
    axs[0].legend(loc = 'lower left', fontsize=16.5, facecolor=None, edgecolor=None, framealpha=0.0)

    spectrum_rates_filtered = matched_filter(spectrum_energy, spectrum_energy_err, spectrum_rates, a, b, c, d)
    spectrum_rates_ratio = []
    for i in range(len(spectrum_energy)):
        spectrum_rates_ratio.append(spectrum_rates_filtered[i]/continuum_model_all[i])


    simulation_rate_median_ratio = []
    simulation_rate_sigma_upp_ratio = []
    simulation_rate_sigma_low_ratio = []
    simulation_rate_3sigma_upp_ratio = []
    simulation_rate_3sigma_low_ratio = []

    for i in range(len(simulated_energy)):
        simulation_rate_median_ratio.append(simulation_rate_median[i]/continuum_model_all[i])
        simulation_rate_sigma_upp_ratio.append(simulation_rate_sigma_upp[i]/continuum_model_all[i])
        simulation_rate_sigma_low_ratio.append(simulation_rate_sigma_low[i]/continuum_model_all[i])
        simulation_rate_3sigma_upp_ratio.append(simulation_rate_3sigma_upp[i]/continuum_model_all[i])
        simulation_rate_3sigma_low_ratio.append(simulation_rate_3sigma_low[i]/continuum_model_all[i])

    linear_interp_sigma_low_ratio = interp1d(simulated_energy, simulation_rate_sigma_low_ratio, kind='linear')
    linear_interp_sigma_upp_ratio = interp1d(simulated_energy, simulation_rate_sigma_upp_ratio, kind='linear')
    linear_interp_3sigma_low_ratio = interp1d(simulated_energy, simulation_rate_3sigma_low_ratio, kind='linear')
    linear_interp_3sigma_upp_ratio = interp1d(simulated_energy, simulation_rate_3sigma_upp_ratio, kind='linear')



    axs[1].plot(spectrum_energy, spectrum_rates_ratio, color="black", linewidth=linewidth_model+0.9, marker="o", markersize=markersize, markerfacecolor="white", markeredgewidth=linewidth_model-0.2, zorder=3, label="Matched filtered spectrum")
    axs[1].plot(simulated_energy, simulation_rate_sigma_low_ratio, color="darkgreen", linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2, label=r"$\pm \sigma$ limit")
    axs[1].plot(simulated_energy, simulation_rate_sigma_upp_ratio, color="darkgreen", linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2)
    axs[1].plot(simulated_energy, simulation_rate_3sigma_low_ratio, color=color2, linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2, label=r"$\pm 3 \sigma$ limit")
    axs[1].plot(simulated_energy, simulation_rate_3sigma_upp_ratio, color=color2, linewidth=linewidth_model-0.6, linestyle="dashed", zorder=2)
    axs[1].fill_between(xd, linear_interp_sigma_low_ratio(xd), linear_interp_sigma_upp_ratio(xd), fc="darkgreen", alpha=0.10)
    axs[1].fill_between(xd, linear_interp_sigma_upp_ratio(xd), linear_interp_3sigma_upp_ratio(xd), fc=color2, alpha=0.10)
    axs[1].fill_between(xd, linear_interp_3sigma_low_ratio(xd), linear_interp_sigma_low_ratio(xd), fc=color2, alpha=0.10)
    axs[1].set_ylabel(r"Ratio", fontsize = label_fontsize+2.5)
    axs[1].set_xlabel(r"Energy ($\mathrm{keV}$)", fontsize = label_fontsize+2.5)
    axs[1].set_ylim(0.7, 1.8)
    axs[1].plot([-10000,100000], [0,0], marker="+", color="black", markersize=0, markeredgecolor="silver", markerfacecolor="dimgray", markeredgewidth=2.0, linewidth=linewidth_data-0.4, linestyle="dashed")
    axs[1].text(0.015, 0.8, "(b)", fontsize=str_fontsize+5, bbox=dict(boxstyle="round",facecolor="None",edgecolor="None",linewidth=1.3,linestyle="-"), transform=axs[1].transAxes, fontweight = "bold", ha="left")
    # axs[1].set_yticks([0.5, 1.0, 1.5])
    axs[1].minorticks_on()
    axs[1].grid(False, which="minor")
    axs[1].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='x', zorder=1)
    axs[1].grid(True, which="major", color="lightgray", linewidth=1.0, linestyle="dashed", axis='y', zorder=1)
    
    
    plt.subplots_adjust(hspace=0, left=0.1, right = 1-0.1, top=0.995, bottom=0.082)
    plt.savefig("Matched_filtered_spectrum.pdf", format="pdf", dpi=1000)



def main():
    ## Path list ## 
    nicer_rmf_path = "/Desktop/SigmaGem/ObsID/observed_spectrum.rmf"
    nicer_arf_path = "/Desktop/SigmaGem/ObsID/observed_spectrum.arf"
    nicer_spectrum_txt_path = "/Desktop/SigmaGem/ObsID/observed_spectrum.txt"
    nicer_continuum_txt_path = "/Desktop/SigmaGem/ObsID/continuum.txt"
    nicer_continuum_xcm_path = "/Desktop/SigmaGem/ObsID/continuum.xcm"

    ## Exposure and Trial Number ## 
    exposure = 1113 #sec
    trial_number = 10000

    ## Energy band of Continuum ##
    continuum_energy_low = "5.0" #keV
    continuum_energy_upp = "8.0" # keV



    nicer_fwhm_function = calc_fwhm(nicer_rmf_path)
    plot_spectrum_with_MC(nicer_rmf_path, nicer_arf_path, nicer_spectrum_txt_path, nicer_continuum_txt_path, nicer_continuum_xcm_path, exposure, trial_number, continuum_energy_low, continuum_energy_upp, nicer_fwhm_function[0], nicer_fwhm_function[1], nicer_fwhm_function[2], nicer_fwhm_function[3])


if __name__ == "__main__":
    main()
