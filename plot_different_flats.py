import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, imp
from glob import glob
from astropy.table import Table
import varconvolve as varcon
from scipy import mgrid

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize, spectra_resample

def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)

os.chdir('/home/klemen/data4_mount/Solar_data_dr53/')
galah_linelist = Table.read('/home/klemen/data4_mount/GALAH_Cannon_linelist_newer.csv')
# telluric_spectrum = np.loadtxt('/home/klemen/data4_mount/telluric_spectra.dat')
# telluric_spectrum = np.loadtxt('/home/klemen/data4_mount/telluric_spectra_conv.dat')
telluric_spectrum_O2 = np.loadtxt('/home/klemen/data4_mount/telluric_O2_conv.dat')
telluric_spectrum_H2O = np.loadtxt('/home/klemen/data4_mount/telluric_H2O_conv.dat')
solar_ref_conv = pd.read_csv('/home/klemen/data4_mount/Solar_data_dr53/solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values

# # convovle telluric spectrum
# kernel_widths = np.linspace(0.07, 0.145, telluric_spectrum_H20.shape[0])
# telluric_spectrum_H20[:, 1] = varcon.varconvolve(telluric_spectrum_H20[:, 0], telluric_spectrum_H20[:, 1], kernel, kernel_widths)
# txt = open('/home/klemen/data4_mount/telluric_H2O_conv.dat', 'w')
# for i_l in range(telluric_spectrum_H20.shape[0]):
#     txt.write(str(telluric_spectrum_H20[i_l, 0])+' '+str(telluric_spectrum_H20[i_l, 1])+'\n')
# txt.close()
# raise SystemExit

# renormalize data to the same
sl = 2.
sh = 3.
st = 15
ord = 1
telluric_spectrum_O2[:, 1] = spectra_normalize(telluric_spectrum_O2[:, 0]-np.mean(telluric_spectrum_O2[:, 0]), telluric_spectrum_O2[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
telluric_spectrum_H2O[:, 1] = spectra_normalize(telluric_spectrum_H2O[:, 0]-np.mean(telluric_spectrum_H2O[:, 0]), telluric_spectrum_H2O[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
solar_ref_conv[:, 1] = spectra_normalize(solar_ref_conv[:, 0]-np.mean(solar_ref_conv[:, 0]), solar_ref_conv[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
solar_ref_conv_res = spectra_resample(solar_ref_conv[:, 1], solar_ref_conv[:, 0], telluric_spectrum_H2O[:, 0], k=1)


for spectrum_file in glob('twilight_spectrum_galah_ext0_date*.txt'):
    spectrum = pd.read_csv(spectrum_file, header=None, delimiter=' ', na_values='nan').values
    spectrum[:, 1] = spectra_normalize(spectrum[:, 0] - np.mean(spectrum[:, 0]), spectrum[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
    plt.plot(spectrum[:, 0], spectrum[:, 1], c='red', lw=1, alpha=0.5)
# add telluric spectrum
plt.plot(solar_ref_conv[:, 0], solar_ref_conv[:, 1], c='black', lw=1.5)
telluric_contr = (1. - telluric_spectrum_H2O[:, 1])*4. + (1. - telluric_spectrum_O2[:, 1])*3.
plt.plot(telluric_spectrum_H2O[:, 0], solar_ref_conv_res - telluric_contr, lw=1.5, c='blue')
for line in galah_linelist:
    plt.axvspan(line['line_start'], line['line_end'], lw=0, color='black', alpha=0.2)
plt.ylim((0, 1.3))
plt.show()
plt.close()

