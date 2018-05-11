import imp
from os import chdir, system
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import varconvolve as varcon
from scipy import mgrid


dr52_dir = '/media/storage/HERMES_REDUCED/dr5.3/'
imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import *

chdir('Asiago')

# RV diagnostics
# rv_res = pd.read_csv('rv_res.txt', sep=' ', header=None).values
# plt.plot(rv_res[:, 0])
# plt.show()
# plt.close()

wvl_range = [4000, 7000]


def trim_data(x, y, r):
    idx = np.logical_and(x >= r[0], x <= r[1])
    return x[idx], y[idx]


def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)


ref_orig = pd.read_csv('/home/klemen/data4_mount/Solar_data_dr53/solar_spectra.txt', sep=' ', header=None).values
x_ref, y_ref = trim_data(ref_orig[:, 0], ref_orig[:, 1], wvl_range)
y_ref1 = varcon.varconvolve(x_ref, y_ref, kernel, np.linspace(0.11, 0.19, len(x_ref)))
# y_ref2 = varcon.varconvolve(x_ref, y_ref, kernel, np.linspace(0.12, 0.22, len(x_ref)))
# y_ref3 = varcon.varconvolve(x_ref, y_ref, kernel, np.linspace(0.13, 0.23, len(x_ref)))
# y_ref4 = varcon.varconvolve(x_ref, y_ref, kernel, np.linspace(0.14, 0.24, len(x_ref)))

# plt.plot(x_ref, y_ref2, c='C2', lw=1)
# plt.plot(x_ref, y_ref3, c='C3', lw=1)
# plt.plot(x_ref, y_ref4, c='C4', lw=1)

# s2, w2 = get_spectra_dr52(str(160524006601258), bands=[1,2,3], root=dr52_dir, extension=4)
# for i_b in range(4):
#     plt.plot(w2[i_b], s2[i_b], c='blue', lw=1)

ref_data = pd.read_csv('/home/klemen/data4_mount/Solar_data_dr53/twilight_spectrum_galah_ext0_dateall.txt', sep=' ', header=None).values

asiago_obs = fits.open('combined_final.0001.fits')
asiago_obs_h = asiago_obs[0].header
flx = asiago_obs[0].data
wvl = asiago_obs_h.get('CRVAL1') + asiago_obs_h.get('CDELT1') * np.arange(len(flx))
rv_shift = -49.5
wvl *= (1. - rv_shift / 299792.458)
wvl, flx = trim_data(wvl, flx, wvl_range)

wvl_step = 200
wvl_beg_multi = np.arange(wvl_range[0], wvl_range[1], wvl_step)
n_plots = len(wvl_beg_multi)
fig, ax = plt.subplots(n_plots, 1, figsize=(8, 10.5))
fig.subplots_adjust(hspace=0.0, wspace=0, left=0.125, right=0.99, top=0.99, bottom=0.01)
for i_p in range(n_plots):
    wvl_beg = wvl_beg_multi[i_p]
    wvl_end = wvl_beg+wvl_step
    x_ref_sub, y_ref1_sub = trim_data(x_ref, y_ref1, [wvl_beg, wvl_end])
    wvl_sub, flx_sub = trim_data(wvl, flx, [wvl_beg, wvl_end])

    # normalize them both
    y_ref1_sub = spectra_normalize(x_ref_sub - np.mean(x_ref_sub), y_ref1_sub,
                                   steps=11, sigma_low=2., sigma_high=3., order=1, n_min_perc=5.,
                                   return_fit=False, func='poly')
    flx_sub = spectra_normalize(wvl_sub - np.mean(wvl_sub), flx_sub,
                                steps=11, sigma_low=2., sigma_high=3., order=1, n_min_perc=5.,
                                return_fit=False, func='poly')

    ax[i_p].plot(x_ref_sub, y_ref1_sub, c='black', lw=.4)
    ax[i_p].plot(wvl_sub, flx_sub, c='red', lw=.4)
    ax[i_p].set(ylim=(0.4, 1.05), xlim=(wvl_beg, wvl_end))
    ax[i_p].set_ylabel(r''+str(wvl_beg)+' - '+str(wvl_end), rotation=0)
    ax[i_p].yaxis.set_label_coords(-0.08, 0.4)
    ax[i_p].xaxis.set_ticklabels([])
    ax[i_p].yaxis.set_ticklabels([])
    # if i_p == n_plots-1:
    #     ax[i_p].set(xlabel=r'Wavelenght [$\AA$]')
# plt.show()
file_out = 'asiago_ref_comp.eps'
plt.savefig(file_out, format='eps')
plt.close()
system('ps2pdf -dEPSCrop ' + file_out + '  ' + file_out.split('.')[0]+'.pdf')
