import imp
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import varconvolve as varcon
from scipy import mgrid

from solar_siblings_functions import plot_spectra_with_difference

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize


def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)


def get_spectra_subset(data, range):
    idx = np.logical_and(data[:, 0] < range[1], data[:, 0] > range[0])
    return data[idx, 1], data[idx, 0]


data_dir = '/home/klemen/data4_mount/'
solar_data_dir = data_dir+'Solar_data_dr53/'
galah_linelist = Table.read(data_dir+'GALAH_Cannon_linelist_newer.csv')

# reference spectra
solar_ref = pd.read_csv(solar_data_dir + 'solar_spectra.txt', header=None, delimiter=' ', na_values='nan').values
solar_ref_conv = pd.read_csv(solar_data_dir + 'solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values

# Galah spectrum
solar_galah = pd.read_csv(solar_data_dir + 'twilight_spectrum_galah_ext0_order10.txt', header=None, delimiter=' ', na_values='nan').values

# # convolve solar spectrum - different, modified and visually modified values
# kernel_widths = np.linspace(0.07, 0.155, solar_ref.shape[0])
# solar_ref[:, 1] = varcon.varconvolve(solar_ref[:, 0], solar_ref[:, 1], kernel, kernel_widths)
# # export convolved spectra
# txt = open(solar_data_dir + 'solar_spectra_conv.txt', 'w')
# for i_l in range(solar_ref.shape[0]):
#     txt.write(str(solar_ref[i_l, 0])+' '+str(solar_ref[i_l, 1])+'\n')
# txt.close()

min_wvl = list([4705, 5640, 6475, 7680])
max_wvl = list([4915, 5885, 6750, 7900])

sl = 2.
sh = 3.
st = 17
ord = 1
flux_offset = [0.04, 0.03, 0.03, 0.03]
flux_offset_amp = [0.06, 0.04, 0.1, 0.01]

for i_b in range(4):
    wvl_range = (min_wvl[i_b], max_wvl[i_b])
    flx_galah, wvl_galah = get_spectra_subset(solar_galah, wvl_range)
    flx_ref, wvl_ref = get_spectra_subset(solar_ref_conv, wvl_range)

    flx_galah = spectra_normalize(wvl_galah-np.mean(wvl_galah), flx_galah, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
    flx_ref = spectra_normalize(wvl_ref-np.mean(wvl_ref), flx_ref, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')

    plot_spectra_with_difference(flx_ref, flx_galah, wvl_ref, x_range=wvl_range, linelist=galah_linelist)
