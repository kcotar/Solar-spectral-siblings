import imp
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import varconvolve as varcon
from scipy import mgrid

from solar_siblings_functions import plot_spectra_with_difference

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize, move_to_dir


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
twilight_spectrum_file = 'twilight_spectrum_galah_ext0_dateall.txt'
solar_galah = pd.read_csv(solar_data_dir + twilight_spectrum_file, header=None, delimiter=' ', na_values='nan').values

# # convolve solar spectrum - different, modified and visually modified values
# kernel_widths = np.linspace(0.07, 0.155, solar_ref.shape[0])
# solar_ref[:, 1] = varcon.varconvolve(solar_ref[:, 0], solar_ref[:, 1], kernel, kernel_widths)
# # export convolved spectra
# txt = open(solar_data_dir + 'solar_spectra_conv.txt', 'w')
# for i_l in range(solar_ref.shape[0]):
#     txt.write(str(solar_ref[i_l, 0])+' '+str(solar_ref[i_l, 1])+'\n')
# txt.close()

min_wvl = list([4705, 5640, 6475, 7680])
max_wvl = list([4905, 5885, 6750, 7900])

sl = 2.
sh = 3.
st = 15
ord = 1

min_wvl_offset = np.array([4715, 5665, 6485, 7700])
max_wvl_offset = np.array([4890, 5865, 6725, 7875])

perform_analysis = False
final_output = True


def get_spectrum_with_offset(flx, wvl, offset, offset_amp, offset_wvl_range):
    d_wvl = offset_wvl_range[1] - offset_wvl_range[0]
    y_off_perwvl = (1. - (wvl-offset_wvl_range[0])/d_wvl) * offset_amp + offset
    return spectra_normalize(wvl - np.mean(wvl), flx - y_off_perwvl, steps=st, sigma_low=sl,
                             sigma_high=sh, order=ord, func='poly')


move_to_dir('Twilight_offset_determine')

if perform_analysis:
    for i_b in range(4):
        wvl_range = (min_wvl[i_b], max_wvl[i_b])
        flux_offset_wvl_range = (min_wvl_offset[i_b], max_wvl_offset[i_b])
        flx_galah, wvl_galah = get_spectra_subset(solar_galah, wvl_range)
        flx_ref, wvl_ref = get_spectra_subset(solar_ref_conv, wvl_range)

        flx_galah = spectra_normalize(wvl_galah-np.mean(wvl_galah), flx_galah, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
        flx_ref = spectra_normalize(wvl_ref-np.mean(wvl_ref), flx_ref, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')

        results = []
        idx_sim = np.logical_and(wvl_galah >= min_wvl_offset[i_b], wvl_galah <= max_wvl_offset[i_b])
        for flux_offset_amp in np.arange(0, 0.11, 0.005):
            for flux_offset in np.arange(0, 0.11, 0.005):
                flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, flux_offset, flux_offset_amp, flux_offset_wvl_range)
                flx_diff = flx_ref-flx_galah_new
                flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_diff,
                                                 steps=st, sigma_low=2., sigma_high=2., order=12, func='poly', return_fit=True)
                plot_sim = np.nansum(np.abs(flx_diff - flx_diff_fit)[idx_sim])
                results.append([flux_offset, flux_offset_amp, min_wvl_offset[i_b], max_wvl_offset[i_b], plot_sim])
        results = np.vstack(results)
        print 'Best for band', i_b
        idx_best = np.argsort(results[:, 4])
        results_best = results[idx_best[:10], :]
        print results_best

        # plot them
        for res in results_best:
            plot_p = str('b{:.0f}_o{:0.3f}_a{:0.3f}.png'.format(i_b, res[0], res[1]))
            flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, res[0], res[1], (res[2], res[3]))
            flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_ref - flx_galah_new,
                                             steps=st, sigma_low=2., sigma_high=2., order=12, func='poly', return_fit=True)
            plot_spectra_with_difference(flx_ref, flx_galah_new, wvl_ref, x_range=wvl_range, linelist=galah_linelist,
                                         title=str(res[4]), path=plot_p, diff_func=flx_diff_fit)

if final_output:
    # output new spectra as a csv format for later use
    flux_offsets = [0.025, 0.0, 0.005, 0.005]
    flux_offset_amps = [0.045, 0.035, 0.015, 0.0]

    wvl_galah = solar_galah[:, 0]
    flx_galah = solar_galah[:, 1]
    flx_galah_offset = np.array(flx_galah)

    for i_b in range(4):
        flux_offset_wvl_range = (min_wvl_offset[i_b], max_wvl_offset[i_b])
        idx_use = np.logical_and(wvl_galah >= min_wvl[i_b], wvl_galah <= max_wvl[i_b])
        flx_galah_new = get_spectrum_with_offset(flx_galah[idx_use], wvl_galah[idx_use],
                                                 flux_offsets[i_b], flux_offset_amps[i_b], flux_offset_wvl_range)
        flx_galah_offset[idx_use] = flx_galah_new

    # output to file
    out_file = solar_data_dir + twilight_spectrum_file[:-4]+'_offset.txt'
    txt = open(out_file, 'w')
    for i_l in range(len(wvl_galah)):
        txt.write(str(wvl_galah[i_l]) + ' ' + str(flx_galah_offset[i_l]) + '\n')
    txt.close()
