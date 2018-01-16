import imp
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import varconvolve as varcon
from scipy import mgrid
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

data_input_dir = '/home/klemen/data4_mount/Solar_data_dr53/'
galah_linelist = Table.read('/home/klemen/data4_mount/GALAH_Cannon_linelist_newer.csv')

solar_ref = pd.read_csv(data_input_dir + 'solar_spectra.txt', header=None, delimiter=' ', na_values='nan').values
suffix = '_ext0_2'
solar_g1 = pd.read_csv(data_input_dir + 'b1_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(data_input_dir + 'b2_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(data_input_dir + 'b3_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(data_input_dir + 'b4_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values

solar_ref_conv = pd.read_csv(data_input_dir + 'solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values
y_conv1 = solar_ref_conv[:, 1]

solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))

x = solar_ref[:, 0]
y = solar_ref[:, 1]

min_wvl = list([4730, 5660, 6490, 7700])
max_wvl = list([4890, 5860, 6725, 7870])

# resolution degradation computation
R_solar_ref = 1e10
R_galah = 28000.
wvl_start = x[0]
wvl_end = x[-1]

# theoretical values for kernel convolution
s_ref_start = wvl_start/R_solar_ref
s_ref_end = wvl_end/R_solar_ref
s_galah_start = wvl_start/R_galah
s_galah_end = wvl_end/R_galah
s_conv_start = np.sqrt(s_galah_start**2 - s_ref_start**2)
s_conv_end = np.sqrt(s_galah_end**2 - s_ref_end**2)
print s_conv_start, s_conv_end

# convolve solar spectrum - different, modified and visually modified values
kernel_widths = np.linspace(0.07, 0.17, len(y))
y_conv1 = varcon.varconvolve(x, y, kernel, kernel_widths)
# export convolved spectra
txt = open('solar_spectra_conv.txt', 'w')
for i_l in range(len(x)):
    txt.write(str(x[i_l])+' '+str(y_conv1[i_l])+'\n')
txt.close()

# # plot everything
# plt.plot(x, y_conv1, c='black', lw=1.5)
# plt.plot(solar_wvl, solar_flx, c='blue', lw=1.5)
# # kernel_widths = np.linspace(0.07, 0.17, len(y))
# # y_conv1 = varcon.varconvolve(x, y, kernel, kernel_widths)
# # plt.plot(x, y_conv1, c='green', lw=1.5)
# plt.show()
# plt.close()

raise SystemExit

sl = 2.
sh = 3.
st = 17
ord = 5
flux_offset = [0.04, 0.03, 0.03, 0.03]
flux_offset_amp = [0.06, 0.04, 0.1, 0.01]
# flux_offset_amp = [0.00, 0.00, 0.0, 0.00]
for i_b in range(4):
    d_wvl = 0
    idx_ref_sub = np.logical_and(x > min_wvl[i_b]-d_wvl, x < max_wvl[i_b]+d_wvl)
    y_conv1_norm = spectra_normalize(x[idx_ref_sub], y_conv1[idx_ref_sub], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
    idx_sol_sub = np.logical_and(solar_wvl > min_wvl[i_b]-d_wvl, solar_wvl < max_wvl[i_b]+d_wvl)
    flux_ref_solar_b = y_conv1[np.in1d(x, solar_wvl[idx_sol_sub])]

    y_off_perwvl = (1. - 1. * np.arange(np.sum(idx_sol_sub)) / np.sum(idx_sol_sub)) * flux_offset_amp[i_b] + flux_offset[i_b]
    solar_flux_norm_offset = spectra_normalize(solar_wvl[idx_sol_sub], solar_flx[idx_sol_sub] - y_off_perwvl, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')

    plt.figure(1, figsize=(12, 7))
    # [left, bottom, width, height]
    axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
    axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])

    axSpectra.plot(x[idx_ref_sub], y_conv1_norm, c='black', lw=2)
    axSpectra.plot(solar_wvl[idx_sol_sub], solar_flux_norm_offset, c='blue', lw=2)

    # # export results
    # txt_out = 'b'+str(i_b+1)+'_solar_galah' + suffix + '_offset.txt'
    # txt = open(txt_out, 'w')
    # for i_l in range(len(solar_flux_norm_offset)):
    #     txt.write(str(solar_wvl[idx_sol_sub][i_l])+' '+str(solar_flux_norm_offset[i_l])+'\n')
    # txt.close()

    d_abs_wvl = 0.0
    for line in galah_linelist:
        if line['line_centre'] <solar_wvl[idx_sol_sub][-1] and line['line_centre'] > solar_wvl[idx_sol_sub][0]:
            axSpectra.axvspan(line['line_start'] - d_abs_wvl, line['line_end'] + d_abs_wvl, lw=0, color='black', alpha=0.2)
            axDiff.axvspan(line['line_start'] - d_abs_wvl, line['line_end'] + d_abs_wvl, lw=0, color='black', alpha=0.2)

    axDiff.axhline(y=0, c='black', lw=1)
    axDiff.plot(solar_wvl[idx_sol_sub], y_conv1_norm - solar_flux_norm_offset, c='blue', lw=1.5)

    # print '-----', i_b, '-----'
    # idx_sum = y_conv1_norm < 0.95
    # dy_offset = 0.04
    # for dy_amp in np.arange(0, 0.2, 0.005):
    #     y_off_perwvl = (1. - 1. * np.arange(np.sum(idx_sol_sub)) / np.sum(idx_sol_sub)) * dy_amp + dy_offset
    #     solar_flux_norm_offset = spectra_normalize(solar_wvl[idx_sol_sub], solar_flx[idx_sol_sub]-y_off_perwvl, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
    #     print dy_amp, np.nansum(np.abs(y_conv1_norm - solar_flux_norm_offset)[idx_sum])
    #     axSpectra.plot(solar_wvl[idx_sol_sub], solar_flux_norm_offset, c='red', lw=1.)
    #     axDiff.plot(solar_wvl[idx_sol_sub], y_conv1_norm - solar_flux_norm_offset, c='red', lw=1.)

    axSpectra.set(ylim=(0.3, 1.1), xlim=(min_wvl[i_b]-d_wvl, max_wvl[i_b]+d_wvl))
    axDiff.set(ylim=(-0.05, 0.05), xlim=(min_wvl[i_b]-d_wvl, max_wvl[i_b]+d_wvl))
    plt.show()
    plt.close()
#
# plt.plot(x, spectra_normalize(x, y_conv1, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly'), c='red', lw=2)
# plt.plot(solar_g1[:, 0], spectra_normalize(solar_g1[:, 0], solar_g1[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly'), c='blue', lw=1)
# plt.plot(solar_g2[:, 0], spectra_normalize(solar_g2[:, 0], solar_g2[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly'), c='blue', lw=1)
# plt.plot(solar_g3[:, 0], spectra_normalize(solar_g3[:, 0], solar_g3[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly'), c='blue', lw=1)
# plt.plot(solar_g4[:, 0], spectra_normalize(solar_g4[:, 0], solar_g4[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly'), c='blue', lw=1)




# plt.plot(solar_wvl, solar_flx/flx_fit, c='green', lw=2)
# plt.ylim(0.3, 1.1)
# plt.show()
# plt.close()

# plt.plot(x, y, c='black', lw=2)
# plt.plot(x, y_conv1, c='red', lw=1)
# plt.plot(x, y_fit, c='red', lw=3)
# plt.scatter(x[~idx_fit], y_conv1[~idx_fit], s=10, lw=1, c='black', marker='x')

# plt.plot(solar_wvl, solar_flx, c='green', lw=1)
# plt.plot(solar_wvl, flx_fit, c='green', lw=3)
# plt.scatter(solar_wvl[~idx_flx_fit], solar_flx[~idx_flx_fit], s=10, lw=1, c='black', marker='x')
# plt.show()
# plt.close()
