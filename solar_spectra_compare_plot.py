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

data_input_dir = '/home/klemen/GALAH_data/Solar_data/'

solar_ref = pd.read_csv(data_input_dir + 'solar_spectra.txt', header=None, delimiter=' ', na_values='nan').values
suffix = '_ext0'
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
idx_use_ref = x > 4700
x = x[idx_use_ref]
y = y[idx_use_ref]
y_conv1 = y_conv1[idx_use_ref]

# resolution degradation computation
R_solar_ref = 1e10
R_galah = 28000.
wvl_start = x[0]
wvl_end = x[-1]

# s_ref_start = wvl_start/R_solar_ref
# s_ref_end = wvl_end/R_solar_ref
# s_galah_start = wvl_start/R_galah
# s_galah_end = wvl_end/R_galah
# s_conv_start = np.sqrt(s_galah_start**2 - s_ref_start**2)
# s_conv_end = np.sqrt(s_galah_end**2 - s_ref_end**2)
# print s_conv_start, s_conv_end

# convolve solar spectrum
# kernel_widths = np.linspace(s_conv_start, s_conv_end, len(y))
# y_conv1 = varcon.varconvolve(x, y, kernel, kernel_widths/2.)
# export convolved spectra
# txt = open('solar_spectra_conv.txt', 'w')
# for i_l in range(len(x)):
#     txt.write(str(x[i_l])+' '+str(y_conv1[i_l])+'\n')
# txt.close()

# idx_g1 = np.isfinite(solar_g1[:, 1])
# g1_y = solar_g1[:, 1][idx_g1]
# g1_x = solar_g1[:, 0][idx_g1]
# g1_y[g1_y > 0.95] = 0.
# y[y > 0.95] = 0.
#
# idx_b1 = np.logical_and(x>=np.min(g1_x), x<=np.max(g1_x))
# corr_res = np.correlate(y[idx_b1], g1_y, mode='same')
# plt.plot(corr_res)
# plt.show()
# plt.close()
# print len(g1_y), len(y[idx_b1])
# print np.argmax(corr_res), np.argmax(corr_res) - len(corr_res)/2.


# plot everything
# plt.plot(x, y_conv1, c='red', lw=2)
# y_fit, idx_fit = spectra_normalize(x, y_conv1, steps=9, sigma_low=1.5, sigma_high=3., order=7, func='spline', return_fit=True, return_idx=True)

# idx_fit_for_solar = np.ndarray(len(solar_wvl), dtype=np.bool)
# idx_fit_for_solar.fill(False)
# for wvl_fit_use in x[idx_fit]:
#     idx_fit_for_solar[np.abs(solar_wvl - wvl_fit_use)<1e-5] = True
# print np.sum(idx_fit_for_solar)

# plt.plot(x, y_conv2, c='green', lw=2)
# plt.plot(x, y_conv3, c='orange', lw=2)
# plt.plot(x, y_conv4, c='yellow', lw=2)

# plt.plot(solar_g1[:, 0], solar_g1[:, 1], c='blue', lw=2)
# plt.plot(solar_g2[:, 0], solar_g2[:, 1], c='blue', lw=2)
# plt.plot(solar_g3[:, 0], solar_g3[:, 1], c='blue', lw=2)
# plt.plot(solar_g4[:, 0], solar_g4[:, 1], c='blue', lw=2)

# plt.plot(x, y_conv1, c='red', lw=1)
# plt.plot(x, spectra_normalize(x, y_conv1, steps=7, sigma_low=3., sigma_high=5, order=9, func='spline', return_fit=True), c='black', lw=2)

plt.plot(solar_g1[:, 0], solar_g1[:, 1], c='blue', lw=1)
plt.plot(solar_g1[:, 0], spectra_normalize(solar_g1[:, 0], solar_g1[:, 1], steps=9, sigma_low=2., sigma_high=5., order=7, func='spline', return_fit=True), c='black', lw=2)

plt.plot(solar_g2[:, 0], solar_g2[:, 1], c='blue', lw=1)
plt.plot(solar_g2[:, 0], spectra_normalize(solar_g2[:, 0], solar_g2[:, 1], steps=9, sigma_low=2., sigma_high=5., order=7, func='spline', return_fit=True), c='black', lw=2)

plt.plot(solar_g3[:, 0], solar_g3[:, 1], c='blue', lw=1)
plt.plot(solar_g3[:, 0], spectra_normalize(solar_g3[:, 0], solar_g3[:, 1], steps=9, sigma_low=2., sigma_high=5., order=7, func='spline', return_fit=True), c='black', lw=2)

# plt.plot(solar_g4[:, 0], spectra_normalize(solar_g4[:, 0], solar_g4[:, 1], steps=9, sigma_low=1.5, sigma_high=3., order=7, func='spline'), c='blue', lw=2)
#
# flx_fit, idx_flx_fit = spectra_normalize(solar_wvl, solar_flx, steps=11, sigma_low=1.5, sigma_high=3., order=7,
#                                          func='spline', return_fit=True, return_idx=True)

# plt.plot(solar_wvl, solar_flx/flx_fit, c='green', lw=2)
plt.ylim(0.3, 1.1)
plt.show()
plt.close()

# plt.plot(x, y, c='black', lw=2)
# plt.plot(x, y_conv1, c='red', lw=1)
# plt.plot(x, y_fit, c='red', lw=3)
# plt.scatter(x[~idx_fit], y_conv1[~idx_fit], s=10, lw=1, c='black', marker='x')

# plt.plot(solar_wvl, solar_flx, c='green', lw=1)
# plt.plot(solar_wvl, flx_fit, c='green', lw=3)
# plt.scatter(solar_wvl[~idx_flx_fit], solar_flx[~idx_flx_fit], s=10, lw=1, c='black', marker='x')
# plt.show()
# plt.close()
