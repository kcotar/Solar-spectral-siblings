import imp
from astropy.table import Table
from socket import gethostname
from scipy.signal import correlate
from lmfit.models import LinearModel, GaussianModel, VoigtModel, LorentzianModel

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# PC hostname
pc_name = gethostname()

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr52_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
    galah_data_input = '/home/klemen/GALAH_data/'
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
else:
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *

# some settings
min_wvl = list([4730, 5710, 6490, 7710])
max_wvl = list([4890, 5860, 6720, 7870])

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix = '_ext0'
solar_input_dir = galah_data_input+'Solar_data/'
solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 250)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# linelist mask
idx_lines_mask = solar_wvl < 0.
for line in galah_linelist:
    idx_lines_mask[np.logical_and(solar_wvl >= line['line_start']-0.1, solar_wvl <= line['line_end']+0.1)] = True
print 'Linelist mask pixels', np.sum(idx_lines_mask)

# find Solar parameters
teff_solar = np.nanmedian(galah_param[idx_rows]['teff_guess'])
teff_solar_std = np.nanstd(galah_param[idx_rows]['teff_guess'])
logg_solar = np.nanmedian(galah_param[idx_rows]['logg_guess'])
logg_solar_std = np.nanstd(galah_param[idx_rows]['logg_guess'])
feh_solar = np.nanmedian(galah_param[idx_rows]['feh_guess'])
print 'Solar parameters:', teff_solar, '+/-', teff_solar_std, ',  ', logg_solar, '+/-', logg_solar_std, ',  ', feh_solar

# Search for objects with similar physical properties
idx_solar_like = np.logical_and(np.abs(galah_param['teff_guess'] - teff_solar) < teff_solar_std,
                                np.abs(galah_param['logg_guess'] - logg_solar) < logg_solar_std)
idx_solar_like = np.logical_and(idx_solar_like, galah_param['red_flag'] == 0)
idx_solar_like = np.logical_and(idx_solar_like, galah_param['snr_c2_iraf'] > 25)
n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

solar_like_sobjects = galah_param['sobject_id'][idx_solar_like]
sim_results = Table(names=('sobject_id', 'similarity'),
                    dtype=('int64', 'float64'))
for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
    # get spectra of all bands for observed objects
    flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=4, individual=False)
    # compute per band spectra similarity
    spectrum_sim = 0.
    spectrum_sim_valpx = 0.
    for i_c in range(4):
        # define subset of spectra to be compared to reference solar spectrum
        idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
        flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=3)
        # filter out possible strange flux values that will cause anomolues distance estimation
        idx_bad = np.logical_or(flux_b_res > 1.2, flux_b_res < 0.01)
        flux_b_res[idx_bad] = np.nan
        # compute similarity/distance estimator
        diff = (solar_flx[idx_ref] - flux_b_res) ** 2
        # mask difference by elements absorption lines
        diff = diff[idx_lines_mask[idx_ref]]
        #
        spectrum_sim += np.nansum(diff)
        spectrum_sim_valpx += np.sum(np.isfinite(diff))
    spectrum_sim = np.sqrt(spectrum_sim / spectrum_sim_valpx)
    sim_results.add_row([s_obj, spectrum_sim])

print sim_results
print ''
sobj_id_like = sim_results[np.argsort(sim_results['similarity'])[:50]]['sobject_id']
print ','.join([str(s) for s in sobj_id_like])

print ''
sobj_id_dislike = sim_results[np.argsort(sim_results['similarity'])[-50:]]['sobject_id']
print ','.join([str(s) for s in sobj_id_dislike])


