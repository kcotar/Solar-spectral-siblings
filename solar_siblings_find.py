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
min_wvl = list([4710, 5640, 6475, 7700])
max_wvl = list([4915, 5885, 6750, 7900])

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix = '_ext0'
solar_input_dir = galah_data_input+'Solar_data/'
solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_wvl = [solar_g1[0], solar_g2[0], solar_g3[0], solar_g4[0]]
solar_flx = [solar_g1[1], solar_g2[1], solar_g3[1], solar_g4[1]]

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'

# select ok objects
print 'Reading and filtering other data'
galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 250)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# find Solar parameters
teff_solar = np.nanmedian(galah_param[idx_rows]['teff_guess'])
teff_solar_std = np.nanstd(galah_param[idx_rows]['teff_guess'])
logg_solar = np.nanmedian(galah_param[idx_rows]['logg_guess'])
logg_solar_std = np.nanstd(galah_param[idx_rows]['logg_guess'])
feh_solar = np.nanmedian(galah_param[idx_rows]['feh_guess'])
print 'Solar parameters:', teff_solar, teff_solar_std, logg_solar, logg_solar_std, feh_solar

# Search for objects with similar physical properties
idx_solar_like = np.logical_and(np.abs(galah_param['teff_guess'] - teff_solar) < teff_solar_std,
                                np.abs(galah_param['logg_guess'] - logg_solar) < logg_solar_std)
print 'Solar like by parameters:', np.sum(idx_solar_like)


