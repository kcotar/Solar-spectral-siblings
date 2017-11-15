import imp, os
from astropy.table import Table
from socket import gethostname

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# PC hostname
pc_name = gethostname()

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    galah_data_input = '/home/klemen/GALAH_data/'
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
else:
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *

line_file = 'GALAH_Cannon_linelist.csv'
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'
abund_param_file = 'Cannon3.0.1_Sp_SMEmasks_trainingset.fits'
spectra_file_list = ['galah_dr52_ccd1_4710_4910_wvlstep_0.04_lin_'+data_date+'.pkl',
                     'galah_dr52_ccd2_5640_5880_wvlstep_0.05_lin_'+data_date+'.pkl',
                     'galah_dr52_ccd3_6475_6745_wvlstep_0.06_lin_'+data_date+'.pkl',
                     'galah_dr52_ccd4_7700_7895_wvlstep_0.07_lin_'+data_date+'.pkl']

galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c2_iraf'] > 200)
to_read_row = np.where(idx_rows)[0]
print 'Number of solar spectra:', len(to_read_row)

for i_b in range(len(spectra_file_list)):
    print 'Band:', i_b+1
    spectra_file = spectra_file_list[i_b]
    collection_param = CollectionParameters(spectra_file)
    wvl_values = collection_param.get_wvl_values()
    # read all selected GALAH solar spectra
    solar_spectra = read_pkl_spectra(galah_data_input + spectra_file,
                                     read_rows=to_read_row, read_cols=None)
    # simple median stacking of read spectra
    print ' Stacking result'
    solar_spectra_median = np.nanmedian(solar_spectra, axis=0)
    # save result as txt file
    print ' Saving'
    out_file = galah_data_input+'Solar_data/'
    out_file += 'b'+str(i_b+1)+'_solar_galah.txt'
    txt = open(out_file, 'w')
    for i_l in range(len(wvl_values)):
        txt.write(str(wvl_values[i_l])+' '+str(solar_spectra_median[i_l])+'\n')
    txt.close()

