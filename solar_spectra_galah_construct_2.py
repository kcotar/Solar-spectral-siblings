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
print 'Read reference Solar spectra'
solar_ref = pd.read_csv(galah_data_input + 'Solar_data/solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values
solar_ref_flux_all = solar_ref[:, 1]
solar_ref_wvl_all = solar_ref[:, 0]

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'

# select ok objects
galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 250)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)
galah_param = galah_param[idx_rows]
to_read_row = np.where(idx_rows)[0]
print 'Number of solar spectra:', len(to_read_row)

# output 100 sobject_ids from random solar spectra
# s_ids = galah_param['sobject_id'][list(np.int64(np.random.rand(250)*len(galah_param)))]
# print ','.join([str(s) for s in s_ids])
# raise SystemExit

# galah_param = galah_param[list(np.int64(np.random.rand(500)*len(galah_param)))]

# do the whole procedure for every spectra and every ccd
for ccd in [1, 2, 3, 4]:
    print 'Band', ccd
    idx_solar_sub = np.logical_and(solar_ref_wvl_all >= min_wvl[ccd-1], solar_ref_wvl_all <= max_wvl[ccd-1])
    solar_ref_flux = solar_ref_flux_all[idx_solar_sub]
    solar_ref_wvl = solar_ref_wvl_all[idx_solar_sub]
    wvl_step = solar_ref_wvl[1] - solar_ref_wvl[0]
    galah_solar = np.ndarray((len(galah_param), len(solar_ref_wvl)))

    for i_s_obj in range(len(galah_param)):
        s_obj = galah_param[i_s_obj]['sobject_id']
        # print s_obj
        flux, wvl = get_spectra_dr52(str(s_obj), bands=[ccd], root=dr52_dir, extension=0, individual=False)

        if ccd == 4:
            idx_wvl_use = wvl[0] >= 7680
            wvl[0] = wvl[0][idx_wvl_use]
            flux[0] = flux[0][idx_wvl_use]

        # init fit
        flux[0] = spectra_normalize(wvl[0], flux[0], steps=5, sigma_low=1.5, sigma_high=2.5, order=1, n_min_perc=5.,
                                    return_fit=False, func='poly')

        flux_fit1 = spectra_normalize(wvl[0], flux[0], steps=3, sigma_low=1.8, sigma_high=3., order=3, n_min_perc=5.,
                                      return_fit=True, func='spline', window=5)
        flux_fit2 = spectra_normalize(wvl[0], flux[0], steps=25, sigma_low=1.8, sigma_high=3., order=17, n_min_perc=5.,
                                      return_fit=True, func='poly')
        plt.plot(wvl[0], flux[0], c='black', lw=1)
        plt.plot(wvl[0], flux_fit1, c='red', alpha=0.5, lw=2)
        plt.plot(wvl[0], flux_fit2, c='green', alpha=0.5, lw=2)
        plt.ylim((0.4, 1.2))
        plt.show()
        plt.close()

        # flux_fit3 = spectra_normalize(wvl[0], flux[0], steps=5, sigma_low=2., sigma_high=3., order=3, n_min_perc=5.,
        #                               return_fit=True, func='spline', window=5)
        # flux_fit4 = spectra_normalize(wvl[0], flux[0], steps=7, sigma_low=2., sigma_high=3., order=3, n_min_perc=5.,
        #                               return_fit=True, func='spline', window=5)

        flux[0] = flux[0]/flux_fit2
        plt.plot(wvl[0], flux[0], c='black', lw=1)
        # plt.plot(wvl[0], flux_fit3, c='red', alpha=0.5, lw=2)
        # plt.plot(wvl[0], flux_fit4, c='green', alpha=0.5, lw=2)
        plt.ylim((0.4, 1.2))
        plt.show()
        continue

        flux[0] = spectra_normalize(wvl[0], flux[0], steps=35, sigma_low=1.5, sigma_high=2.8, order=29, n_min_perc=3.,
                                    return_fit=False, func='poly')
        # resample read spectra to the common wvl step found in reference solar spectra
        flux_res = spectra_resample(flux[0], wvl[0], solar_ref_wvl, k=1)
        # store to the array
        # galah_solar[i_s_obj, :] = flux_res

        # correlation and stuff
        # get a valid subset of data
        idx_valid = np.isfinite(flux_res)
        flux_res = flux_res[idx_valid]
        solar_ref_flux_sub = solar_ref_flux[idx_valid]
        wvl_valid = solar_ref_wvl[idx_valid]
        # convert spectra sampling to logspace
        flux_res_log, wvl_valid_log = spectra_logspace(flux_res, wvl_valid)
        solar_ref_flux_sub_log, _ = spectra_logspace(solar_ref_flux_sub, wvl_valid)
        # correlate the two spectra
        min_flux = 0.97
        solar_ref_flux_sub_log[solar_ref_flux_sub_log > min_flux] = 0.
        flux_res_log[flux_res_log > min_flux] = 0.
        corr_res = correlate(solar_ref_flux_sub_log, flux_res_log, mode='same', method='fft')

        # create a correlation subset that will actually be analysed
        corr_w_size = 550
        corr_c_off = np.int64(len(corr_res) / 2.)
        corr_pos_min = corr_c_off - corr_w_size
        corr_pos_max = corr_c_off + corr_w_size
        corr_res_sub = corr_res[corr_pos_min:corr_pos_max]
        corr_res_sub -= np.median(corr_res_sub)
        corr_res_sub_x = np.arange(len(corr_res_sub))
        # analyze correlation function by fitting gaussian/voigt/lorentzian distribution to it
        fit_model = VoigtModel()
        parameters = fit_model.guess(corr_res_sub, x=corr_res_sub_x)
        corr_fit_res = fit_model.fit(corr_res_sub, parameters, x=corr_res_sub_x)
        corr_center = corr_fit_res.params['center'].value

        # determine the actual shift
        idx_no_shift = np.int32(len(corr_res) / 2.)
        idx_center = corr_c_off - corr_w_size + corr_center
        log_shift_px = idx_no_shift - idx_center
        log_shift_wvl = log_shift_px * wvl_step

        # move spectra for the determined shift
        flux_res_log, wvl_valid_log = spectra_logspace(flux_res, wvl_valid)
        wvl_valid_log -= log_shift_wvl
        flux_rv_shifted = np.interp(wvl_valid, wvl_valid_log, flux_res_log)

        # store to the array
        print s_obj, log_shift_px, log_shift_wvl
        galah_solar[i_s_obj, idx_valid] = flux_rv_shifted

        """
        # print corr_res
        plt.plot(corr_res_sub, c='black')
        plt.plot(corr_fit_res.best_fit, c='red')
        plt.axvline(x=corr_center, c='red', lw=1)
        plt.show()
        plt.close()

        # plot some spectra
        plt.plot(flux_res, c='red')
        plt.plot(solar_ref_flux_sub, c='black')
        plt.plot(flux_rv_shifted, c='green')
        plt.show()
        plt.close()
        """

    #     plt.plot(wvl_valid, flux_rv_shifted, c='blue', alpha=0.01, lw=1)
    # plt.show()
    # plt.close()

    # filter out outlying spectra (based on supplied reference solar spectra)
    eucl_match = np.sqrt(np.nansum((galah_solar - solar_ref_flux)**2, axis=1)/np.sum(np.isfinite(galah_solar), axis=1))
    eucl_match_thr = np.nanpercentile(eucl_match, 95.)
    idx_bad_match = eucl_match >= eucl_match_thr
    galah_solar[np.where(idx_bad_match)[0], :] = np.nan
    print 'Removed rows:', np.sum(idx_bad_match)

    # filter out any outlying data points aka spikes in the dataset (based on per pixel statistics)
    std_outliers = 2.
    galah_solar_median = np.nanmedian(galah_solar, axis=0)
    galah_solar_std = np.nanstd(galah_solar, axis=0)
    idx_outliers = np.abs(galah_solar - galah_solar_median) > std_outliers * galah_solar_std
    galah_solar[idx_outliers] = np.nan
    print 'Removed outliers:', np.sum(idx_outliers)

    # remove zero (or very small) values in the final array - might be the result of final resampling after rv shift
    galah_solar[galah_solar < 0.005] = np.nan

    # save result
    print 'Stacking'
    galah_solar_median = np.nanmedian(galah_solar, axis=0)

    out_file = galah_data_input + 'Solar_data/'
    out_file += 'b' + str(ccd) + '_solar_galah_ext0.txt'
    txt = open(out_file, 'w')
    for i_l in range(len(solar_ref_wvl)):
        txt.write(str(solar_ref_wvl[i_l]) + ' ' + str(galah_solar_median[i_l]) + '\n')
    txt.close()


