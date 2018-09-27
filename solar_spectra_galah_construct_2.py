from solar_siblings_functions import *
from multiprocessing import Pool

# PC hostname
pc_name = gethostname()

# Galah dir
dr_ver = 'dr5.3'
dr_num = dr_ver[2]+dr_ver[4]

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr_dir = '/media/storage/HERMES_REDUCED/'+dr_ver+'/'
    galah_data_input = '/home/klemen/data4_mount/'
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
else:
    galah_data_input = '/data4/cotar/'
    dr_dir = galah_data_input + dr_ver + '/'
from helper_functions import get_spectra_dr52, spectra_normalize, spectra_resample

# some settings
merge_all_dates = True

# ---------------------
# FUNCTIONS
# ---------------------

# reference solar spectra
print 'Read reference Solar spectra'
galah_data_input_solar = galah_data_input + 'Solar_data_dr'+dr_num+'/'
solar_ref = pd.read_csv(galah_data_input_solar + 'solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values
solar_ref_flux_all = solar_ref[:, 1]
solar_ref_wvl_all = solar_ref[:, 0]

# data-table settings
data_date = '20180327'
galah_param_file = 'sobject_iraf_'+dr_num+'_reduced_'+data_date+'.fits'

# select ok objects
galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64,  # also removes wavelength calibration and molecfit problems
                          galah_param['snr_c2_guess'] > 210)  # get only spectra of good quality
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)
galah_param = galah_param[idx_rows]
to_read_row = np.where(idx_rows)[0]
print 'Number of solar spectra that will be used:', len(to_read_row)

# galah_param = galah_param[list(np.int64(np.random.rand(500)*len(galah_param)))]

read_ext = 0
reduce_bands = list([1, 2, 3, 4])

if merge_all_dates:
    possible_unique = ['all']
else:
    print 'Unique dates:', len(np.unique(galah_param['date']))
    galah_param['cob_id'] = np.int64(galah_param['sobject_id'] / 1e5)
    print 'Unique fields:', len(np.unique(galah_param['cob_id']))
    unique_merge_field = 'cob_id'  # field_id or date or something else
    possible_unique = np.unique(galah_param[unique_merge_field])

for date_sel in possible_unique:
    if merge_all_dates:
        flats_sobjects = galah_param['sobject_id']
    else:
        print 'For unique field (' + str(date_sel) + '): ' + str(len(galah_param[unique_merge_field] == date_sel))
        galah_param_sub = galah_param[galah_param[unique_merge_field] == date_sel]
        # select n with best snr per observation date if there is many spectra with good snr
        flats_sobjects = galah_param_sub[np.argsort(galah_param_sub['snr_c2_guess'])[::-1]]['sobject_id'][:100]

    # galah_twilight_spectra_list = list()
    # for s_obj in flats_sobjects:
    def read_and_prepare_twilight_spectra(s_obj):
        print s_obj
        # read all spectral bands
        flux, wvl = get_spectra_dr52(str(s_obj), bands=reduce_bands, root=dr_dir, extension=read_ext, individual=False)

        if read_ext == 0:
            # renormalize them if needed
            for i_b in range(len(flux)):
                norm_ok_mask = determine_norm_mask(wvl[i_b], norm_bad_ranges)
                flux[i_b] = spectra_normalize(wvl[i_b]-np.mean(wvl[i_b]), flux[i_b], fit_mask=norm_ok_mask,
                                              steps=11, sigma_low=2., sigma_high=3., order=7, n_min_perc=5.,
                                              return_fit=False, func='cheb')

                # fit_res = spectra_normalize(wvl[i_b] - np.mean(wvl[i_b]), flux[i_b], fit_mask=norm_ok_mask,
                #                               steps=15, sigma_low=2., sigma_high=3., order=11, n_min_perc=5.,
                #                               return_fit=True, func='cheb')
                # plt.plot(wvl[i_b], fit_res, c='red', lw=2)
                # plt.plot(wvl[i_b], flux[i_b], c='black', lw=1)
                # plt.show()
                # plt.close()
                # continue

            # determine radial velocity of the object and shift spectra accordingly
            wvl_log_shift = get_wvl_log_shift(flux, wvl, solar_ref_flux_all, solar_ref_wvl_all, reduce_bands)
            # apply shift
            flux, wvl = do_wvl_shift_log_and_resample(flux, wvl, wvl_log_shift, solar_ref_wvl_all, reduce_bands)

        else:
            # only resampling is needed
            for i_r in range(len(reduce_bands)):
                ccd = reduce_bands[i_r] - 1
                idx_ref_sub = np.logical_and(solar_ref_wvl_all >= min_wvl[ccd], solar_ref_wvl_all <= max_wvl[ccd])
                ref_wvl_sub = solar_ref_wvl_all[idx_ref_sub]

                obs_flux_res = spectra_resample(flux[i_r], wvl[i_r], ref_wvl_sub, k=1)
                # save results
                flux[i_r] = obs_flux_res
                wvl[i_r] = ref_wvl_sub
        # store resampled and corrected data to the array or list
        # galah_twilight_spectra_list.append(np.hstack(flux))
        return np.hstack(flux)

    # multiprocessing for faster processing on cluster
    n_multi = 10
    pool = Pool(processes=n_multi)
    galah_twilight_spectra_list = pool.map(read_and_prepare_twilight_spectra, flats_sobjects)
    pool.close()

    # galah_twilight_spectra_list = list([])
    # for f_id in flats_sobjects:
    #     galah_twilight_spectra_list.append(read_and_prepare_twilight_spectra(f_id))

    # define final wavelengths of the twilight spectrum
    galah_twilight_spectra_list = np.vstack(galah_twilight_spectra_list)
    galah_twilight_spectra_wvl = get_wvl_regions(solar_ref_wvl_all, reduce_bands)

    print 'Stacking'
    # filter out outlying spectra (based on supplied reference solar spectra)
    solar_ref_flux_galah_sub = solar_ref_flux_all[np.in1d(solar_ref_wvl_all, galah_twilight_spectra_wvl)]
    eucl_match = np.sqrt(np.nansum((galah_twilight_spectra_list - solar_ref_flux_galah_sub) ** 2, axis=1) / np.sum(np.isfinite(galah_twilight_spectra_list), axis=1))
    eucl_match_thr = np.nanpercentile(eucl_match, 97.)
    idx_bad_match = eucl_match >= eucl_match_thr
    galah_twilight_spectra_list[np.where(idx_bad_match)[0], :] = np.nan
    print 'Removed rows:', np.sum(idx_bad_match)

    # filter out any outlying data points aka spikes in the dataset (based on per pixel statistics)
    std_outliers = 2.5
    galah_solar_median = np.nanmedian(galah_twilight_spectra_list, axis=0)
    galah_solar_std = np.nanstd(galah_twilight_spectra_list, axis=0)
    idx_outliers = np.abs(galah_twilight_spectra_list - galah_solar_median) > std_outliers * galah_solar_std
    galah_twilight_spectra_list[idx_outliers] = np.nan
    print 'Removed outliers:', np.sum(idx_outliers)

    # remove zero (or very small) values in the final array - might be the result of final resampling after rv shift
    galah_twilight_spectra_list[galah_twilight_spectra_list < 0.005] = np.nan

    # creation of final solar spectrum per Galah band
    galah_solar_median = np.nanmedian(galah_twilight_spectra_list, axis=0)
    galah_solar_std = np.nanstd(galah_twilight_spectra_list, axis=0)

    # plt.plot(solar_ref_wvl_all, solar_ref_flux_all)
    # plt.plot(galah_twilight_spectra_wvl, galah_solar_median)
    # plt.show()
    # plt.close()

    out_file = galah_data_input_solar
    out_file += 'twilight_spectrum_galah_ext'+str(read_ext)+'_date'+str(date_sel)+'_p7.txt'
    txt = open(out_file, 'w')
    for i_l in range(len(galah_twilight_spectra_wvl)):
        txt.write(str(galah_twilight_spectra_wvl[i_l]) + ' ' + str(galah_solar_median[i_l])+ ' ' + str(galah_solar_std[i_l]) + '\n')
    txt.close()
