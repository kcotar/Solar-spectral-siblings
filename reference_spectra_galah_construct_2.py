from solar_siblings_functions import *
from multiprocessing import Pool
import os

# PC hostname
pc_name = gethostname()

# Galah dir
dr_ver = 'dr5.3'
dr_num = dr_ver[2]+dr_ver[4]

dr_dir = '/shared/ebla/cotar/'+dr_ver+'/'
galah_data_input = '/shared/ebla/cotar/'
results_dir = '/shred/data-camelot/cotar/'
# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E' or pc_name=='new-gigli':
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import get_spectra_dr52, spectra_normalize, spectra_resample


# ---------------------
# FUNCTIONS
# ---------------------
def read_and_prepare_spectra(s_obj):
    print s_obj  # , dr_dir, reduce_bands
    # read all spectral bands
    flux, wvl = get_spectra_dr52(str(s_obj), bands=reduce_bands, root=dr_dir, extension=read_ext, individual=False)

    if len(flux) == 0:
        flux_out = list([])
        for i_r in range(len(reduce_bands)):
            ccd = reduce_bands[i_r] - 1
            idx_ref_sub = np.logical_and(solar_ref_wvl_all >= min_wvl[ccd], solar_ref_wvl_all <= max_wvl[ccd])
            flux_out.append(np.full(np.sum(idx_ref_sub), np.nan))
        return np.hstack(flux_out)

    if read_ext == 0:
        # renormalize them if needed
        for i_b in range(len(flux)):
            norm_ok_mask = determine_norm_mask(wvl[i_b], norm_bad_ranges)
            flux[i_b] = spectra_normalize(wvl[i_b]-np.mean(wvl[i_b]), flux[i_b], fit_mask=norm_ok_mask,
                                          steps=11, sigma_low=2., sigma_high=3., order=7, n_min_perc=5.,
                                          return_fit=False, func='cheb')

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


# ---------------------
# SETTINGS
# ---------------------
i_man_sel = 9
# iDR2
# teff_solar_c_MANUAL = [5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000][i_man_sel]
# logg_solar_c_MANUAL = [4.54, 4.53, 4.51, 4.48, 4.45, 4.41, 4.36, 4.31, 4.25, 4.18][i_man_sel]
# iDR3
teff_solar_c_MANUAL = [5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000][i_man_sel]
logg_solar_c_MANUAL = [4.55, 4.53, 4.51, 4.48, 4.46, 4.43, 4.40, 4.37, 4.34, 4.30][i_man_sel]

feh_solar_c_MANUAL = 0.0
d_teff = 60
d_logg = 0.05
d_feh = 0.05
unlike_ref_suffix = '_{:04.0f}_{:01.2f}_{:01.2f}'.format(teff_solar_c_MANUAL, logg_solar_c_MANUAL, feh_solar_c_MANUAL)

# data-table settings
print 'Reading GALAH and Cannon parameters and data'
data_date = '20180327'
galah_param_file = 'sobject_iraf_'+dr_num+'_reduced_'+data_date+'.fits'
# cannon_date = '180325'
# cannon_param_file = 'sobject_iraf_iDR2_'+cannon_date+'_cannon.fits'
cannon_date = '181221'
cannon_param_file = 'GALAH_iDR3_v1_'+cannon_date+'_cannon.fits'


# reference solar spectra
print 'Read reference Solar spectra'
galah_data_input_solar = galah_data_input + 'Solar_data_dr'+dr_num+'/'
solar_ref = pd.read_csv(galah_data_input_solar + 'solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values
solar_ref_flux_all = solar_ref[:, 1]
solar_ref_wvl_all = solar_ref[:, 0]
galah_data_output_solar = galah_data_input + 'Galah_ref_spectra_dr'+dr_num+'/'
os.system('mkdir '+galah_data_output_solar)

# select ok objects
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)

idx_r1 = np.logical_and(cannon_param['flag_cannon'] == 0,
                        np.abs(cannon_param['Teff_cannon'] - teff_solar_c_MANUAL) <= d_teff)
idx_r2 = np.logical_and(np.abs(cannon_param['Logg_cannon'] - logg_solar_c_MANUAL) <= d_logg,
                        np.abs(cannon_param['Fe_H_cannon'] - feh_solar_c_MANUAL) <= d_feh)

idx_cannon_read = np.logical_and(idx_r1, idx_r2)

print 'Number of parameter-alike spectra that will be used:', np.sum(idx_cannon_read)

read_ext = 4
reduce_bands = list([1, 2, 3, 4])

# multiprocessing for faster processing on cluster
n_multi = 10
pool = Pool(processes=n_multi)
sobjects = cannon_param[idx_cannon_read]['sobject_id']
galah_twilight_spectra_list = pool.map(read_and_prepare_spectra, sobjects)
pool.close()

# define final wavelengths of the twilight spectrum
galah_twilight_spectra_list = np.vstack(galah_twilight_spectra_list)
galah_twilight_spectra_wvl = get_wvl_regions(solar_ref_wvl_all, reduce_bands)
print 'Shape spectra stack:', galah_twilight_spectra_list.shape

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


out_file = galah_data_output_solar
out_file += 'galah_ref_spectrum'+unlike_ref_suffix+'_ext'+str(read_ext)+'_cannon_'+cannon_date+'.txt'
txt = open(out_file, 'w')
for i_l in range(len(galah_twilight_spectra_wvl)):
    txt.write(str(galah_twilight_spectra_wvl[i_l]) + ' ' + str(galah_solar_median[i_l])+ ' ' + str(galah_solar_std[i_l]) + '\n')
txt.close()
