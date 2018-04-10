from solar_siblings_functions import *
import sys, getopt
from os.path import isfile
from time import time

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
# parse inputs from command line

process_bands = np.array([1])  # in range 1...4
read_ext = 0
process_obj_begin = 0
process_obj_end = -1
out_dir_suffix = ''

argv = sys.argv
if len(argv) > 1:
    # parse input options
    opts, args = getopt.getopt(argv[1:], '', ['bands=', 'ext=', 'obj_beg=', 'obj_end=', 'dir_suffix='])
    # set parameters, depending on user inputs
    print opts
    for o, a in opts:
        if o == '--bands':
            process_bands = np.array([np.int32(b) for b in a.split(',')])
            print 'Command line selected bands: ' + ','.join([str(pb) for pb in process_bands])
        if o == '--ext':
            read_ext = np.int8(a)
        if o == '--obj_beg':
            process_obj_begin = np.int64(a)
        if o == '--obj_end':
            process_obj_end = np.int64(a)
        if o == '--dir_suffix':
            out_dir_suffix = str(a)

d_wvl = 0.0
save_plots = True
output_differences = False  # so far available only for the first analysis step
min_wvl = min_wvl[process_bands-1]
max_wvl = max_wvl[process_bands-1]

GP_compute = True
GP_per_element_analysis = True
save_gp_params = True
save_gp_median_spectra = True
save_gp_params_read_append = True
save_gp_params_read_interpol = False
n_threads = 11
n_walkers = np.array([32, 32, 32, 32])[process_bands-1]
n_steps = np.array([80, 80, 90, 130])[process_bands-1]

# evaluate spectrum
n_noise_samples = 1000
noise_power = 0

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'
# solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix, every_nth=8)

# data-table settings
data_date = '20180327'
galah_param_file = 'sobject_iraf_53_reduced_'+data_date+'.fits'
cannon_param_file = 'sobject_iraf_iDR2_180325_cannon.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)
# join datasets and add some information to cannon parameters
cannon_param = join(cannon_param, galah_param['sobject_id','snr_c1_guess','snr_c2_guess','snr_c3_guess','snr_c4_guess'],
                    keys='sobject_id')

idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c2_guess'] > 200)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# same for Cannon
idx_row_cannon = np.in1d(cannon_param['sobject_id'], galah_param[idx_rows]['sobject_id'])
idx_row_cannon = np.logical_and(idx_row_cannon, cannon_param['flag_cannon'] == 0)  # only unflagged flats for parameters

teff_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Teff_cannon'])
teff_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Teff_cannon'])
logg_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Logg_cannon'])
logg_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Logg_cannon'])
feh_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Fe_H_cannon'])
feh_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Fe_H_cannon'])
print 'Solar parameters - cannon:', teff_solar_c, '+/-', teff_solar_std_c, ',  ', logg_solar_c, '+/-', logg_solar_std_c, ',  ', feh_solar_c, '+/-', feh_solar_std_c

# manual parameter selection
idx_solar_like = (np.abs(cannon_param['Teff_cannon'] - teff_solar_c) <= 250) & \
                 (np.abs(cannon_param['Logg_cannon'] - logg_solar_c) <= 0.4) & \
                 (np.abs(cannon_param['Fe_H_cannon'] - feh_solar_c) <= 0.3)
# preform flag filtering if needed - later selection is currently implemented
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['flag_cannon'] >= 0)  # no flagging at this point
idx_solar_like = np.logical_and(idx_solar_like, np.bitwise_and(cannon_param['red_flag'], 64) == 0)  # only flats are taken out at this point in processing
# snr selection
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_guess'] >= 0)  # no snr limits
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['sobject_id'] > 140301000000000)  # leave out comissoning phase

n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = cannon_param['sobject_id'][idx_solar_like]
sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski', 'wchebyshev', 'sqeuclidean', 'euclidean', 'chi2', 'EW', 'median_sep', 'sum', 'px_over', 'px_under']
sim_metrices_std = [m+'_std' for m in sim_metrices]
if GP_compute:
    sim_metrices_min = [m+'_min' for m in sim_metrices]
    sim_metrices_max = [m+'_max' for m in sim_metrices]
    sim_dtypes = ['float64' for i in range(4*len(sim_metrices))]
    sim_results = Table(names=np.hstack(('sobject_id', 'snr_spectrum', 'median_cont', 'median_cont_after', sim_metrices, sim_metrices_std, sim_metrices_min, sim_metrices_max)),
                        dtype=(np.hstack(('int64', 'float64', 'float64', 'float64', sim_dtypes))))
else:
    sim_dtypes = ['float64' for i in range(2*len(sim_metrices))]
    sim_results = Table(names=np.hstack(('sobject_id', 'snr_spectrum', 'median_cont', sim_metrices, sim_metrices_std)),
                        dtype=(np.hstack(('int64', 'float64', 'float64', sim_dtypes))))

bands_suffix = '_b'+''.join([str(b) for b in process_bands])
print bands_suffix
dir_suffix = '_p'+str(noise_power)+'_SNRsamples'+str(n_noise_samples)
dir_suffix += '_ext'+str(read_ext)
if OK_LINES_ONLY:
    dir_suffix += '_oklinesonly'
if not USE_SUBSAMPLE:
    dir_suffix += '_origsamp'
dir_suffix += out_dir_suffix
if GP_compute:
    move_to_dir(out_dir + 'Distances_Step2' + dir_suffix)
    gp_param_labels = ['amp_noise', 'rad_noise', 'amp_cont', 'rad_cont', 'cont_norm']
    txt_out_gp = 'GP_fit_res.txt'
    if isfile(txt_out_gp) and not save_gp_params_read_append or not isfile(txt_out_gp):
        # do not rewrite original file
        txt = open(txt_out_gp, 'w')
        txt_out_header_str = 'sobject_id'
        for i_p_b in process_bands:
            txt_out_header_str += ',' + ','.join([p_l+'_b'+str(i_p_b) for p_l in gp_param_labels])
        txt.write(txt_out_header_str + '\n')
        txt.close()
else:
    move_to_dir(out_dir + 'Distances_Step1' + dir_suffix)

file_out_fits = 'solar_similarity'+bands_suffix+'.fits'
file_out_diff = 'solar_spectral_diff'+bands_suffix+'.csv'

if GP_compute:
    file_out_csv_gp = 'solar_similarity' + bands_suffix + '_gp.csv'
    header_csv_sim = ','.join(np.hstack(('sobject_id', 'snr_spectrum', sim_metrices)))+'\n'
    txt_gp = open(file_out_csv_gp, 'w')
    txt_gp.write(header_csv_sim)
    txt_gp.close()
    if GP_per_element_analysis:
        # same header for per eleemnt analysis
        for elem in get_used_elements():
            txt_gp = open(file_out_csv_gp[:-4] + '_' + elem + '.csv', 'w')
            txt_gp.write(header_csv_sim)
            txt_gp.close()

if save_gp_median_spectra:
    gp_wvl_out_root = 'gp_median_wvl_b'
    gp_flux_out_root = 'gp_median_flx_b'
    for p_b in process_bands:
        txt_gp = open(gp_wvl_out_root+str(p_b)+'.csv', 'w')
        txt_gp.close()
        txt_gp = open(gp_flux_out_root+str(p_b)+'.csv', 'w')
        txt_gp.close()

# predetermined objects in a text file
list_file = 'final_selection_0.10.txt'
solar_like_sobjects = []
txt_o = open(list_file, 'r')
for line in txt_o:
    line_split = line.split(',')
    sobj_vals = np.int64([l_s for l_s in line_split if len(l_s) > 5])
    solar_like_sobjects.append(sobj_vals)
txt_o.close()
solar_like_sobjects = list(np.hstack(solar_like_sobjects))
print 'Number of pre-selected objects:', len(solar_like_sobjects)

# random subset objects from parameters selection
# n_rand = 25
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(n_rand)*len(solar_like_sobjects))]

# first erase all results from previous processing runs
if output_differences:
    csv_diff = open(file_out_diff, 'w')
    csv_diff.close()

for s_obj in solar_like_sobjects[process_obj_begin:process_obj_end]:
    print 'Evaluating', s_obj
    time_ss = time()
    galah_object = galah_param[galah_param['sobject_id'] == s_obj]
    # get spectra of all bands for observed objects
    # flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext, individual=False)
    flux, wvl, flux_std = get_spectra_dr52(str(s_obj), bands=process_bands, root=dr52_dir, extension=read_ext,
                                           individual=False, read_sigma=True)
    if len(flux) <= 0:
        continue
    if read_ext == 0:
        # normalize flux
        try:
            for i_c in range(len(process_bands)):
                # ------ NORM v1 - high order polynomial, many steps
                # flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                #                               n_min_perc=3.,  return_fit=False, func='poly')
                # ------ NORM v2 - the same as used in the process of reference Solar spectra construction
                norm_ok_mask = determine_norm_mask(wvl[i_c], norm_bad_ranges)
                flux[i_c] = spectra_normalize(wvl[i_c] - np.mean(wvl[i_c]), flux[i_c], fit_mask=norm_ok_mask,
                                              steps=15, sigma_low=2., sigma_high=3., order=5, n_min_perc=5.,
                                              return_fit=False, func='cheb')
                # # additional normalization step with symmetric sigma rejection intervals to cancel out noise
                # flux[i_c] = spectra_normalize(wvl[i_c] - np.mean(wvl[i_c]), flux[i_c], fit_mask=norm_ok_mask,
                #                               steps=15, sigma_low=2.5, sigma_high=2.5, order=1, n_min_perc=5.,
                #                               return_fit=False, func='poly')
            # apply computed rv shift to the spectrum
            rv_shift = galah_object['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        except:
            print ' -> Something wrong with spectra or reading'
            continue

    # compute guess like snr for particular spectrum and observed region
    # get absorption features indices
    idx_lines_mask = get_linelist_mask(np.hstack(wvl))
    wvl_all_abs = np.hstack(wvl)#[idx_lines_mask]
    flx_all_abs = np.hstack(flux)#[idx_lines_mask]
    # median signal at selected abundance wavelength pixels
    snr_signal = np.nanmedian(flx_all_abs)
    # determine actual snr of generated noise at selected pixels - guess like
    snr_noise = 1.4826 / np.sqrt(2) * np.nanmedian(np.abs(flx_all_abs[1:] - flx_all_abs[:-1]))
    snr_guesslike = snr_signal / snr_noise
    # print 'SNRs:', galah_object['snr_c' + str(i_c + 1) + '_guess'].data[0], snr_guesslike

    # determine continuum-like pixels
    min_cont_level = 0.98
    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=4)
    idx_cont_px = solar_flx > min_cont_level
    flx_all_abs_res = spectra_resample(flx_all_abs, wvl_all_abs, solar_wvl, k=1)
    idx_cont_px = np.logical_and(idx_cont_px, np.isfinite(flx_all_abs_res))
    cont_median_solar = np.nanmedian(solar_flx[idx_cont_px])
    cont_median_flx = np.nanmedian(flx_all_abs_res[idx_cont_px])
    cont_median_dif = cont_median_solar - cont_median_flx

    pix_ref = list([])
    pix_ref_wvl = list([])
    pix_gp_predicted = list([])
    pix_ref_noise = list([])
    pix_ref_cont = list([])  # continuum pixels for continuum offset determination
    pix_spec = list([])
    pix_std = list([])
    if GP_compute:
        gp_final_res = list([])
        # Start GP process for every band in spectrum independently

        for i_c in range(len(process_bands)):
            evaluate_band = process_bands[i_c]
            # first prepare reference data
            solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref,
                                                  every_nth=every_nth_solar_pixel[evaluate_band - 1])
            # band wvl mask
            idx_ref = get_band_mask(solar_wvl, evaluate_band)
            # generate mask of pixels used in comparison
            idx_lines_mask = get_linelist_mask(solar_wvl)

            # define subset of spectra to be compared to reference solar spectrum
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]

            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            # correct flux values if needed
            flux_b_res[flux_b_res > 1.2] = 1.2
            flux_b_res[flux_b_res < 0] = 0.

            # determine spectrum difference and its variance
            diff = (solar_flx[idx_ref] - flux_b_res)
            diff_var = np.nanvar(diff)

            skip_gp_computation = False
            # check if GP results are already available:
            if (isfile(txt_out_gp) and save_gp_params_read_append) or save_gp_params_read_interpol:
                try:
                    gp_res_read_cols = [gp_p_l + '_b' + str(evaluate_band) for gp_p_l in gp_param_labels]
                    if save_gp_params_read_interpol:
                        gp_precom_train = Table.read('GP_fit_res_train.txt', format='ascii.csv')
                        idx_sel = np.abs(gp_precom_train['snr_spectrum'] - snr_guesslike) <= 20
                        if np.sum(idx_sel) > 0:
                            kernel_fit = gp_precom_train[idx_sel][gp_res_read_cols].to_pandas().values
                            kernel_fit = np.median(kernel_fit, axis=0)
                            skip_gp_computation = True
                        print ' GP emcee parameters intepolated'
                    else:
                        gp_precom_fit = Table.read(txt_out_gp, format='ascii.csv')
                        idx_line = np.where(gp_precom_fit['sobject_id'] == s_obj)[0]
                        if len(idx_line) == 1:
                            # we found a match, read it for current band
                            skip_gp_computation = True
                            kernel_fit = gp_precom_fit[idx_line][gp_res_read_cols].to_pandas().values[0]
                            print ' GP emcee parameters restored'
                except:
                    print ' Problem restoring GP emcee parameters'

            if not skip_gp_computation:
                # determine kernel parameters trough emcee fit
                print ' Running emcee'
                # print 'diff_var/2:', diff_var/2.
                rad_noise_init = [0.003, 0.005, 0.007, 0.009][evaluate_band - 1]  # start process with different initial values for every band
                init_guess_l = [diff_var-diff_var/4., rad_noise_init-0.001, 1e-5,  5., 0.98]
                init_guess_h = [diff_var+diff_var/4., rad_noise_init+0.001, 1e-4, 15., 1.02]
                sampler, fit_res, fit_prob = fit_gp_kernel(init_guess_l, init_guess_h,
                                                           solar_flx[idx_ref], flux_b_res, solar_wvl[idx_ref],
                                                           nwalkers=n_walkers[i_c], n_threds=n_threads, n_burn=n_steps[i_c],
                                                           exit_lnp=10, normal_dist_guess=False, n_per_burn=n_steps[i_c])
                # walker prob plot
                if save_plots:
                    print(" Plotting walker probabilities")
                    walkers_prob = sampler.lnprobability/len(flux_b_res)
                    for i_w in range(walkers_prob.shape[0]):
                        plt.plot(walkers_prob[i_w, :], lw=0.3)
                    walkers_prob = walkers_prob.flatten()                   # without this correction numpy
                    walkers_prob = walkers_prob[np.isfinite(walkers_prob)]  # percentile may return incorrect -inf value
                    plt.ylim((np.percentile(walkers_prob, 1), np.percentile(walkers_prob, 99.9)))
                    plt.savefig(str(s_obj) + '_gp-lnprob_b' + str(evaluate_band) + '.png', dpi=250)
                    # plt.show()
                    plt.close()

                last_n_steps = 20
                sampler_chain_vals = sampler.flatchain
                kernel_fit = np.median(sampler_chain_vals, axis=0)  # flatchain holds parameters of all emcee steps
                kernel_fit_last_n = np.median(sampler_chain_vals[-last_n_steps*n_walkers[i_c]:, :], axis=0)
                kernel_fit_best = fit_res[np.argmax(fit_prob)]

                # print out different solutions
                print 'Median all   :', kernel_fit
                print 'Median n last:', kernel_fit_last_n
                print 'Best in last :', kernel_fit_best

                # corner plot of parameters
                if save_plots:
                    c_fig = corner.corner(sampler_chain_vals, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                                          labels=gp_param_labels, bins=30)
                    c_fig.savefig(str(s_obj) + '_corner_b' + str(evaluate_band) + '.png', dpi=200)
                    plt.close(c_fig)

                    c_fig = corner.corner(sampler_chain_vals[-last_n_steps*n_walkers[i_c]:, :], truths=kernel_fit_last_n,
                                          quantiles=[0.16, 0.5, 0.84], labels=gp_param_labels, bins=30)
                    c_fig.savefig(str(s_obj) + '_corner_b' + str(evaluate_band) + '_last'+str(last_n_steps)+'.png', dpi=200)
                    plt.close(c_fig)

            # add fitted values to the resulting table
            gp_final_res.append(kernel_fit)

            # create a gaussian process that will be used for the whole spectra
            print ' Creating conditional GP samples'
            # solar_flx_corr = spectrum_offset_norm(kernel_fit[-1:], solar_flx[idx_ref])
            solar_flx_corr = solar_flx[idx_ref]
            solar_flx_class = mean_flux_class(solar_flx_corr)
            gp = george.GP(get_kernel(kernel_fit[:-1]), mean=solar_flx_class)
            gp.compute(solar_wvl[idx_ref], yerr=flux_std_b_res*flux_b_res,)
            # renormalization added to the stellar spectrum itself not solar fo consistency
            flux_b_res_norm = spectrum_offset_norm(kernel_fit[-1:], flux_b_res)
            flux_gp_pred = gp.sample_conditional(flux_b_res_norm, solar_wvl[idx_ref], size=n_noise_samples)

            if save_plots:
                plt.plot(solar_flx_corr, c='red', alpha=0.8, lw=0.5)
                plt.plot(flux_b_res, c='blue', alpha=0.8, lw=0.5)
                plt.plot(np.median(flux_gp_pred, axis=0), c='black', alpha=0.8, lw=0.5)
                for i_pred in range(20):
                    plt.plot(flux_gp_pred[i_pred, :], c='black', alpha=0.1, lw=0.3)
                plt.ylim((0.4, 1.1))
                plt.tight_layout()
                # plt.show()
                plt.savefig(str(s_obj) + '_gp-sample_b' + str(evaluate_band) + '.png', dpi=400)
                plt.close()

            if save_gp_median_spectra:
                txt_gp = open(gp_wvl_out_root + str(evaluate_band) + '.csv', 'w')
                txt_gp.write(','.join([str(sw) for sw in solar_wvl[idx_ref]]))
                txt_gp.close()
                txt_gp = open(gp_flux_out_root + str(evaluate_band) + '.csv', 'a')
                txt_gp.write(','.join([str(sf) for sf in np.median(flux_gp_pred, axis=0)])+'\n')
                txt_gp.close()

            # determine continuum-like pixels after GP fitting was performed
            pix_ref_cont.append(flux_gp_pred[:, np.where(solar_flx_corr > min_cont_level)[0]])

            # store results for current band
            pix_ref.append(solar_flx_corr[abs_lines_cols])
            pix_ref_wvl.append(solar_wvl[idx_ref][abs_lines_cols])
            pix_gp_predicted.append(flux_gp_pred[:, abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])

        # determine continuum median difference after GP fitting was performed
        cont_median_dif_after = np.median(np.hstack(pix_ref_cont)) - cont_median_flx

        # save fit res if they are not already in it
        if save_gp_params and not skip_gp_computation:
            txt = open(txt_out_gp, 'a')
            gp_res_string = str(s_obj) + ',' + ','.join([str(v) for v in np.array(gp_final_res).flatten()])
            txt.write(gp_res_string + '\n')
            txt.close()

    else:
        for i_c in range(len(process_bands)):
            evaluate_band = process_bands[i_c]
            # first prepare reference data
            solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref,
                                                  every_nth=every_nth_solar_pixel[evaluate_band - 1])
            # band wvl mask
            idx_ref = get_band_mask(solar_wvl, evaluate_band)
            # generate mask of pixels used in comparison
            idx_lines_mask = get_linelist_mask(solar_wvl)

            # define subset of spectra to be compared to reference solar spectrum
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
            # print flux[i_c], wvl[i_c], solar_wvl[idx_ref]
            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            if n_noise_samples <= 0:
                # do not add any noise to used reference spectra
                snr_noise_pred = np.zeros((1, len(flux_b_res)))
            else:
                # generate poissonian noise to make a spectrum with snr into a spectrum with target snr
                snr_ref = np.inf
                snr_spectrum = galah_object['snr_c' + str(i_c + 1) + '_guess'].data
                snr_sigma = np.sqrt((1.0 / snr_spectrum) ** 2)  # - (1.0 / snr_ref) ** 2)
                snr_noise_pred = np.random.poisson((1.0 / snr_sigma)**2, size=(n_noise_samples, len(flux_b_res)))
                snr_noise_pred = snr_noise_pred / ((1.0 / snr_sigma)**2) - 1.

            # store results for current band
            pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])
            pix_ref_noise.append(snr_noise_pred[:, abs_lines_cols])

    # compute different distance measurements
    pix_ref = np.hstack(pix_ref)
    pix_spec = np.hstack(pix_spec)
    pix_std = np.hstack(pix_std)
    if GP_compute:
        pix_gp_predicted = np.hstack(pix_gp_predicted)
        pix_ref_wvl = np.hstack(pix_ref_wvl)

    if not evaluate_spectrum(pix_spec, flux_std):
        continue

    # similarity between median gp predicted signal and solar reference
    if GP_compute:
        gp_median_flx = np.nanmean(pix_gp_predicted, axis=0)
        spectrum_distances_gp = compute_distances(gp_median_flx, 0., pix_ref, d=noise_power)
        txt_gp = open(file_out_csv_gp, 'a')
        gp_res_sim_string = str(s_obj)+',' + str(snr_guesslike)+',' + ','.join([str(val) for val in spectrum_distances_gp])
        txt_gp.write(gp_res_sim_string + '\n')
        txt_gp.close()
        # perform the same analysis for all/selected chemical elements
        if GP_per_element_analysis:
            print '  Per element similarity'
            for elem in get_used_elements():
            # for elem in ['Fe']:
                idx_element_mask = get_linelist_mask(pix_ref_wvl, d_wvl=0., element=elem)
                spectrum_element_distances_gp = compute_distances(gp_median_flx[idx_element_mask], 0., pix_ref[idx_element_mask], d=noise_power)
                txt_gp = open(file_out_csv_gp[:-4]+'_'+elem+'.csv', 'a')
                gp_res_sim_string = str(s_obj) + ',' + str(snr_guesslike) + ',' + ','.join([str(val) for val in spectrum_element_distances_gp])
                txt_gp.write(gp_res_sim_string + '\n')
                txt_gp.close()


    # iterate and add noise to observed spectrum
    n_distances_compute = np.max([n_noise_samples, 1])
    spectrum_distances = np.zeros((n_distances_compute, len(sim_metrices)))
    if save_plots:
        plt.figure(1, figsize=(12, 7))
        axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
        axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
    for i_snr in range(n_distances_compute):
        # determine weights for the distance computation (different approaches)
        if GP_compute:
            spectrum_distances[i_snr, :] = compute_distances(pix_gp_predicted[i_snr, :], 0., pix_ref, d=noise_power)
        else:
            spectrum_distances[i_snr, :] = compute_distances(pix_spec, 0., pix_ref, d=noise_power)

        if save_plots:
            axSpectra.plot(pix_gp_predicted[i_snr, :], lw=0.2, alpha=0.01, c='blue')
        if output_differences and not GP_compute:
            csv_diff = open(file_out_diff, 'a')
            if os.path.getsize(file_out_diff) == 0:  # size of file is zero -> add wavelength header info
                csv_diff.write('0,'+','.join([str(sw) for sw in solar_wvl[idx_ref][abs_lines_cols]])+'\n')
            diff_csv_string = ','.join([str(pf) for pf in (pix_spec - pix_ref)])
            csv_diff.write(str(s_obj)+','+diff_csv_string+'\n')
            csv_diff.close()

    # add agregated results to final table
    if GP_compute:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, cont_median_dif_after, np.nanmean(spectrum_distances, axis=0),
                            np.nanstd(spectrum_distances, axis=0), np.nanmin(spectrum_distances, axis=0), np.nanmax(spectrum_distances, axis=0)]))
    else:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, np.nanmean(spectrum_distances, axis=0), 
                            np.nanstd(spectrum_distances, axis=0)]))

    if save_plots:
        axSpectra.plot(pix_ref, c='black', lw=0.5)
        axSpectra.plot(pix_spec, c='blue', lw=0.5)
        axDiff.axhline(y=0, c='black', lw=0.5)
        axDiff.plot(pix_ref - pix_spec, c='green', lw=0.5)
        axDiff.plot(pix_ref - np.nanmean(pix_gp_predicted, axis=0), c='red', lw=0.5)
        axSpectra.set(ylim=(0.3, 1.1))
        axDiff.set(ylim=(-0.05, 0.05))
        axSpectra.set_title('SNR from abs lines: {:.2f}'.format(snr_guesslike))
        plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_guess'].data[0])+bands_suffix+'.png', dpi=300)
        plt.close()

    time_ee = time()
    print 'Total time: '+str((time_ee-time_ss)/60.)+' min'
    print ''

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)

