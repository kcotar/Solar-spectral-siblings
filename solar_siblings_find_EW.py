from solar_siblings_functions import *

# PC hostname
pc_name = gethostname()

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr52_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
    galah_data_input = '/home/klemen/data4_mount/'
    out_dir = ''
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
else:
    dr52_dir = '/data4/cotar/dr5.2/'
    out_dir = '/data4/cotar/'
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
d_wvl = 0.0
save_plots = False
min_wvl = list([4725, 5665, 6485, 7700])
max_wvl = list([4895, 5865, 6725, 7875])

GP_compute = False
save_gp_params = True
n_threads = 20
n_walkers = [2*n_threads, 2*n_threads, 2*n_threads, 2*n_threads]
n_steps = [40, 40, 40, 40]

# evaluate spectrum
n_noise_samples = 100
noise_power = 0

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix = '_ext0_2_offset'
solar_input_dir = galah_data_input+'Solar_data/'
solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix)

# downscale Solar spectra for faster processing
every_nth_pixel = 8
solar_wvl = solar_wvl[::every_nth_pixel]
solar_flx = solar_flx[::every_nth_pixel]

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'
cannon_param_file = 'sobject_iraf_iDR2_171103_cannon.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 50)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# linelist subset if nedeed
# galah_linelist = galah_linelist[galah_linelist['Element'] == 'Fe']

# linelist mask
idx_lines_mask = solar_wvl < 0.
for line in galah_linelist:
    idx_lines_mask[np.logical_and(solar_wvl >= line['line_start']-d_wvl, solar_wvl <= line['line_end']+d_wvl)] = True
print 'Linelist mask pixels', np.sum(idx_lines_mask)

# find Solar parameters
teff_solar = np.nanmedian(galah_param[idx_rows]['teff_guess'])
teff_solar_std = np.nanstd(galah_param[idx_rows]['teff_guess'])
logg_solar = np.nanmedian(galah_param[idx_rows]['logg_guess'])
logg_solar_std = np.nanstd(galah_param[idx_rows]['logg_guess'])
feh_solar = np.nanmedian(galah_param[idx_rows]['feh_guess'])
feh_solar_std = np.nanstd(galah_param[idx_rows]['feh_guess'])
print 'Solar parameters - guess:', teff_solar, '+/-', teff_solar_std, ',  ', logg_solar, '+/-', logg_solar_std, ',  ', feh_solar, '+/-', feh_solar_std

# same for Cannon
idx_row_cannon = np.in1d(cannon_param['sobject_id'],galah_param[idx_rows]['sobject_id'])
teff_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Teff_cannon'])
teff_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Teff_cannon'])
logg_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Logg_cannon'])
logg_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Logg_cannon'])
feh_solar_c = np.nanmedian(cannon_param[idx_row_cannon]['Feh_cannon'])
feh_solar_std_c = np.nanstd(cannon_param[idx_row_cannon]['Feh_cannon'])
print 'Solar parameters - cannon:', teff_solar_c, '+/-', teff_solar_std_c, ',  ', logg_solar_c, '+/-', logg_solar_std_c, ',  ', feh_solar_c, '+/-', feh_solar_std_c

# Search for objects with similar physical properties
# atomatic parameter selection
# idx_solar_like = np.logical_and(np.abs(galah_param['teff_guess'] - teff_solar) < teff_solar_std*1.5,
#                                 np.abs(galah_param['logg_guess'] - logg_solar) < logg_solar_std*1.5)
# manual parameter selection
idx_solar_like = (np.abs(cannon_param['Teff_cannon']-5605) < 200) & (np.abs(cannon_param['Logg_cannon']-4.21) < 0.3) & (np.abs(cannon_param['Feh_cannon']-(-0.14)) < 0.3)
#
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['red_flag'] == 0)
# snr selection
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_iraf'] > 15)
# idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_iraf'] <= 25)
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['sobject_id'] > 140301000000000)

n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = cannon_param['sobject_id'][idx_solar_like]
sim_metrices = list(np.unique(galah_linelist['Element']))
sim_metrices_std = [m+'_std' for m in sim_metrices]
sim_dtypes = ['float64' for i in range(2*len(sim_metrices))]
sim_results = Table(names=np.hstack(('sobject_id', sim_metrices, sim_metrices_std)),
                    dtype=(np.hstack(('int64', sim_dtypes))))

file_out_fits = 'solar_similarity_narrow.fits'

dir_suffix = '_weighted-p'+str(noise_power)+'_EW_all'
if GP_compute:
    move_to_dir(out_dir + 'Distances_GP' + dir_suffix)
else:
    move_to_dir(out_dir + 'Distances_SNRadded' + dir_suffix)


# solar_like_sobjects = [
#
# ]

# n_rand = 2000
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(n_rand)*len(solar_like_sobjects))]

txt_out = 'GP_fit_res.txt'
for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
    galah_object = galah_param[galah_param['sobject_id'] == s_obj]
    # get spectra of all bands for observed objects
    read_ext = 0
    # flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext, individual=False)
    flux, wvl, flux_std = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext,
                                           individual=False, read_sigma=True)
    if len(flux) <= 0:
        continue
    if read_ext == 0:
        # normalize flux
        try:
            for i_c in range(4):
                # ------ NORM v1 - high order polynomial, many steps
                # flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                #                               n_min_perc=3.,  return_fit=False, func='poly')
                # ------ NORM v2 - the same as used in the process of reference Solar spectra construction
                flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=5, sigma_low=1.5, sigma_high=2.5, order=1,
                                              n_min_perc=5., return_fit=False, func='poly')
                flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=25, sigma_low=1.8, sigma_high=3., order=17,
                                              n_min_perc=5., return_fit=False, func='poly')
            # apply computed rv shift to the spectrum
            rv_shift = galah_object['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        except:
            print ' -> Something wrong with spectra or reading'
            continue

    pix_wvl = list([])
    pix_ref = list([])
    pix_ref_noise = list([])
    pix_spec = list([])
    pix_std = list([])
    if GP_compute:
        gp_final_res = list([])
        # Start GP process for every band in spectrum independently

        for i_c in range(4):
            # define subset of spectra to be compared to reference solar spectrum
            idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]

            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            # correct flux values in needed
            flux_b_res[flux_b_res > 1.15] = 1.15
            flux_b_res[flux_b_res < 0] = 0.

            # determine spectrum difference and its variance
            diff = (solar_flx[idx_ref] - flux_b_res)
            diff_var = np.nanvar(diff)

            # determine kernel parameters trough emcee fit
            print ' Running emcee'
            # emcee_fit_px = 100
            sampler, fit_res, fit_prob = fit_gp_kernel([diff_var/2., 0.0025, 1e-5, 15],
                                                       diff, solar_wvl[idx_ref], data_std=None,  # data_std=flux_std_b_res,
                                                       nwalkers=n_walkers[i_c], n_threds=n_threads, n_burn=n_steps[i_c])

            # walker prob plot
            if save_plots:
                print(" Plotting walker probabilities")
                walkers_prob = sampler.lnprobability
                for i_w in range(walkers_prob.shape[0]):
                    plt.plot(walkers_prob[i_w, :])
                plt.savefig(str(s_obj) + '_gp-lnprob_b' + str(i_c + 1) + '.png', dpi=400)
                # plt.show()
                plt.close()

            sampler_chain_vals = sampler.flatchain
            kernel_fit = np.median(sampler_chain_vals, axis=0)  # flatchain holds parameters of all emcee steps
            # kernel_fit = np.median(fit_res, axis=0)  # fit_res holds only the parameters of the last step
            # kernel_fit = fit_res[np.argmax(fit_prob)]
            # print 'Median val:', np.median(sampler_chain_vals, axis=0)
            # print 'Max lnprob:', fit_res[np.argmax(fit_prob)]
            # print diff_var
            # print kernel_fit
            gp_final_res.append(kernel_fit)

            # corner plot of parameters
            if save_plots:
                c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                                      labels=['amp_noise', 'rad_noise', 'amp_cont', 'rad_cont'], bins=30)
                c_fig.savefig(str(s_obj)+'_corner_b'+str(i_c+1)+'.png', dpi=400)
                plt.close(c_fig)

            # create a gaussian process that will be used for the whole spectra
            gp = george.GP(get_kernel(kernel_fit))
            gp.compute(solar_wvl[idx_ref])
            gp_noise_pred = gp.sample(size=n_noise_samples)

            if save_plots:
                plt.plot(solar_flx[idx_ref], c='red', lw=0.5)
                for i_pred in range(20):
                    plt.plot(solar_flx[idx_ref] + gp_noise_pred[i_pred, :], c='black', alpha=0.15, lw=0.3)
                plt.plot(flux_b_res, c='blue', lw=0.5)
                plt.ylim((0.4, 1.1))
                # plt.show()
                plt.savefig(str(s_obj)+'_gp_b'+str(i_c+1)+'.png', dpi=550)
                plt.close()

            pix_ref.append(solar_flx[idx_ref])
            pix_ref_noise.append(gp_noise_pred)
            pix_spec.append(flux_b_res)
            pix_std.append(flux_std_b_res)

        # save fit res
        if save_gp_params:
            txt = open(txt_out, 'a')
            gp_res_string = str(s_obj) + ',' + str(galah_object['snr_c2_iraf'].data[0]) + ','.join([str(v) for v in np.array(gp_final_res).flatten()])
            txt.write(gp_res_string + '\n')
            txt.close()

    else:
        for i_c in range(4):
            idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
            abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
            # print flux[i_c], wvl[i_c], solar_wvl[idx_ref]
            flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
            flux_std_b_res = spectra_resample(flux_std[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

            # generate poissonian noise to make a spectrum with snr into a spectrum with target snr
            snr_ref = np.inf
            snr_spectrum = galah_object['snr_c'+str(i_c+1)+'_iraf'].data

            snr_sigma = np.sqrt((1.0 / snr_spectrum)**2)  #- (1.0 / snr_ref) ** 2)
            snr_noise_pred = np.random.poisson((1.0 / snr_sigma)**2, size=(n_noise_samples, len(flux_b_res)))
            snr_noise_pred = snr_noise_pred / ((1.0 / snr_sigma)**2) - 1.

            # generate noise to observed spectrum based on given snr value of the spectrum
            pix_wvl.append(solar_wvl[idx_ref])
            pix_ref.append(solar_flx[idx_ref])
            pix_spec.append(flux_b_res)
            pix_std.append(flux_std_b_res)
            pix_ref_noise.append(snr_noise_pred)

    # compute different distance measurements
    pix_wvl = np.hstack(pix_wvl)
    pix_ref = np.hstack(pix_ref)
    pix_ref_noise = np.hstack(pix_ref_noise)
    pix_spec = np.hstack(pix_spec)
    pix_std = np.hstack(pix_std)

    if not evaluate_spectrum(pix_spec, flux_std):
        continue

    # iterate and add noise to observed spectrum
    spectrum_distances = np.zeros((n_noise_samples, len(sim_metrices)))
    for i_snr in range(n_noise_samples):
        # determine weights for the distance computation (different approaches)
        spectrum_distances[i_snr, :] = compute_ews(pix_spec, pix_wvl, pix_ref_noise[i_snr, :] + pix_ref, sim_metrices, galah_linelist)
        if save_plots:
            plt.plot(pix_ref_noise[i_snr, :] + pix_ref, lw=0.2, alpha=0.01, c='blue')
    # add agregated results to final table
    sim_results.add_row(np.hstack([s_obj, np.nanmean(spectrum_distances, axis=0), np.nanstd(spectrum_distances, axis=0)]))
    if save_plots:
        plt.plot(pix_spec, lw=0.5, alpha=1, c='black')
        plt.ylim((0.6, 1.2))
        plt.xlim(0, len(pix_spec))
        plt.tight_layout()
        # plt.show()
        plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_iraf'].data[0])+'.png', dpi=550)
        plt.close()

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)


