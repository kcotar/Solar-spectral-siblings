from solar_siblings_functions import *
import sys, getopt

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
# parse inputs from command line

process_bands = np.array([1])  # in range 1...4
read_ext = 0

argv = sys.argv
if len(argv) > 1:
    # parse input options
    opts, args = getopt.getopt(argv[1:], '', ['bands=', 'ext='])
    # set parameters, depending on user inputs
    print opts
    for o, a in opts:
        if o == '--bands':
            process_bands = np.array([np.int32(b) for b in a.split(',')])
            print 'Command line selected bands: ' + ','.join([str(pb) for pb in process_bands])
        if o == '--ext':
            read_ext = np.int8(a)

d_wvl = 0.0
save_plots = False
output_differences = True  # so far available only for the first analysis step
min_wvl = min_wvl[process_bands-1]
max_wvl = max_wvl[process_bands-1]

GP_compute = True
save_gp_params = True
n_threads = 20
n_walkers = np.array([2*n_threads, 2*n_threads, 2*n_threads, 2*n_threads])[process_bands-1]
n_steps = np.array([40, 40, 40, 40])[process_bands-1]

# evaluate spectrum
n_noise_samples = 500
noise_power = 0

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'
# solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix, every_nth=8)

# data-table settings
data_date = '20180214'
galah_param_file = 'sobject_iraf_53_reduced_'+data_date+'.fits'
cannon_param_file = 'sobject_iraf_iDR2_180108_cannon.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)
# join datasets and add some information to cannon parameters
cannon_param = join(cannon_param, galah_param['sobject_id','snr_c1_guess','snr_c2_guess','snr_c3_guess','snr_c4_guess'],
                    keys='sobject_id')

idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c2_guess'] > 100)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# same for Cannon
idx_row_cannon = np.in1d(cannon_param['sobject_id'], galah_param[idx_rows]['sobject_id'])
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
idx_solar_like = (np.abs(cannon_param['Teff_cannon'] - teff_solar_c) < 200) & \
                 (np.abs(cannon_param['Logg_cannon'] - logg_solar_c) < 0.3) & \
                 (np.abs(cannon_param['Feh_cannon'] - feh_solar_c) < 0.3)
#
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['red_flag'] == 0)
# snr selection
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_guess'] > 0)
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['sobject_id'] > 140301000000000)

n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = cannon_param['sobject_id'][idx_solar_like]
sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski','wchebyshev','sqeuclidean','euclidean','chi2', 'EW', 'median_sep']
sim_metrices_std = [m+'_std' for m in sim_metrices]
if GP_compute:
    sim_metrices_min = [m+'_min' for m in sim_metrices]
    sim_metrices_max = [m+'_max' for m in sim_metrices]
    sim_dtypes = ['float64' for i in range(4*len(sim_metrices))]
    sim_results = Table(names=np.hstack(('sobject_id', 'snr_spectrum', 'median_cont', sim_metrices, sim_metrices_std, sim_metrices_min, sim_metrices_max)),
                        dtype=(np.hstack(('int64', 'float64', 'float64', sim_dtypes))))
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
# TEMP suffix only:
dir_suffix += ''
if GP_compute:
    txt_out = 'GP_fit_res.txt'
    move_to_dir(out_dir + 'Distances_Step2' + dir_suffix)
else:
    move_to_dir(out_dir + 'Distances_Step1' + dir_suffix)

file_out_fits = 'solar_similarity'+bands_suffix+'.fits'
file_out_diff = 'solar_spectral_diff'+bands_suffix+'.csv'

# predetermined objects
solar_like_sobjects = [
140312003501087,140312004501092,140413003201328,140611004001251,140711002401243,140711002401316,140713004001112,140713004001188,140805002601115,140805003101303,140805003101338,140806002301392,140806002901279,140806004101159,140807005001277,140808002701299,140808003701261,140809003101322,140809003701146,141102002401182,141102002701267,141103003101379,141231003001209,150102003201119,150103002701017,150103003001048,150204002901346,150207002101161,150208003201286,150211004701179,150408004101278,150408005301184,150409003101129,150409004101294,150427004301259,150427004801275,150504003001111,150602003301072,150602003901271,150607005601279,150703005101366,150703005601062,150703005601082,150705005401363,150705005401364,150705005401377,150705006401129,150705006401314,150827005201080,150828003701040,150828003701062,150828004201034,150828004201212,150828004701096,150828005701069,150829002601283,150829003101225,150830004001040,150830004601175,150830005101144,150830005101337,150830005101391,151009001601351,151009001601363,151109002101254,151110003101178,151111002101236,151219001601078,151219002601248,151219003601245,151219003601298,151225002701158,151225003801118,151230003201236,151231002601074,160110003601039,160112001601056,160123002601062,160123002601079,160125004501038,160125004501256,160125004501389,160130005201357,160130006301234,160325003201071,160326002101077,160327004101056,160327004101248,160327004101343,160327004601337,160327006101355,160328003201282,160330002601095,160402004601106,160402005101147,160402005601084,160402006101388,160403005201158,160422004501162,160424002101194,160424003101290,160424004701369,160426004501264,160426006101338,160426006701393,160513001601022,160522006101314,160522006601193,160524004901098,160524005501270,160524006101090,160524006601258,160525003201154,160530005001359,160531001601307,160531004601362,160531005601362,160531006101153,161008003001092,161009002601018,161009004801144,161011003401169,161104002801361,161104003801323,161105003101019,161105003101345,161105004601015,161107001601132,161107001601133,161116002201296,161116003801308,161118002601176,161118004701023,161118004701221,161119004701207,161210004201073,161210004201315,161213004101187,161217002601138,161217004101075,161217004601271,161219003101122,161219004101351,161219005101191,161219005101228,170102001901072,170105003101138,170108003301041,170108004601041,170109002101106,170109002801065,170112002101294,170112002601348,170112003101027,170112003601298,170114004101329,170115001601273,170115002201381,170119002601393,170121002801292,170130003101184,170130003601019,170219003601351,170220004101215,170508004801312,170509004701096,170509005201063,170510007301226,170514003001180,170514003301001,170514003301011,170515003101035,170515003101036,170515006101289,170614004101220,170614004101381,170614004601055,170614004601061,170614005101328,170710003201082,170710003201145,170711002001243,170711005101185,170711005801034,170713005101388,170723003601118,170723005101385,170724003601147,170724004601055,170801002801345,170801004001010,170801004001215,170805003601318,170829001901039,170905002601072,170905003101147,170906003601147,170906003601159,170906003601210,170906004101067,170908001601149,170909002101136,170909002601291,170910001801313,170910003101074,170910003101081,170910005101082,170911002101145,170911004201366
]
print 'Number of pre-selected objects:', len(solar_like_sobjects)

# random subset objects from parameters selection
# n_rand = 25
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(n_rand)*len(solar_like_sobjects))]

# first erase all results from previous processing runs
if output_differences:
    csv_diff = open(file_out_diff, 'w')
    csv_diff.close()

for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
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
                                              steps=15, sigma_low=2., sigma_high=3., order=11, n_min_perc=5.,
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
    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=4)
    idx_cont_px = solar_flx > 0.98
    flx_all_abs_res = spectra_resample(flx_all_abs, wvl_all_abs, solar_wvl, k=1)
    idx_cont_px = np.logical_and(idx_cont_px, np.isfinite(flx_all_abs_res))
    cont_median_solar = np.nanmedian(solar_flx[idx_cont_px])
    cont_median_flx = np.nanmedian(flx_all_abs_res[idx_cont_px])
    cont_median_dif = cont_median_solar - cont_median_flx

    pix_ref = list([])
    pix_ref_noise = list([])
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

            # determine kernel parameters trough emcee fit
            print ' Running emcee'
            # emcee_fit_px = 100
            sampler, fit_res, fit_prob = fit_gp_kernel([diff_var/2., 0.0025, 1e-5, 15],
                                                       diff, solar_wvl[idx_ref], data_std=None,  # data_std=flux_std_b_res,
                                                       nwalkers=n_walkers[i_c], n_threds=n_threads, n_burn=n_steps[i_c],
                                                       exit_lnp=10)

            # walker prob plot
            if save_plots:
                print(" Plotting walker probabilities")
                walkers_prob = sampler.lnprobability/len(flux_b_res)
                for i_w in range(walkers_prob.shape[0]):
                    plt.plot(walkers_prob[i_w, :])
                plt.ylim((np.nanpercentile(walkers_prob, 1), np.nanpercentile(walkers_prob, 99)))
                plt.savefig(str(s_obj) + '_gp-lnprob_b' + str(i_c + 1) + '.png', dpi=400)
                # plt.show()
                plt.close()

            sampler_chain_vals = sampler.flatchain
            kernel_fit = np.median(sampler_chain_vals, axis=0)  # flatchain holds parameters of all emcee steps
            gp_final_res.append(kernel_fit)

            # corner plot of parameters
            if save_plots:
                c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                                      labels=['amp_noise', 'rad_noise', 'amp_cont', 'rad_cont'], bins=30)

                # # add something to the plots on diagonal of corner plot
                # # extract the axes
                # n_dim_plot = len(kernel_fit)
                # axes = np.array(c_fig.axes).reshape((n_dim_plot,n_dim_plot))
                # # loop over the axes of plots on diagonal
                # for i in range(n_dim_plot):
                #     ax = axes[i, i]
                #     ax.plot()

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

            pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
            pix_ref_noise.append(gp_noise_pred[:, abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])

        # save fit res
        if save_gp_params:
            txt = open(txt_out, 'a')
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

            # generate noise to observed spectrum based on given snr value of the spectrum
            pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
            pix_spec.append(flux_b_res[abs_lines_cols])
            pix_std.append(flux_std_b_res[abs_lines_cols])
            pix_ref_noise.append(snr_noise_pred[:, abs_lines_cols])

    # compute different distance measurements
    pix_ref = np.hstack(pix_ref)
    pix_ref_noise = np.hstack(pix_ref_noise)
    pix_spec = np.hstack(pix_spec)
    pix_std = np.hstack(pix_std)

    if not evaluate_spectrum(pix_spec, flux_std):
        continue

    # iterate and add noise to observed spectrum
    n_distances_compute = np.max([n_noise_samples, 1])
    spectrum_distances = np.zeros((n_distances_compute, len(sim_metrices)))
    if save_plots:
        plt.figure(1, figsize=(12, 7))
        axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
        axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
    for i_snr in range(n_distances_compute):
        # determine weights for the distance computation (different approaches)
        spectrum_distances[i_snr, :] = compute_distances(pix_spec, pix_std, pix_ref_noise[i_snr, :] + pix_ref, d=noise_power)
        if save_plots:
            axSpectra.plot(pix_ref_noise[i_snr, :] + pix_ref, lw=0.2, alpha=0.01, c='blue')
        if output_differences:
            csv_diff = open(file_out_diff, 'a')
            if os.path.getsize(file_out_diff) == 0:  # size of file is zero -> add wavelength header info
                csv_diff.write('0,'+','.join([str(sw) for sw in solar_wvl[idx_ref][abs_lines_cols]])+'\n')
            diff_csv_string = ','.join([str(pf) for pf in (pix_spec - (pix_ref_noise[i_snr, :] + pix_ref))])
            csv_diff.write(str(s_obj)+','+diff_csv_string+'\n')
            csv_diff.close()

    # add agregated results to final table
    if GP_compute:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, np.nanmean(spectrum_distances, axis=0), 
				np.nanstd(spectrum_distances, axis=0), np.nanmin(spectrum_distances, axis=0), np.nanmax(spectrum_distances, axis=0)]))
    else:
        sim_results.add_row(np.hstack([s_obj, snr_guesslike, cont_median_dif, np.nanmean(spectrum_distances, axis=0), 
				np.nanstd(spectrum_distances, axis=0)]))

    if save_plots:
        axSpectra.plot(pix_ref, c='black', lw=0.5)
        axSpectra.plot(pix_spec, c='blue', lw=0.5)
        axDiff.axhline(y=0, c='black', lw=0.5)
        axDiff.plot(pix_ref-pix_spec, c='blue', lw=0.5)
        axSpectra.set(ylim=(0.3, 1.15))
        axDiff.set(ylim=(-0.05, 0.05))
        axSpectra.set_title('SNR from abs lines: {:.2f}'.format(snr_guesslike))
        plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_guess'].data[0])+bands_suffix+'.png', dpi=350)
        plt.close()

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)

'''
print sim_results5
print ''
sobj_id_like = sim_results[np.argsort(sim_results['chi2'])[:75]]['sobject_id']
print ','.join([str(s) for s in sobj_id_like])

print ''
sobj_id_dislike = sim_results[np.argsort(sim_results['chi2'])[-75:]]['sobject_id']
print ','.join([str(s) for s in sobj_id_dislike])

# output a plot of the most solar like spectra
for i_b in range(1, 5):
    for s_obj in sobj_id_like:
        flux, wvl = get_spectra_dr52(str(s_obj), bands=[i_b], root=dr52_dir, extension=read_ext)
        if read_ext == 0:
            # apply the same normalization as in the process of creation of master solar spectrum
            flux[0] = spectra_normalize(wvl[0], flux[0], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                                          n_min_perc=3., return_fit=False, func='poly')
            # apply computed rv shift to the spectrum
            rv_shift = galah_param[galah_param['sobject_id'] == s_obj]['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        plt.plot(wvl[0], flux[0], lw=0.2, c='blue', alpha=0.02)
    plt.plot(solar_wvl, solar_flx, lw=0.2, c='black')
    plt.xlim((min_wvl[i_b-1], max_wvl[i_b-1]))
    plt.ylim((0.4, 1.1))
    plt.savefig('similar_spectra_b'+str(i_b)+'.png', dpi=1000)
    plt.close()
'''
