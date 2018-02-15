from solar_siblings_functions import *


# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
process_bands = np.array([1,2,3,4])  # in range 1...4

d_wvl = 0.0
save_plots = False
min_wvl = min_wvl[process_bands-1]
max_wvl = max_wvl[process_bands-1]

GP_compute = False
save_gp_params = True
n_threads = 20
n_walkers = np.array([2*n_threads, 2*n_threads, 2*n_threads, 2*n_threads])[process_bands-1]
n_steps = np.array([40, 40, 40, 40])[process_bands-1]

# evaluate spectrum
n_noise_samples = 0
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
idx_solar_like = (np.abs(cannon_param['Teff_cannon'] - teff_solar_c) < 200) & \
                 (np.abs(cannon_param['Logg_cannon'] - logg_solar_c) < 0.3) & \
                 (np.abs(cannon_param['Feh_cannon'] - feh_solar_c) < 0.3)
#
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['red_flag'] == 0)
# snr selection
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_guess'] > 20)
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['sobject_id'] > 140301000000000)

n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = cannon_param['sobject_id'][idx_solar_like]
sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski','wchebyshev','sqeuclidean','euclidean','chi2', 'EW']
sim_metrices_std = [m+'_std' for m in sim_metrices]
sim_dtypes = ['float64' for i in range(2*len(sim_metrices))]
sim_results = Table(names=np.hstack(('sobject_id', sim_metrices, sim_metrices_std)),
                    dtype=(np.hstack(('int64', sim_dtypes))))

bands_suffix = '_b'+''.join([str(b) for b in process_bands])
dir_suffix = '_p'+str(noise_power)+'_SNRsamples'+str(n_noise_samples)
if GP_compute:
    txt_out = 'GP_fit_res.txt'
    move_to_dir(out_dir + 'Distances_Step2' + dir_suffix)
else:
    move_to_dir(out_dir + 'Distances_Step1' + dir_suffix)

file_out_fits = 'solar_similarity'+bands_suffix+'.fits'

# predetermined objects
# solar_like_sobjects = [
# 140312004501092,140413003201328,140415002401342,140609002101076,140710000101111,140710000101230,140710000101252,140710003901284,140711002401169,140713004001112,140805002601115,140805003101236,140806003501169,140807000601229,140808003701254,140808003701261,140808004701129,140809003101032,140812003201380,140813003001168,150102003201049,150103003001048,150205003401268,150207002601073,150208005201146,150211003701046,150211004701179,150330002601242,150409002101311,150409003101073,150411003601395,150412002101135,150412002601280,150427004301259,150427004801275,150602002701341,150607005101180,150703005601193,150705001901052,150705001901154,150705004901230,150705005401363,150705006401314,150827005201160,150828003201199,150828003701040,150828004701096,150828005201075,150828005701069,150829002601283,150829003101225,150830002801072,150830005101026,150830005101118,150830005101337,151008003501085,151111002101170,151219003601204,151219003601245,151225002701158,151225003801144,151227004701247,151230003201236,160108003601124,160123002601079,160125001601112,160125003001258,160125004501038,160129005201070,160130003601121,160130004101385,160130004601081,160130005201356,160130006301234,160325004201164,160326000101035,160327003601167,160327003601386,160327004101056,160330001601064,160330002101179,160330002601230,160330002601251,160331002701324,160331005301201,160331005301347,160401003901032,160401004401168,160401005401170,160402006101109,160402006101236,160422003501107,160424003601095,160424004701331,160426004501395,160426005501291,160426006101272,160513001601022,160513001601328,160513001601330,160513002101129,160513002601073,160514003301310,160520004201135,160520004901083,160522005101386,160524006601258,160530003901327,160530005001259,160530005501111,160531001601307,160531003101325,160613001801073,160723002601365,160812002601097,160812003601188,160812003601268,160813003101162,160817002601053,160817002601198,160916003801105,160916004301242,160919001601330,160919004001328,160919004601215,160919004601278,160919005101297,160923003001123,160923003701363,160923005201077,161006004901385,161007002801208,161007002801316,161007002801397,161008003001154,161009004801199,161011002801251,161011003401169,161011004001038,161013004401003,161013005401277,161105004601245,161106003601199,161108002101205,161109003901144,161109004901196,161118002601176,161211002601193,161212002601330,161213003101115,161213004101187,161213004601193,161213004601373,161217002101218,161217002101270,161217002101291,161217002601138,161217004601261,161217005101123,161217006101046,161218003601197,161219001801223,161219003601287,161219005101228,161219005101286,170102001901072,170102001901154,170108004601041,170109002801065,170111001601329,170113002601194,170114004101329,170115002201246,170117003101044,170118003801037,170118003801289,170119003601394,170122002101017,170122002601124,170202001201178,170202001201286,170202001701249,170205004401235,170205005401275,170206004701096,170206005701164,170206005701276,170220002101322,170220004101045,170403001601115,170404001601158,170407002101137,170407002101395,170407004101131,170407004601313,170407005201023,170417003201330,170418002701169,170418003701125,170506003401249,170506003901391,170510003301024,170510004301213,170510007301226,170511000101097,170511004001234,170512000101230,170514002401099,170514003001180,170514003301001,170514003301011,170514003301312,170515003101035,170515003101036,170515003101127,170516002101273,170517001801249,170517001801323,170530002601252,170531002801277,170601002601293,170614005101328,170710002701354,170710003201145,170711004001385,170711005101185,170713005101151,170724004601055,170805005101295,170829001901204,170905001601188,170905002101238,170905003101178,170906001601022,170906003101214,170906003601178,170906003601210,170906003601259,170907004601162,170907005201317,170908001601351,170909001601114,170909001601310,170909002601291,170910003101074,170910005101263,170911002101388,170911003601297
# ]

# random subset objects from parameters selection
# n_rand = 25
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(n_rand)*len(solar_like_sobjects))]

for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
    galah_object = galah_param[galah_param['sobject_id'] == s_obj]
    # get spectra of all bands for observed objects
    read_ext = 0
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
                flux[i_c] = spectra_normalize(wvl[i_c]-np.mean(wvl[i_c]), flux[i_c], fit_mask=norm_ok_mask,
                                              steps=15, sigma_low=2., sigma_high=3., order=11, n_min_perc=5.,
                                              return_fit=False, func='cheb')
            # apply computed rv shift to the spectrum
            rv_shift = galah_object['rv_guess_shift']
            wvl *= (1 - rv_shift / 299792.458)
        except:
            print ' -> Something wrong with spectra or reading'
            continue

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
            gp_res_string = str(s_obj) + ',' + str(galah_object['snr_c2_guess'].data[0]) + ','.join([str(v) for v in np.array(gp_final_res).flatten()])
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
    for i_snr in range(n_distances_compute):
        # determine weights for the distance computation (different approaches)
        spectrum_distances[i_snr, :] = compute_distances(pix_spec, pix_std, pix_ref_noise[i_snr, :] + pix_ref, d=noise_power)
        if save_plots:
            plt.figure(1, figsize=(12, 7))
            axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
            axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
            axSpectra.plot(pix_ref_noise[i_snr, :] + pix_ref, lw=0.2, alpha=0.01, c='blue')

    # add agregated results to final table
    sim_results.add_row(np.hstack([s_obj, np.nanmean(spectrum_distances, axis=0), np.nanstd(spectrum_distances, axis=0)]))
    if save_plots:
        axSpectra.plot(pix_ref, c='black', lw=0.5)
        axSpectra.plot(pix_spec, c='blue', lw=0.5)
        axDiff.axhline(y=0, c='black', lw=0.5)
        axDiff.plot(pix_ref-pix_spec, c='blue', lw=0.5)
        axSpectra.set(ylim=(0.3, 1.1))
        axDiff.set(ylim=(-0.05, 0.05))
        plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_guess'].data[0])+bands_suffix+'_flat-new.png', dpi=550)
        plt.close()

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)

'''
print sim_results
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
