from solar_siblings_functions import *
from astropy.modeling import models, fitting

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

# evaluate spectrum
n_noise_samples = 100
noise_power = 0

# reference solar spectra
print 'Read reference GALAH Solar spectra'

suffix_solar_ref = '_ext0_2_offset'
solar_input_dir = galah_data_input+'Solar_data/'

move_to_dir(out_dir + 'Distances_SNR-functions_multioffset_allbands_subsample')

observe_bands = list([1, 2, 3, 4])
observe_snr = range(5, 180, 1)
observe_flux = np.arange(0, 0.21, 0.02)  # float percents  0...1

for evaluate_band in observe_bands:
    print 'Evaluationg band', evaluate_band, every_nth_solar_pixel[evaluate_band-1]
    # solar wvl and flux data
    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=every_nth_solar_pixel[evaluate_band-1])
    # band wvl mask
    idx_band_mask = get_band_mask(solar_wvl, evaluate_band)
    # generate mask of pixels used in comparison
    idx_lines_mask = get_linelist_mask(solar_wvl)
    # print number of pixels
    print ' Band mask pixels:', np.sum(idx_band_mask)
    print ' Initial linelist mask pixels', np.sum(idx_lines_mask)

    sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski',
                    'wchebyshev', 'sqeuclidean', 'euclidean', 'chi2', 'EW']
    sim_metrices_std = [m + '_std' for m in sim_metrices]
    sim_dtypes = ['float64' for i in range(2 * len(sim_metrices))]

    list_fit_results = list([])

    observe_masks = [len(idx_lines_mask)]  # whole range
    # observe_masks = np.linspace(0, len(idx_lines_mask), 10)[1:]  # multiple masks

    spectrum_flux_offsets = list([])
    for offset_percent in observe_flux:
        print ' Flux offset:', offset_percent
        filename_suffix = '_b{:.0f}_flux{:.2f}'.format(evaluate_band, offset_percent)
        flux_offset = (1. - solar_flx) * offset_percent
        spectrum_flux_offsets.append(flux_offset)

        file_out_fits = 'solar_similarity_narrow'+filename_suffix+'.fits'
        fit_results_path = 'metrices_snr_function'+filename_suffix+'.csv'
        list_fit_results.append(fit_results_path)

        sim_results = Table(names=np.hstack(('snr', sim_metrices, sim_metrices_std, 'n_px_eval')),
                            dtype=(np.hstack(('float64', sim_dtypes, 'int64'))))

        # generate noise and evaluate distances
        if not os.path.isfile(file_out_fits):
            for l_mask in observe_masks:
                idx_band_lines_mask = np.logical_and(idx_lines_mask, idx_band_mask)
                # idx_band_lines_mask[np.int64(l_mask):] = False
                print '  Lines mask pixels for this band:', l_mask, np.sum(idx_band_lines_mask)
                for snr_spectrum in observe_snr:
                    print '  SNR generated:', snr_spectrum

                    if n_noise_samples == 1:
                        snr_noise_pred = np.zeros((n_noise_samples, len(solar_wvl)))
                    else:
                        snr_sigma = np.sqrt((1.0 / snr_spectrum)**2)
                        snr_noise_pred = np.random.poisson((1.0 / snr_sigma)**2, size=(n_noise_samples, len(solar_wvl)))
                        snr_noise_pred = snr_noise_pred / ((1.0 / snr_sigma)**2) - 1.

                    # iterate and add noise to observed spectrum
                    spectrum_distances = np.zeros((n_noise_samples, len(sim_metrices)))
                    for i_snr in range(n_noise_samples):
                        # determine weights for the distance computation (different approaches)
                        spectrum_distances[i_snr, :] = compute_distances(solar_flx[idx_band_lines_mask], None,
                                                                         (snr_noise_pred[i_snr, :] + solar_flx + flux_offset)[idx_band_lines_mask],
                                                                         d=noise_power)
                        if save_plots:
                            pass
                    # add agregated results to final table
                    sim_results.add_row(np.hstack([snr_spectrum, np.nanmean(spectrum_distances, axis=0),
                                                   np.nanstd(spectrum_distances, axis=0), np.sum(idx_band_lines_mask)]))

                # save results
                sim_results.write(file_out_fits, overwrite=True)
        else:
            sim_results = Table.read(file_out_fits)

        # generate plots and fitted functions
        txt_out = open(fit_results_path, 'w')
        txt_out.write('metric,amplitude,x_0,alpha,y_const,n_mask_px\n')
        for col in sim_metrices:
            # fit function
            f_init = models.PowerLaw1D() + models.Const1D(amplitude=0)
            fitter = fitting.LevMarLSQFitter()
            gg_fit = fitter(f_init, sim_results['snr'], sim_results[col])
            plt.plot(observe_snr, gg_fit(observe_snr), alpha=0.5, label='Fitted curve')
            # plot
            # plt.scatter(sim_results['snr'], sim_results[col],
            #             lw=0, s=3, c='black', alpha=0.5)
            plt.errorbar(sim_results['snr'], sim_results[col], yerr=sim_results[col+'_std'], c='black', alpha=0.75, errorevery=2,
                         capsize=0, elinewidth=0.75, linewidth=0.5, fmt='o', ms=1.5, label='Simulations')  #, lw=0, s=3, c='black', alpha=0.5)
            # output(s)
            plt.xlabel('Simulated SNR')
            plt.ylabel('Similarity meassure')
            plt.legend()
            plt.tight_layout()
            plt.savefig(col + filename_suffix+'.png', dpi=500)
            plt.close()
            # print gg_fit
            # save fitted values
            txt_out.write(col+','+str(gg_fit.amplitude_0.value)+','+str(gg_fit.x_0_0.value)+','+str(gg_fit.alpha_0.value)+','+str(gg_fit.amplitude_1.value)+'\n')
        txt_out.close()

    # plot observed spectra
    plt.figure(figsize=(14, 4))
    for i_f in range(len(observe_flux)):
        plt.plot(solar_wvl, solar_flx + spectrum_flux_offsets[i_f], label='{:.2f}'.format(observe_flux[i_f]), lw=0.5)
    plt.xlim((min_wvl[evaluate_band-1], max_wvl[evaluate_band-1]))
    # plt.ylim((0.65, 1.02))
    # plt.xlim((5701, 5716))
    plt.tight_layout()
    plt.legend()
    plt.savefig('solar_spectra_b'+str(evaluate_band)+'.png', dpi=300)
    plt.close()

    print ' Final fit results comparison'
    # plot all curves with different flux offsets
    for metric in sim_metrices:
        for i_f in range(len(list_fit_results)):
            y_vals = metric_by_snr(Table.read(list_fit_results[i_f]), metric, observe_snr)
            plt.plot(observe_snr, y_vals, lw=1, label='{:.2f}'.format(observe_flux[i_f]))
        plt.legend()
        plt.savefig(metric+'_all_b'+str(evaluate_band)+'.png', dpi=500)
        plt.close()
