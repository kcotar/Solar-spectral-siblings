from solar_siblings_functions import *
from os.path import isfile

from multiprocessing import Pool
from time import time

# PC hostname
pc_name = gethostname()

# input data
dr52_dir = '/shared/ebla/cotar/dr5.3/'
out_dir = '/shared/data-camelot/cotar/'
galah_data_input = '/shared/ebla/cotar/'
imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')

from helper_functions import *
from spectra_collection_functions import *


# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
d_wvl = 0.0
save_plots = False

run_unlike_solar = False

teff_solar_c_MANUAL = [6000, 5900, 5800, 5700, 5600, 5500, 5400, 5300, 5200, 5100][2]
logg_solar_c_MANUAL = [4.18, 4.25, 4.31, 4.36, 4.41, 4.45, 4.48, 4.51, 4.53, 4.54][2]
feh_solar_c_MANUAL = 0.0
unlike_ref_suffix = '_{:04.0f}_{:01.2f}_{:01.2f}'.format(teff_solar_c_MANUAL, logg_solar_c_MANUAL, feh_solar_c_MANUAL)
unlike_input_dir = galah_data_input + 'Galah_ref_spectra_dr53/'


# evaluate spectrum
n_noise_samples = 1000
noise_power = 0
test_spectrum_is_line = False

# added noise and evaluation settings
normalize_noisy_spectrum = False
gauss_noise = True
compute_guesslike_snr = True
guesslike_all_lines = True

# reference solar spectra
print 'Read reference GALAH Solar spectra'

suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'

new_dir = 'Distances_SNR_models'
if run_unlike_solar:
    new_dir += unlike_ref_suffix
if compute_guesslike_snr:
    new_dir += '_guesslike'
    if guesslike_all_lines:
        new_dir += '-alllines'
    else:
        new_dir += '-abslines'
if gauss_noise:
    new_dir += '_gauss'
if normalize_noisy_spectrum:
    new_dir += '_renorm'
if test_spectrum_is_line:
    new_dir += '_spectrumisline'
if OK_LINES_ONLY:
    new_dir += '_oklinesonly'
if not USE_SUBSAMPLE:
    new_dir += '_origsamp'
if n_noise_samples == 1:
    new_dir += '_nonoise'

new_dir += '_withH'
move_to_dir(out_dir + new_dir)

observe_bands = list([1, 2, 3, 4])
all_bands_at_once = False
observe_snr = range(5, 170, 2)
observe_flux = [0., 0.05, 0.1, 0.15, 0.2]  #  float percents  0...1
# observe_flux = np.arange(-0.20, 0.21, 0.02)

for evaluate_band in observe_bands:
    print 'Evaluating band', evaluate_band, every_nth_solar_pixel[evaluate_band-1]
    # solar wvl and flux data
    if not run_unlike_solar:
        solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=every_nth_solar_pixel[evaluate_band-1])
    else:
        solar_wvl, solar_flx = get_unlikesolar_data(unlike_input_dir, unlike_ref_suffix, every_nth=every_nth_solar_pixel[evaluate_band-1])

    if test_spectrum_is_line:
        solar_flx[:] = 1.
    # band wvl mask
    if all_bands_at_once:
        idx_band_mask = np.full(len(solar_wvl), False)
        for e_b_mask in observe_bands:
            idx_band_mask = np.logical_or(idx_band_mask, get_band_mask(solar_wvl, e_b_mask))
        b_suffix = ''.join([str(o_b) for o_b in observe_bands])
    else:
        idx_band_mask = get_band_mask(solar_wvl, evaluate_band)
        b_suffix = str(evaluate_band)
    # generate mask of pixels used in comparison
    idx_lines_mask = get_linelist_mask(solar_wvl)
    # print number of pixels
    print ' Band mask pixels:', np.sum(idx_band_mask)
    print ' Initial linelist mask pixels', np.sum(idx_lines_mask)

    sim_metrices = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski',
                    'wchebyshev', 'sqeuclidean', 'euclidean', 'chi2', 'EW', 'median_sep', 'sum', 'px_over', 'px_under']
    sim_metrices_std = [m + '_std' for m in sim_metrices]
    sim_dtypes = ['float64' for i in range(2 * len(sim_metrices))]

    list_fit_results = list([])

    observe_masks = [len(idx_lines_mask)]  # whole range
    # observe_masks = np.linspace(0, len(idx_lines_mask), 10)[1:]  # multiple masks

    spectrum_flux_offsets = list([])
    for offset_percent in observe_flux:
        print ' Flux offset:', offset_percent
        filename_suffix = '_b'+b_suffix+'_flux{:.2f}'.format(offset_percent)
        flux_offset = (1. - solar_flx) * offset_percent
        spectrum_flux_offsets.append(flux_offset)

        file_out_fits = 'solar_similarity_narrow'+filename_suffix+'.fits'
        fit_results_path = 'metrices_snr_function'+filename_suffix+'.csv'
        list_fit_results.append(fit_results_path)

        sim_results = Table(names=np.hstack(('snr', 'snr_guesslike', sim_metrices, sim_metrices_std, 'n_px_eval')),
                            dtype=(np.hstack(('float64', 'float64', sim_dtypes, 'int64'))))

        # generate noise and evaluate distances
        if not isfile(file_out_fits):
            for l_mask in observe_masks:
                idx_band_lines_mask = np.logical_and(idx_lines_mask, idx_band_mask)
                print 'Number of abs line features:', np.sum(idx_band_lines_mask)
                # idx_band_lines_mask[np.int64(l_mask):] = False
                print '  Lines mask pixels for this band:', l_mask, np.sum(idx_band_lines_mask)
                for snr_spectrum in observe_snr:
                    print '  SNR generated:', snr_spectrum

                    if n_noise_samples == 1:
                        snr_noise_pred = np.zeros((n_noise_samples, len(solar_wvl)))
                    else:
                        snr_sigma = 1.0/snr_spectrum
                        if not gauss_noise:
                            # Poissonian noise
                            snr_noise_pred = np.random.poisson((1.0 / snr_sigma)**2, size=(n_noise_samples, len(solar_wvl)))
                            snr_noise_pred = snr_noise_pred / ((1.0 / snr_sigma)**2) - 1.
                        else:
                            # Gaussian noise
                            snr_noise_pred = np.random.normal(loc=0, scale=snr_sigma, size=(n_noise_samples, len(solar_wvl)))

                    if compute_guesslike_snr:
                        if guesslike_all_lines:
                            # median signal at all wavelength pixels
                            signal_med = np.nanmedian(solar_flx)
                            # determine actual snr of generated noise at selected pixels - guess like
                            signal_noisy = (solar_flx + snr_noise_pred)
                        else:
                            # median signal at selected abundance wavelength pixels
                            signal_med = np.nanmedian(solar_flx[idx_band_lines_mask])
                            # determine actual snr of generated noise at selected pixels - guess like
                            signal_noisy = (solar_flx[idx_band_lines_mask] + snr_noise_pred)[:, idx_band_lines_mask]
                        noise_med = 1.4826 / np.sqrt(2) * np.nanmedian(np.abs(signal_noisy[:, 1:] - signal_noisy[:, :-1]))
                        # plt.plot(signal_noisy[:, :-1][0])
                        # plt.plot(signal_noisy[:, 1:][0])
                        # plt.show()
                        snr_guess_like = signal_med / noise_med
                        print '   "Guess like" snr {:.2f}:'.format(snr_guess_like)
                    else:
                        snr_guess_like = 0.

                    # iterate and add noise to observed spectrum
                    spectrum_distances = np.zeros((n_noise_samples, len(sim_metrices)))

                    for i_snr in range(n_noise_samples):
                        # determine weights for the distance computation (different approaches)
                        solar_flx_noisy = snr_noise_pred[i_snr, :] + solar_flx + flux_offset
                        if normalize_noisy_spectrum:
                            norm_ok_mask = determine_norm_mask(solar_wvl[idx_band_mask], norm_bad_ranges)
                            solar_flx_noisy[idx_band_mask] = spectra_normalize(solar_wvl[idx_band_mask] - np.mean(solar_wvl[idx_band_mask]),
                                                                               solar_flx_noisy[idx_band_mask], fit_mask=norm_ok_mask,
                                                                               steps=15, sigma_low=2.5, sigma_high=2.5, order=1, n_min_perc=5.,
                                                                               return_fit=False, func='poly')

                        spectrum_distances[i_snr, :] = compute_distances(solar_flx[idx_band_lines_mask], None,
                                                                         solar_flx_noisy[idx_band_lines_mask],
                                                                         d=noise_power)

                    # add agregated results to final table
                    sim_results.add_row(np.hstack([snr_spectrum, snr_guess_like, np.nanmean(spectrum_distances, axis=0),
                                                   np.nanstd(spectrum_distances, axis=0), np.sum(idx_band_lines_mask)]))

                # save results
                sim_results.write(file_out_fits, overwrite=True)
        else:
            sim_results = Table.read(file_out_fits)

        # generate plots and fitted functions
        if n_noise_samples != 1:
            txt_out = open(fit_results_path, 'w')
            txt_out.write('metric,x_shift,pow_multi,pow_amp,pow_x_0,pow_alpha,lin_slop,lin_inter\n')
            for col in sim_metrices:

                # determine initial fit parameters
                sim_data = sim_results[col]
                if sim_data[0] > sim_data[-1]:
                    pow_multi = 1.
                else:
                    pow_multi = -1.
                sim_med = np.median(sim_data)
                if sim_med > 100:
                    pow_amp = 250.
                else:
                    pow_amp = 1.

                # fit function
                f_init = models.Shift(offset=0.) | (models.Const1D(amplitude=pow_multi)*models.PowerLaw1D(amplitude=pow_amp, x_0=1.0, alpha=1.0) + models.Linear1D(slope=0, intercept=sim_med))
                fitter = fitting.LevMarLSQFitter()
                gg_fit = fitter(f_init, sim_results['snr_guesslike'], sim_results[col])
                plt.plot(observe_snr, gg_fit(observe_snr), alpha=0.5, label='Fitted curve', c='black', lw=0.5)
                # plot
                plt.errorbar(sim_results['snr_guesslike'], sim_results[col], yerr=sim_results[col+'_std'], c='black', alpha=0.75, errorevery=2,
                             capsize=0, elinewidth=0.75, linewidth=0.5, fmt='o', ms=1.5, label='Simulations')  #, lw=0, s=3, c='black', alpha=0.5)
                # output(s)
                plt.xlabel('Simulated SNR')
                plt.ylabel('Similarity measure')
                plt.legend()
                plt.tight_layout()
                plt.savefig(col + filename_suffix+'.png', dpi=500)
                plt.close()
                # print gg_fit
                # save fitted values
                txt_out.write(col+','+str(gg_fit.offset_0.value)+','+str(gg_fit.amplitude_1.value)+','+str(gg_fit.amplitude_2.value)+','+str(gg_fit.x_0_2.value)+','+str(gg_fit.alpha_2.value)+','+str(gg_fit.slope_3.value)+','+str(gg_fit.intercept_3.value)+'\n')
            txt_out.close()

    # plot observed spectra
    plt.figure(figsize=(14, 4))
    for i_f in range(len(observe_flux)):
        plt.plot(solar_wvl, solar_flx + spectrum_flux_offsets[i_f], label='{:.2f}'.format(observe_flux[i_f]), lw=0.5)
    plt.xlim((min_wvl[evaluate_band-1], max_wvl[evaluate_band-1]))
    # plt.ylim((0.6, 1.03))
    # plt.xlim((4795, 4815))
    plt.tight_layout()
    plt.legend()
    # plt.show()
    if not run_unlike_solar:
        plt.savefig('solar_spectra_b'+b_suffix+'.png', dpi=300)
    else:
        plt.savefig('ref_spectra'+unlike_ref_suffix+'_b'+b_suffix+'.png', dpi=300)
    plt.close()

    print ' Final fit results comparison'
    # plot all curves with different flux offsets
    for metric in sim_metrices:
        for i_f in range(len(list_fit_results)):
            y_vals = metric_by_snr(Table.read(list_fit_results[i_f]), metric, observe_snr)
            plt.plot(observe_snr, y_vals, lw=1, label='{:.2f}'.format(observe_flux[i_f]))
        plt.legend()
        plt.savefig(metric+'_all'+unlike_ref_suffix+'_b'+b_suffix+'.png', dpi=500)
        plt.close()

    # only one iteration is needed
    if all_bands_at_once:
        break


