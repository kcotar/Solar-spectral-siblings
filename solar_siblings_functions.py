from socket import gethostname
# PC hostname
pc_name = gethostname()
if pc_name != 'gigli':
    import matplotlib
    matplotlib.use('Agg')
import imp
from astropy.table import Table, join
import matplotlib.pyplot as plt
import george, emcee, corner
from george import kernels
from george.modeling import Model
import pandas as pd
import numpy as np
from time import time
from scipy.spatial.distance import *
np.seterr(all='ignore')
from scipy.signal import correlate
from lmfit.models import GaussianModel, VoigtModel

OK_LINES_ONLY = True
USE_SUBSAMPLE = False
REF_SPECTRUM_PREPARE = False

if REF_SPECTRUM_PREPARE:
    # use for reference spectrum preparation
    min_wvl = np.array([4705, 5640, 6470, 7680])
    max_wvl = np.array([4915, 5885, 6750, 7900])
else:
    # use for normal processing and comparison
    min_wvl = np.array([4725, 5665, 6485, 7700])
    max_wvl = np.array([4895, 5865, 6725, 7875])

rv_weights = np.array([0.75, 1., 1., 0.25])

if USE_SUBSAMPLE:
    every_nth_solar_pixel = np.array([4, 5, 6, 7])  # double sub-sampling of spectra
else:
    every_nth_solar_pixel = np.array([8, 10, 12, 14])  # almost original sampling of the data - GP works much faster with this

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr52_dir = '/media/storage/HERMES_REDUCED/dr5.3/'
    galah_data_input = '/home/klemen/data4_mount/'
    out_dir = ''
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
else:
    dr52_dir = '/data4/cotar/dr5.3/'
    out_dir = '/data4/cotar/'
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *

galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')
# define lines with clearly visible telluric contamination between different night
problematic_lines = ['Ti5689.46',
                     'Si5690.425',
                     'Fe5696.0892',
                     'Sc5717.307',
                     'Cr5719.815',
                     'Ti5720.4359',
                     'Sc5724.107',
                     'Cr5787.919',
                     'Ba5853.668',
                     'Ti5866.4513',
                     'Ca5867.562',
                     'Ni6482.7983',
                     'Ca6508.8496',
                     'Fe6516.0766',
                     'Fe6518.3657',
                     'Ni6532.873',
                     'Ti6599.108'
                     ]  # TODO: might even add absorption lines in the beginning of the red arm - lines are too shallow there

# define regions without continuum or heavily influenced by residual telluric lines
# prepared specifically for Sun/flats and solar like spectra
norm_bad_ranges = [[4727.0, 4737.0],  # region without continuum
                   [4845.0, 4875.0],  # h-beta
                   [4902.0, 4905.0],  # useless part with problems
                   [5787.5, 5795.0],  # atmospheric absorption bands
                   [6492.5, 6498.0],  # region without continuum
                   [6542.0, 6575.0],  # h-alpha + nearby region without continuum
                   [7590.0, 7670.0]   # atmospheric absorption bands
                   ]

# filter read linelist
if OK_LINES_ONLY:
    remove_linelist_rows = list([])
    for i_l_r in range(len(galah_linelist)):
        line = galah_linelist[i_l_r]
        line_name = line['Element'] + str(line['line_centre'])
        if line_name in problematic_lines:
            remove_linelist_rows.append(i_l_r)
    galah_linelist.remove_rows(remove_linelist_rows)


def get_solar_data(solar_input_dir, suffix, every_nth=1, bands_together=True):
    if bands_together:
        solar_all = pd.read_csv(solar_input_dir + 'twilight_spectrum_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
        solar_wvl = solar_all[:, 0]
        solar_flx = solar_all[:, 1]
    else:
        solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
        solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
        solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
        solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
        solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
        solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))
    idx_finite = np.isfinite(solar_flx)
    if every_nth <= 1:
        return solar_wvl[idx_finite], solar_flx[idx_finite]
    else:
        return solar_wvl[idx_finite][::every_nth], solar_flx[idx_finite][::every_nth]


def get_linelist_mask(wvl_values, d_wvl=0):
    idx_lines_mask = wvl_values < 0.
    for line in galah_linelist:
        idx_lines_mask[np.logical_and(wvl_values >= line['line_start'] - d_wvl, wvl_values <= line['line_end'] + d_wvl)] = True
    return idx_lines_mask


def get_band_mask(wvl_values, evaluate_band):
    return np.logical_and(wvl_values >= min_wvl[evaluate_band - 1], wvl_values <= max_wvl[evaluate_band - 1])


def kernel_params_ok(p):
    amp, rad, amp2, rad2, amp_off, cont_norm = p
    if not 1e-9 < amp < 1e-2:
        return False
    if not 1e-8 < rad < 0.05:
        return False
    if not 1e-9 < amp2 < 1e-3:
        return False
    if not 0.1 < rad2 < 40:
        return False
    if not -0.2 < amp_off < 0.2:
        return False
    if not 0.8 < cont_norm < 1.2:
        return False
    return True


def kernel_noise(amp, rad):
    return amp * kernels.ExpSquaredKernel(rad)


def kernel_cont(amp, rad):
    return amp * kernels.Matern52Kernel(rad)


def spectrum_offset_norm(params, f):
    amp_off, cont_norm = params
    f_new = f + (1. - f) * amp_off
    f_new /= cont_norm
    return f_new


def get_kernel(p, add_cont=True):
    amp, rad, amp2, rad2 = p
    kernel = kernel_noise(amp, rad)
    if add_cont:
        kernel += kernel_cont(amp2, rad2)
    return kernel


def lnprob_gp(params, f_ref, f_obs, wvl, data_std, spectrum_off_norm=True):
    # evaluate selected parameters
    if kernel_params_ok(params):
        gp = george.GP(get_kernel(params[:-2]))
        if data_std is not None:
            gp.compute(wvl, data_std)
        else:
            gp.compute(wvl)
        # draw a sample noise from kernel distribution - not needed
        # gp_noise = gp.sample(size=1)
        # normalization and abs lines amplitude fitting
        if spectrum_off_norm:
            f_ref_new = spectrum_offset_norm(params[-2:], f_ref)
            return gp.lnlikelihood(f_ref_new - f_obs, wvl)
        else:
            return gp.lnlikelihood(f_ref - f_obs, wvl)
    else:
        return -np.inf


def fit_gp_kernel(init_guess, ref_data, obs_data, wvl, nwalkers=32, n_threds=1, n_burn=75, data_std=None,
                  n_per_burn=10, exit_lnp=1.5, normal_dist_guess=True):
    # compute number of burn steps
    n_burn_steps = np.int32(n_burn/n_per_burn)
    # data to be fitted handling
    ndim = len(init_guess)
    data_diff = ref_data - obs_data

    given_guess = np.array(init_guess)
    # add random amount of noise to the data
    if normal_dist_guess:
        # standard normal distribution of initial values
        perc_rand = 25.
        p0 = [given_guess + given_guess * np.random.randn(ndim) * perc_rand/100. for i_w in range(nwalkers)]
    else:
        # uniform distribution of initial values
        perc_rand = 100.
        p0 = list([])
        for i_w in range(nwalkers):
            guess_max_offset = given_guess * perc_rand/100.
            p0_new = given_guess + guess_max_offset * np.random.rand(ndim) - guess_max_offset / 2.
            # correction for spectrum amp and offset levels
            p0_new[-2] = 0. + np.random.rand(1) * 0.2 - 0.1
            p0_new[-1] = 1. + np.random.rand(1) * 0.2 - 0.1
            p0.append(p0_new)


    # initialize emcee sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, threads=n_threds, args=(ref_data, obs_data, wvl, data_std))

    print(' Running burn-in')
    time_1 = time()
    for i_b in range(n_burn_steps):
        print '  Run:', i_b+1
        if i_b == 0:
            p0, lnp, _ = sampler.run_mcmc(p0, n_per_burn)
        else:
            p0, lnp, _ = sampler.run_mcmc(None, n_per_burn)
        # test exit conditions
        if (lnp/len(data_diff) > exit_lnp).all():
            break

    time_2 = time()
    print '  {:.1f} min'.format((time_2-time_1)/60.)

    # print(' Running production')
    # p0 = [p + 1e-8 * np.random.randn(ndim) for i_w in range(nwalkers)]
    # p0, ln
    # p, _ = sampler.run_mcmc(p0, 500)
    # time_3 = time()
    # print '  {:.1f} min'.format((time_3 - time_2) / 60.)

    if n_threds > 1:
        sampler.pool.close()
    return sampler, p0, lnp


def compute_distances(obs, obs_std, orig, d=1., norm=True):
    # dist_weight = 1. / np.std(pix_ref - pix_spec)
    if d <= 0 or obs_std is None:
        dist_weight = np.ones(len(obs))
    else:
        dist_weight = (1. / obs_std) ** d

    spec_diff = np.abs(obs - orig)
    if norm:
        n_data = len(spec_diff)
    else:
        n_data = 1.
    results_list = [braycurtis(obs, orig, w=dist_weight),  # ze normaliziran sam po sebi
                    canberra(obs, orig, w=dist_weight)/n_data,
                    np.max(spec_diff),  # chebyshev(pix_ref, pix_spec),  # primerja samo element, ni potrebno normalizirat
                    np.sum(spec_diff * dist_weight)/n_data,  # cityblock(pix_ref, pix_spec, w=dist_weight),
                    correlation(obs, orig, w=dist_weight),
                    cosine(obs, orig, w=dist_weight),  # meri samo kot, ni potrebe po normalizaciji
                    np.sum((spec_diff**3. * dist_weight)**1./3.)/n_data,  # minkowski(pix_ref, pix_spec, 3., w=dist_weight),
                    np.max(spec_diff * dist_weight),  # weighted chebyshev  # spet se primerja samo en element od vseh
                    sqeuclidean(obs, orig, w=dist_weight)/n_data,  #
                    np.sum(np.sqrt(spec_diff**2. * dist_weight))/n_data,  # euclidean(pix_ref, pix_spec, w=dist_weight),
                    np.sum(spec_diff**2. * dist_weight)/n_data,
                    np.abs(np.sum(1. - obs) - np.sum(1. - orig))/n_data,
                    np.nanmedian(orig - obs)]  # difference in median levels, aka median spectral difference
    return results_list


def compute_ews(obs, wvl, orig, elements, linelist, d_wvl=0.0):
    results_list = list([])
    for elem in elements:
        linelist_sub = linelist[linelist['Element'] == elem]
        idx_lines_mask = wvl < 0.
        for line in linelist_sub:
            idx_lines_mask[np.logical_and(wvl >= line['line_start'] - d_wvl,
                                          wvl <= line['line_end'] + d_wvl)] = True
        if np.sum(idx_lines_mask) == 0:
            results_list.append(0.)
        else:
            ew_obs = np.sum(1 - obs[idx_lines_mask])
            ew_orig = np.sum(1 - orig[idx_lines_mask])
            results_list.append(100.*np.abs(ew_obs - ew_orig)/ew_orig)
    return results_list


def evaluate_spectrum(pix_spec, pix_std):
    if np.sum(pix_std < 0) > 0:
        print ' -> Probably bad spectrum'
        return False
    elif np.sum(~np.isfinite(pix_spec)) > 0:
        print ' -> Portion of spectra lies outside bounds'
        return False
    else:
        return True


def metric_by_snr(data, metric, snr_val):
    func_data = data[data['metric'] == metric]
    if len(func_data) == 1:
        return func_data['amplitude'] * (snr_val/func_data['x_0']) ** (-1*func_data['alpha']) + func_data['y_const']
    else:
        return np.full(len(snr_val), np.nan)


def plot_spectra(s_ids, spectra_dir_2, g_data, solar_flx, solar_wvl,
                 galah_linelist=None, save_path='plot.png', band=2, ext=0, y_range=None, x_range=None, d_wvl=0.):
    plt.figure(figsize=(13, 7))
    for s_id in s_ids:
        s2, w2 = get_spectra_dr52(str(s_id), bands=[band], root=spectra_dir_2, individual=False, extension=ext)
        if len(s2) == 0:
            continue
        if ext == 4:
            plt.plot(w2[0], s2[0], alpha=0.5, lw=1)
        if ext == 0:
            try:
                s2[0] = spectra_normalize(w2[0], s2[0], steps=5, sigma_low=1.5, sigma_high=2.5, order=1,
                                          n_min_perc=5., return_fit=False, func='poly')
                s2[0] = spectra_normalize(w2[0], s2[0], steps=25, sigma_low=1.8, sigma_high=3., order=17,
                                          n_min_perc=5., return_fit=False, func='poly')
                rv_shift = g_data[g_data['sobject_id'] == s_id]['rv_guess_shift']
                w2 *= (1 - (rv_shift) / 299792.458)
            except:
                continue
            plt.plot(w2[0], s2[0], alpha=0.2, lw=0.3)
    for line in galah_linelist:
        if line['line_centre'] < w2[0][-1] and line['line_centre'] > w2[0][0]:
            plt.axvspan(line['line_start'] - d_wvl, line['line_end'] + d_wvl, lw=0, color='black', alpha=0.2)
            plt.text(line['line_centre'], 1.05, line['Element'])
    plt.plot(solar_wvl, solar_flx, lw=1, c='black', alpha=0.5)
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


def plot_spectra_with_difference(flux1, flux2, wvl, flux3=None, diff_func=None,
                                 x_range=None, linelist=None, path=None, title=None):
    plt.figure(1, figsize=(12, 7))
    # [left, bottom, width, height]
    axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
    axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
    axSpectra.plot(wvl, flux1, c='black', lw=1.)
    axSpectra.plot(wvl, flux2, c='blue', lw=1.)
    if flux3 is not None:
        axSpectra.plot(wvl, flux3, c='red', lw=1.)
    axDiff.axhline(y=0, c='black', lw=1)
    axDiff.plot(wvl, flux1-flux2, c='blue', lw=1.)
    if diff_func is not None:
        axDiff.plot(wvl, diff_func, c='red', lw=1)
    axSpectra.set(ylim=(0.3, 1.1), xlim=x_range)
    axDiff.set(ylim=(-0.04, 0.04), xlim=x_range)
    if linelist is not None:
        d_abs_wvl = 0.0
        for line in linelist:
            if line['line_centre'] < wvl[-1] and line['line_centre'] > wvl[0]:
                axSpectra.axvspan(line['line_start'] - d_abs_wvl, line['line_end'] + d_abs_wvl, lw=0, color='black', alpha=0.2)
                axDiff.axvspan(line['line_start'] - d_abs_wvl, line['line_end'] + d_abs_wvl, lw=0, color='black', alpha=0.2)
    if title is not None:
        axSpectra.set_title(title)
    if path is not None:
        plt.savefig(path, dpi=400)
    else:
        plt.show()
    plt.close()


def get_wvl_log_shift(obs_flux, obs_wvl, ref_flux, ref_wvl, reduce_bands):
    log_shifts_res = list([])
    weights_res = list([])
    for i_r in range(len(reduce_bands)):
        ccd = reduce_bands[i_r] - 1
        idx_ref_sub = np.logical_and(ref_wvl >= min_wvl[ccd], ref_wvl <= max_wvl[ccd])
        ref_flux_sub = ref_flux[idx_ref_sub]
        ref_wvl_sub = ref_wvl[idx_ref_sub]
        wvl_step = ref_wvl_sub[1] - ref_wvl_sub[0]

        obs_flux_res = spectra_resample(obs_flux[i_r], obs_wvl[i_r], ref_wvl_sub, k=1)

        # correlation and stuff
        # get a valid subset of data
        idx_valid = np.isfinite(obs_flux_res)
        obs_flux_res = obs_flux_res[idx_valid]
        ref_flux_sub = ref_flux_sub[idx_valid]
        ref_wvl_sub = ref_wvl_sub[idx_valid]

        # convert spectra sampling to logspace
        obs_flux_res_log, wvl_valid_log = spectra_logspace(obs_flux_res, ref_wvl_sub)
        ref_flux_sub_log, _ = spectra_logspace(ref_flux_sub, ref_wvl_sub)

        # correlate the two spectra
        min_flux = 0.95
        ref_flux_sub_log[ref_flux_sub_log > min_flux] = 0.
        obs_flux_res_log[obs_flux_res_log > min_flux] = 0.
        corr_res = correlate(ref_flux_sub_log, obs_flux_res_log, mode='same', method='fft')

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

        # store to the array
        log_shifts_res.append(log_shift_wvl)
        weights_res.append(rv_weights[ccd])

    # compute weighted mean of results
    shift_res = np.sum(np.array(log_shifts_res)*np.array(weights_res))/np.sum(weights_res)
    print ' Log shifts:', log_shifts_res, shift_res
    return shift_res


def do_wvl_shift_log_and_resample(obs_flux, obs_wvl, log_shift, ref_wvl, reduce_bands):
    for i_r in range(len(reduce_bands)):
        ccd = reduce_bands[i_r] - 1
        idx_ref_sub = np.logical_and(ref_wvl >= min_wvl[ccd], ref_wvl <= max_wvl[ccd])
        ref_wvl_sub = ref_wvl[idx_ref_sub]

        obs_flux_log, obs_wvl_log = spectra_logspace(obs_flux[i_r], obs_wvl[i_r])
        obs_wvl_log -= log_shift
        obs_flux_res = spectra_resample(obs_flux_log, obs_wvl_log, ref_wvl_sub, k=1)
        # save results
        obs_flux[i_r] = obs_flux_res
        obs_wvl[i_r] = ref_wvl_sub
    # return results back
    return obs_flux, obs_wvl


def determine_norm_mask(in_wvl, ranges):
    mask = np.isfinite(in_wvl)
    for w_min, w_max in ranges:
        idx_maskout = np.logical_and(in_wvl >= w_min, in_wvl <= w_max)
        if np.sum(idx_maskout) > 0:
            mask[idx_maskout] = False
    return mask


def get_wvl_regions(ref_wvl, reduce_bands):
    wvl_list_out = list([])
    for i_r in range(len(reduce_bands)):
        ccd = reduce_bands[i_r] - 1
        idx_ref_sub = np.logical_and(ref_wvl >= min_wvl[ccd], ref_wvl <= max_wvl[ccd])
        wvl_list_out.append(ref_wvl[idx_ref_sub])
    return np.hstack(wvl_list_out)
