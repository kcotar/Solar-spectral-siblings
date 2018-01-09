import matplotlib
matplotlib.use('Agg')
import os, imp
from astropy.table import Table, join
from socket import gethostname
import matplotlib.pyplot as plt
import george, emcee, corner
from george import kernels
import pandas as pd
import numpy as np
from time import time
from scipy.spatial.distance import *
# braycurtis(u, v[, w]) 	Computes the Bray-Curtis distance between two 1-D arrays.
# canberra(u, v[, w]) 	Computes the Canberra distance between two 1-D arrays.
# chebyshev(u, v) 	Computes the Chebyshev distance.
# cityblock(u, v[, w]) 	Computes the City Block (Manhattan) distance.
# correlation(u, v[, w, centered]) 	Computes the correlation distance between two 1-D arrays.
# cosine(u, v[, w]) 	Computes the Cosine distance between 1-D arrays.
# euclidean(u, v[, w]) 	Computes the Euclidean distance between two 1-D arrays.
# mahalanobis(u, v, VI) 	Computes the Mahalanobis distance between two 1-D arrays.
# minkowski(u, v[, p, w]) 	Computes the Minkowski distance between two 1-D arrays.
# seuclidean(u, v, V) 	Returns the standardized Euclidean distance between two 1-D arrays.
# sqeuclidean(u, v[, w]) 	Computes the squared Euclidean distance between two 1-D arrays.
# wminkowski(u, v, p, w) 	Computes the weighted Minkowski distance between two 1-D arrays.
#
# Distance functions between two boolean vectors (representing sets) u and v. As in the case of numerical vectors, pdist is more efficient for computing the distances between all pairs.
# dice(u, v[, w]) 	Computes the Dice dissimilarity between two boolean 1-D arrays.
# hamming(u, v[, w]) 	Computes the Hamming distance between two 1-D arrays.
# jaccard(u, v[, w]) 	Computes the Jaccard-Needham dissimilarity between two boolean 1-D arrays.
# kulsinski(u, v[, w]) 	Computes the Kulsinski dissimilarity between two boolean 1-D arrays.
# rogerstanimoto(u, v[, w]) 	Computes the Rogers-Tanimoto dissimilarity between two boolean 1-D arrays.
# russellrao(u, v[, w]) 	Computes the Russell-Rao dissimilarity between two boolean 1-D arrays.
# sokalmichener(u, v[, w]) 	Computes the Sokal-Michener dissimilarity between two boolean 1-D arrays.
# sokalsneath(u, v[, w]) 	Computes the Sokal-Sneath dissimilarity between two boolean 1-D arrays.
# yule(u, v[, w]) 	Computes the Yule dissimilarity between two boolean 1-D arrays.

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


def get_solar_data(solar_input_dir, suffix):
    solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
    solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
    solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
    solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah' + suffix + '.txt', header=None, delimiter=' ', na_values='nan').values
    solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
    solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))
    idx_finite = np.isfinite(solar_flx)
    return solar_wvl[idx_finite], solar_flx[idx_finite]


def kernel_params_ok(p):
    amp, rad, amp2, rad2 = p
    if not 1e-8 < amp < 1e-3:
        return False
    if not 1e-8 < rad < 0.02:
        return False
    if not 1e-8 < amp2 < 1e-4:
        return False
    if not 1 < rad2 < 35:
        return False
    return True


def kernel_noise(amp, rad):
    return amp * kernels.ExpSquaredKernel(rad)


def kernel_cont(amp, rad):
    return amp * kernels.Matern52Kernel(rad)


def get_kernel(p, add_cont=True):
    amp, rad, amp2, rad2 = p
    kernel = kernel_noise(amp, rad)
    if add_cont:
        kernel += kernel_cont(amp2, rad2)
    return kernel


def lnprob_gp(params, data, wvl, data_std):
    # evaluate selected parameters
    if kernel_params_ok(params):
        gp = george.GP(get_kernel(params))
        if data_std is not None:
            gp.compute(wvl, data_std)
        else:
            gp.compute(wvl)
        return gp.lnlikelihood(data, wvl)
    else:
        return -np.inf


def fit_gp_kernel(init_guess, data, wvl, nwalkers=32, n_threds=1, n_burn=75, data_std=None,
                  n_per_burn=10, exit_lnp=1.5):
    n_burn_steps = np.int32(n_burn/n_per_burn)
    ndim = len(init_guess)

    given_guess = np.array(init_guess)
    # add random amount of noise to the data
    # p0 = [given_guess + 1e-4 * np.random.randn(ndim) for i_w in range(nwalkers)]
    # multiply by the random amount of noise and add to the data - better for parameters of unequal values
    perc_rand = 20.
    p0 = [given_guess + given_guess * np.random.randn(ndim) * perc_rand/100. for i_w in range(nwalkers)]  #

    # initialize emcee sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, threads=n_threds, args=(data, wvl, data_std))

    print(' Running burn-in')
    time_1 = time()
    for i_b in range(n_burn_steps):
        print '  Run:', i_b+1
        if i_b == 0:
            p0, lnp, _ = sampler.run_mcmc(p0, n_per_burn)
        else:
            p0, lnp, _ = sampler.run_mcmc(None, n_per_burn)
        # test exit conditions
        if (lnp/len(data) > exit_lnp).all():
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
                    np.abs(np.sum(1. - obs) - np.sum(1. - orig))/n_data]
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
