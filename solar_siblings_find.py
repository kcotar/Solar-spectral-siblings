import os, imp
import george, emcee, corner
from george import kernels
from astropy.table import Table
from socket import gethostname
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# -----------------------------------
# --------- Functions ---------------
# -----------------------------------


def kernel_params_ok(p):
    amp, rad = p
    if not 1e-8 < amp < 1:
        return False
    if not 1e-8 < rad < 0.02:
        return False
    return True


def get_kernel(p):
    amp, rad = p
    return amp * kernels.Matern52Kernel(rad)
    # + kernels.Matern32Kernel(0.01)) # (kernels.ExpSquaredKernel(.5) + kernels.Matern32Kernel(.5) + kernels.RationalQuadraticKernel)


def lnprob_gp(params, data, wvl):
    # evaluate selected parameters
    if kernel_params_ok(params):
        gp = george.GP(get_kernel(params))
        gp.compute(wvl)
        return gp.lnlikelihood(data, wvl)
    else:
        return -np.inf


def fit_gp_kernel(init_guess, data, wvl, nwalkers=32, n_threds=1):
    ndim = len(init_guess)
    p0 = [np.array(init_guess) + 1e-4 * np.random.randn(ndim) for i_w in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, threads=n_threds, args=(data, wvl))

    print(' Running burn-in')
    time_1 = time()
    p0, lnp, _ = sampler.run_mcmc(p0, 100)
    p = p0[np.argmax(lnp)]
    time_2 = time()
    print '  {:.1f} min'.format((time_2-time_1)/60.)

    # print(' Running production')
    # p0 = [p + 1e-8 * np.random.randn(ndim) for i_w in range(nwalkers)]
    # p0, lnp, _ = sampler.run_mcmc(p0, 500)
    # time_3 = time()
    # print '  {:.1f} min'.format((time_3 - time_2) / 60.)

    if n_threds > 1:
        sampler.pool.close()
    return sampler, p0, lnp

# PC hostname
pc_name = gethostname()

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr52_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
    galah_data_input = '/home/klemen/GALAH_data/'
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
else:
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
n_gp_samples = 120
save_plots = True
min_wvl = list([4730, 5670, 6490, 7705])
max_wvl = list([4890, 5860, 6720, 7870])

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix = '_ext0'
solar_input_dir = galah_data_input+'Solar_data/'
solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))

# downscale Solar spectra for faster processing
solar_wvl = solar_wvl[::4]
solar_flx = solar_flx[::4]

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 50)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

# linelist mask
idx_lines_mask = solar_wvl < 0.
for line in galah_linelist:
    idx_lines_mask[np.logical_and(solar_wvl >= line['line_start']-0.1, solar_wvl <= line['line_end']+0.1)] = True
print 'Linelist mask pixels', np.sum(idx_lines_mask)

# find Solar parameters
teff_solar = np.nanmedian(galah_param[idx_rows]['teff_guess'])
teff_solar_std = np.nanstd(galah_param[idx_rows]['teff_guess'])
logg_solar = np.nanmedian(galah_param[idx_rows]['logg_guess'])
logg_solar_std = np.nanstd(galah_param[idx_rows]['logg_guess'])
feh_solar = np.nanmedian(galah_param[idx_rows]['feh_guess'])
print 'Solar parameters:', teff_solar, '+/-', teff_solar_std, ',  ', logg_solar, '+/-', logg_solar_std, ',  ', feh_solar

# Search for objects with similar physical properties
idx_solar_like = np.logical_and(np.abs(galah_param['teff_guess'] - teff_solar) < teff_solar_std,
                                np.abs(galah_param['logg_guess'] - logg_solar) < logg_solar_std)
idx_solar_like = np.logical_and(idx_solar_like, galah_param['red_flag'] == 0)
idx_solar_like = np.logical_and(idx_solar_like, galah_param['snr_c2_iraf'] > 25)
n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = galah_param['sobject_id'][idx_solar_like]
sim_results = Table(names=('sobject_id', 'dist_mean', 'dist_std'),
                    dtype=('int64', 'float64', 'float64'))

file_out_fits = 'solar_similarity_gp.fits'

move_to_dir('GaussianProcess_spectra')
for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
    # get spectra of all bands for observed objects
    read_ext = 0
    flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext, individual=False)
    if read_ext == 0:
        # normalize flux
        for i_c in range(4):
            # apply the same normalization as in the process of creation of master solar spectrum
            flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
                                          n_min_perc=3.,  return_fit=False, func='poly')
        # apply computed rv shift to the spectrum
        rv_shift = galah_param[galah_param['sobject_id'] == s_obj]['rv_guess_shift']
        wvl *= (1 - rv_shift / 299792.458)

    # compute per band spectra similarity
    spectra_similarity = np.zeros(n_gp_samples, dtype=np.float64)
    for i_c in range(4):
        # define subset of spectra to be compared to reference solar spectrum
        idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
        flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)

        # determine spectrum difference and its variance
        diff = (solar_flx[idx_ref] - flux_b_res)
        diff_var = np.nanvar(diff)

        # determine kernel parameters trough emcee fit
        print ' Running emcee'
        emcee_fit_px = 350
        sampler, fit_res, fit_prob = fit_gp_kernel([diff_var, 0.008],
                                                   diff[:emcee_fit_px], solar_wvl[idx_ref][:emcee_fit_px],
                                                   nwalkers=10, n_threds=20)

        kernel_fit = fit_res[np.argmax(fit_prob)]
        print diff_var
        print kernel_fit

        # corner plot of parameters
        if save_plots:
            c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, labels=['amp', 'rad'], bins=30)
            c_fig.savefig(str(s_obj)+'_corner_b'+str(i_c+1)+'.png', dpi=400)
            plt.close(c_fig)

        # create a gaussian process that will be used for the whole spectra
        gp = george.GP(get_kernel(kernel_fit))
        gp.compute(solar_wvl[idx_ref])
        gp_noise_pred = gp.sample(size=n_gp_samples)

        if save_plots:
            plt.plot(solar_flx[idx_ref], c='red', lw=0.5)
            for i_pred in range(15):
                plt.plot(solar_flx[idx_ref] + gp_noise_pred[i_pred, :], c='black', alpha=0.2, lw=0.3)
            plt.plot(flux_b_res, c='blue', lw=0.5)
            plt.ylim((0.4, 1.1))
            # plt.show()
            plt.savefig(str(s_obj)+'_gp_b'+str(i_c+1)+'.png', dpi=400)
            plt.close()

        # filter out possible strange flux values that will cause anomalous distance estimation
        # idx_bad = np.logical_or(flux_b_res > 1.3, flux_b_res < 0.01)
        # idx_bad = flux_b_res > 1.3
        # flux_b_res[idx_bad] = np.nan
        flux_b_res[flux_b_res > 1.15] = 1.15
        flux_b_res[flux_b_res < 0.] = 0.

        # fill possible missing values with ones
        idx_missing = ~np.isfinite(flux_b_res)
        if np.sum(idx_missing) > 0:
            print '  filled pixels', np.sum(idx_missing)
            flux_b_res[idx_missing] == 1.

        # compute similarity for the whole bathc of data
        spectra_diff = (solar_flx[idx_ref] + gp_noise_pred) - flux_b_res

        # mask difference by elements absorption lines and compute similarity/distance estimator
        abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
        spectra_eucl_dist = np.sum(spectra_diff[:, abs_lines_cols]**2, axis=1)
        print ' pixels evaluated:', len(abs_lines_cols)
        spectra_similarity += spectra_eucl_dist

    spectra_similarity = np.sqrt(spectra_similarity)
    sim_results.add_row([s_obj, np.mean(spectra_similarity), np.std(spectra_similarity)])

    # check results
    if os.path.isfile(file_out_fits):
        os.remove(file_out_fits)
    sim_results.write(file_out_fits)

    if save_plots:
        plt.hist(spectra_similarity, bins=20)
        plt.savefig(str(s_obj) + '_sim_hist.png', dpi=400)
        plt.close()

print sim_results
print ''
sobj_id_like = sim_results[np.argsort(sim_results['dist_mean'])[:50]]['sobject_id']
print ','.join([str(s) for s in sobj_id_like])

print ''
sobj_id_dislike = sim_results[np.argsort(sim_results['dist_mean'])[-50:]]['sobject_id']
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

