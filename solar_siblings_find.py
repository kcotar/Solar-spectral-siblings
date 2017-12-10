import os, imp
import george, emcee, corner
from george import kernels
from astropy.table import Table
from socket import gethostname
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

# -----------------------------------
# --------- Functions ---------------
# -----------------------------------


def kernel_params_ok(p):
    amp, rad, amp2, rad2 = p
    if not 1e-8 < amp < 1e-3:
        return False
    if not 1e-8 < rad < 0.02:
        return False
    if not 1e-8 < amp2 < 1e-4:
        return False
    if not 1 < rad2 < 50:
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


def lnprob_gp(params, data, wvl):
    # evaluate selected parameters
    if kernel_params_ok(params):
        gp = george.GP(get_kernel(params))
        gp.compute(wvl)
        return gp.lnlikelihood(data, wvl)
    else:
        return -np.inf


def fit_gp_kernel(init_guess, data, wvl, nwalkers=32, n_threds=1, n_burn=75):
    ndim = len(init_guess)

    given_guess = np.array(init_guess)
    # add random amount of noise to the data
    p0 = [given_guess + 1e-4 * np.random.randn(ndim) for i_w in range(nwalkers)]
    # multiply by the random amount of noise and add to the data - better for parameters of unequal values
    perc_rand = 20
    p0 = [given_guess + given_guess * np.random.randn(ndim) * perc_rand/100. for i_w in range(nwalkers)]  #

    # initialize emcee sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, threads=n_threds, args=(data, wvl))

    print(' Running burn-in')
    time_1 = time()
    p0, lnp, _ = sampler.run_mcmc(p0, n_burn)
    p = p0[np.argmax(lnp)]
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

# PC hostname
pc_name = gethostname()

# input data
if pc_name == 'gigli' or pc_name == 'klemen-P5K-E':
    dr52_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
    galah_data_input = '/home/klemen/GALAH_data/'
    imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
    imp.load_source('spectra_collection_functions', '../Carbon-Spectra/spectra_collection_functions.py')
    imp.load_source('distances', '../tSNE_test/distances.py')
else:
    galah_data_input = '/data4/cotar/'
from helper_functions import *
from spectra_collection_functions import *
from distances import *

# -----------------------------------
# --------- Settings ----------------
# -----------------------------------
d_wvl = 0.1
n_gp_samples = 150
save_plots = True
min_wvl = list([4730, 5670, 6490, 7705])
max_wvl = list([4890, 5860, 6720, 7870])

# reference solar spectra
print 'Read reference GALAH Solar spectra'
suffix = '_ext0_2_offset'
solar_input_dir = galah_data_input+'Solar_data/'
solar_g1 = pd.read_csv(solar_input_dir + 'b1_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(solar_input_dir + 'b2_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(solar_input_dir + 'b3_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(solar_input_dir + 'b4_solar_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_wvl = np.hstack((solar_g1[:, 0], solar_g2[:, 0], solar_g3[:, 0], solar_g4[:, 0]))
solar_flx = np.hstack((solar_g1[:, 1], solar_g2[:, 1], solar_g3[:, 1], solar_g4[:, 1]))

# downscale Solar spectra for faster processing
every_nth_pixel = 6
solar_wvl = solar_wvl[::every_nth_pixel]
solar_flx = solar_flx[::every_nth_pixel]

# data-table settings
data_date = '20171111'
galah_param_file = 'sobject_iraf_52_reduced_'+data_date+'.fits'
cannon_param_file = 'sobject_iraf_iDR2_171103_cannon.fits'

# select ok objects
print 'Reading and determining usable twilight flats for Solar statistics determination'
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist.csv')
galah_param = Table.read(galah_data_input + galah_param_file)
cannon_param = Table.read(galah_data_input + cannon_param_file)
idx_rows = np.logical_and(galah_param['red_flag'] == 64, galah_param['snr_c1_iraf'] > 50)
idx_rows = np.logical_and(idx_rows, galah_param['flag_guess'] == 0)
idx_rows = np.logical_and(idx_rows, galah_param['sobject_id'] > 140301000000000)

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
idx_solar_like = np.logical_and(idx_solar_like, cannon_param['snr_c2_iraf'] > 25)
n_solar_like = np.sum(idx_solar_like)
print 'Solar like by parameters:', n_solar_like

# -----------------------------------
# --------- Main program ------------
# -----------------------------------

solar_like_sobjects = galah_param['sobject_id'][idx_solar_like]
sim_results = Table(names=('sobject_id', 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'mahalanobis', 'minkowski','seuclidean','sqeuclidean','wminkowski','hamming'),
                    dtype=('int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64'))

file_out_fits = 'solar_similarity.fits'

move_to_dir('MultiDistance_solar_spectra')

# solar_like_sobjects = [160426005501042,170112001601216,170801004001010,170509004701028,171027002801397,170909002601089,160916004301242,160513002601086,161106002601077,150211004701110,160402004101021,151219003601245,150901000601193,170509005201063,161115002701019,161107001601133,150103004501206,170509004701096,160108003601136,150427004801275,160327006101097,170602003701147,170517001801299,160401004401064,170906003601002,170206005701110,160513002101209,161118004701221,170107004801249,140209001701187,140806002301191,160328003201333,150411006101130,160520002601397,170602006201184,160524006101090,170514002401099,160817002601198,161118004001114,170510004801244,170515006101072,160424004201234,140415002401342,170907002601114,161007002801397,161213004101187,171102005001020,151111002101170,160813003601074,170510001801331,170121002201177,170601003101229,170911002101145,161104002301194,160530005501014,161106003101031,150830005101337,170313001601013,170910004101092,170601003101179,170109002801065,140808004701129,160327003601122,150411004101331,150211004701179,170711004001148,160919001601330,160815004301329,160401004401168,150409004601199,140607000701025,150103003501108,170510002301079,161105003101025,170122002101017,140409003001378,141103003601053,160520003101160,151110002601078,170512000101244,140711002901298,170516003601221,161007003301369,140807005001217,150103004001346,170724004601092,160402004101104,170312001601218,141102003201127,150207005101248,170511001101086,170220003601186,160814000101228,150101004001039,150601004801369,170912002401246,140812003801356,170117002601213,170219003101017,160522006101021,170512000701198,150108001001167,170122002101290,161105003101238,160130005801117,170802002101225,160524004201389,170118001701346,160130004101262,150412003101062,150830005601079,150112002501171,140711003401047,170217004001356,150108002801008,161009002601357,170106004601076,151009003101098,170506004401238,150208002701189,150903002401220,150103004501019,161013003201073,140409003601076,170829001901077,161013003801334,150831004001118,140708000601373,140310003301075,170912001901323,150829004301108,140710006101115,160109003301105,150602004901396,150412002101315,160524004201163,160531002601078,160108002601310,171102003301055,160110003601106,150412003601054,160527002101312,170128003401173,140316002301154,140311009101067,140811005001213,131123003501359,140713004001103,170128002101119,140814004301124]
# solar_like_sobjects = solar_like_sobjects[np.int64(np.random.rand(50)*len(solar_like_sobjects))]

for s_obj in solar_like_sobjects:
    print 'Evaluating', s_obj
    # get spectra of all bands for observed objects
    read_ext = 0
    flux, wvl = get_spectra_dr52(str(s_obj), bands=[1, 2, 3, 4], root=dr52_dir, extension=read_ext, individual=False)
    if read_ext == 0:
        # normalize flux
        for i_c in range(4):
            # ------ NORM v1 - high order polynomial, many steps
            # flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=35, sigma_low=1.5, sigma_high=2.8, order=29,
            #                               n_min_perc=3.,  return_fit=False, func='poly')
            # ------ NORM v2 - the same as used in the process of reference Solar spectra construction
            flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=5, sigma_low=1.5, sigma_high=2.5, order=1, n_min_perc=5.,
                                          return_fit=False, func='poly')
            flux[i_c] = spectra_normalize(wvl[i_c], flux[i_c], steps=25, sigma_low=1.8, sigma_high=3., order=17,
                                          n_min_perc=5., return_fit=False, func='poly')
        # apply computed rv shift to the spectrum
        rv_shift = galah_param[galah_param['sobject_id'] == s_obj]['rv_guess_shift']
        wvl *= (1 - rv_shift / 299792.458)

    # # compute per band spectra similarity
    # spectra_similarity = np.zeros(n_gp_samples, dtype=np.float64)
    # for i_c in range(4):
    #     # define subset of spectra to be compared to reference solar spectrum
    #     idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
    #     flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
    #
    #     # determine spectrum difference and its variance
    #     diff = (solar_flx[idx_ref] - flux_b_res)
    #     diff_var = np.nanvar(diff)
    #
    #     # determine kernel parameters trough emcee fit
    #     print ' Running emcee'
    #     # emcee_fit_px = 100
    #     emcee_fit_px = len(diff)
    #     sampler, fit_res, fit_prob = fit_gp_kernel([diff_var, 0.008, 1e-5, 10],
    #                                                diff[:emcee_fit_px], solar_wvl[idx_ref][:emcee_fit_px],
    #                                                nwalkers=20, n_threds=20, n_burn=45)
    #
    #     # walker prob plot
    #     print("Plotting walker probabilities")
    #     walkers_prob = sampler.lnprobability
    #     for i_w in range(walkers_prob.shape[0]):
    #         plt.plot(walkers_prob[i_w, :])
    #     plt.savefig(str(s_obj) + '_gp-lnprob_b' + str(i_c + 1) + '.png', dpi=400)
    #     # plt.show()
    #     plt.close()
    #
    #     sampler_chain_vals = sampler.flatchain
    #     kernel_fit = np.median(sampler_chain_vals, axis=0)  # flatchain holds parameters of all emcee steps
    #     # kernel_fit = np.median(fit_res, axis=0)  # fit_res holds only the parameters of the last step
    #     # kernel_fit = fit_res[np.argmax(fit_prob)]
    #     print 'Median val:', np.median(sampler_chain_vals, axis=0)
    #     print 'Max lnprob:', fit_res[np.argmax(fit_prob)]
    #     print diff_var
    #     print kernel_fit
    #
    #     # corner plot of parameters
    #     if save_plots:
    #         c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
    #                               labels=['amp_noise', 'rad_noise', 'amp_cont', 'rad_cont'], bins=30)
    #         c_fig.savefig(str(s_obj)+'_corner_b'+str(i_c+1)+'.png', dpi=400)
    #         plt.close(c_fig)
    #
    #     # create a gaussian process that will be used for the whole spectra
    #     gp = george.GP(get_kernel(kernel_fit))
    #     gp.compute(solar_wvl[idx_ref])
    #     gp_noise_pred = gp.sample(size=n_gp_samples)
    #
    #     if save_plots:
    #         plt.plot(solar_flx[idx_ref], c='red', lw=0.5)
    #         for i_pred in range(15):
    #             plt.plot(solar_flx[idx_ref] + gp_noise_pred[i_pred, :], c='black', alpha=0.2, lw=0.3)
    #         plt.plot(flux_b_res, c='blue', lw=0.5)
    #         plt.ylim((0.4, 1.1))
    #         # plt.show()
    #         plt.savefig(str(s_obj)+'_gp_b'+str(i_c+1)+'.png', dpi=400)
    #         plt.close()
    #
    #     # filter out possible strange flux values that will cause anomalous distance estimation
    #     # idx_bad = np.logical_or(flux_b_res > 1.3, flux_b_res < 0.01)
    #     # idx_bad = flux_b_res > 1.3
    #     # flux_b_res[idx_bad] = np.nan
    #     flux_b_res[flux_b_res > 1.15] = 1.15
    #     flux_b_res[flux_b_res < 0.] = 0.
    #
    #     # fill possible missing values with ones
    #     idx_missing = ~np.isfinite(flux_b_res)
    #     if np.sum(idx_missing) > 0:
    #         print '  filled pixels', np.sum(idx_missing)
    #         flux_b_res[idx_missing] == 1.
    #
    #     # compute similarity for the whole bathc of data
    #     spectra_diff = (solar_flx[idx_ref] + gp_noise_pred) - flux_b_res
    #
    #     # mask difference by elements absorption lines and compute similarity/distance estimator
    #     abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
    #     spectra_eucl_dist = np.sum(spectra_diff[:, abs_lines_cols]**2, axis=1)
    #     print ' pixels evaluated:', len(abs_lines_cols)
    #     spectra_similarity += spectra_eucl_dist
    #
    # spectra_similarity = np.sqrt(spectra_similarity)
    # sim_results.add_row([s_obj, d.np.mean(spectra_similarity), np.std(spectra_similarity)])
    #
    # if save_plots:
    #     plt.hist(spectra_similarity, bins=20)
    #     plt.savefig(str(s_obj) + '_sim_hist.png', dpi=400)
    #     plt.close()

    pix_ref = list([])
    pix_spec = list([])
    for i_c in range(4):
        idx_ref = np.logical_and(solar_wvl >= min_wvl[i_c], solar_wvl <= max_wvl[i_c])
        abs_lines_cols = np.where(idx_lines_mask[idx_ref])[0]
        # print flux[i_c], wvl[i_c], solar_wvl[idx_ref]
        flux_b_res = spectra_resample(flux[i_c], wvl[i_c], solar_wvl[idx_ref], k=1)
        spectra_diff = solar_flx[idx_ref] - flux_b_res

        pix_ref.append(solar_flx[idx_ref][abs_lines_cols])
        pix_spec.append(flux_b_res[abs_lines_cols])

    # compute different distance meassurements
    pix_ref = np.hstack(pix_ref)
    pix_spec = np.hstack(pix_spec)

    sim_results.add_row([s_obj,
                         braycurtis(pix_ref, pix_spec),
                         canberra(pix_ref, pix_spec),
                         chebyshev(pix_ref, pix_spec),
                         cityblock(pix_ref, pix_spec),
                         correlation(pix_ref, pix_spec),
                         np.nan,  # mahalanobis(pix_ref, pix_spec),
                         minkowski(pix_ref, pix_spec, 1.),
                         np.nan,  # seuclidean(pix_ref, pix_spec),
                         sqeuclidean(pix_ref, pix_spec),
                         np.nan,  # wminkowski(pix_ref, pix_spec),
                         hamming(pix_ref, pix_spec)])

# check output file with results
if os.path.isfile(file_out_fits):
    os.remove(file_out_fits)
sim_results.write(file_out_fits)

print sim_results
print ''
sobj_id_like = sim_results[np.argsort(sim_results['dist_mean'])[:75]]['sobject_id']
print ','.join([str(s) for s in sobj_id_like])

print ''
sobj_id_dislike = sim_results[np.argsort(sim_results['dist_mean'])[-75:]]['sobject_id']
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

