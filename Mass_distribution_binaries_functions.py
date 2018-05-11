import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imp
import numpy as np
import emcee, corner
import thecannon as tc
from astropy.table import Table, join
from copy import deepcopy
from lmfit.models import SkewedGaussianModel
from os import chdir, path
from time import time

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import *

# directories and data
dr53_dir = '/media/storage/HERMES_REDUCED/dr5.3/'
galah_data_input = '/home/klemen/data4_mount/'
isochrones_dir = galah_data_input+'isochrones/padova_Gaia_DR2/'
isochrones_data = Table.read(isochrones_dir + 'isochrones_all.fits')

# Canonn linelist
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')


# preload and prepare everything - drugace se pojavlja nek cuden error ce to delam znotraj sample procedure
# load cannon spectral model created by Gregor
cannon_model = tc.CannonModel.read('cannon_model_noflats_nobinaries.dat')
thetas = cannon_model.theta
vectorizer = cannon_model.vectorizer
fid = cannon_model._fiducials
sca = cannon_model._scales


def plot_lnprob(chain_vals, fit_sel_vals, path):
    c_fig = corner.corner(chain_vals, truths=fit_sel_vals,
                          quantiles=[0.16, 0.5, 0.84],
                          labels=['Teff 1', 'Teff 2'], bins=75)
    c_fig.savefig(path, dpi=200)
    plt.close(c_fig)


def plot_walkers(walkers_prob, path, title=''):
    for i_w in range(walkers_prob.shape[0]):
        plt.plot(walkers_prob[i_w, :], lw=0.3)
    walkers_prob = walkers_prob.flatten()  # without this correction numpy
    walkers_prob = walkers_prob[np.isfinite(walkers_prob)]  # percentile may return incorrect -inf value
    plt.title(title)
    plt.ylim((np.percentile(walkers_prob, 0.1), np.percentile(walkers_prob, 99.9)))
    plt.savefig(path, dpi=200)
    plt.close()


def add_mag(m1, m2):
    return -2.5 * np.log10(10**(-0.4*m1) + 10**(-0.4*m2))


def select_isochrone(iso_all, mh, age):
    meh_uniq = np.unique(iso_all['MHini'])
    isochrone_meh = meh_uniq[np.argmin(np.abs(meh_uniq - mh))]

    age_uniq = np.unique(iso_all['Age'])
    isochrone_age = age_uniq[np.argmin(np.abs(age_uniq - age))]

    iso_sub = iso_all[np.logical_and(iso_all['MHini'] == isochrone_meh,
                                     iso_all['Age'] == isochrone_age)]
    # return iso_sub[iso_sub['logg'] > 4]
    mass_loss = iso_sub['Mini'] - iso_sub['Mass']
    return iso_sub[mass_loss < 0.1]


def get_gmag_from_teff_MS(iso, teff):
    # determine teff upper and lower data point in isochrone
    for i_i in range(len(iso)):
        if iso[i_i]['teff'] > teff:
            break
    # determine Gmag for selected Teff
    teff_l = iso[i_i - 1]['teff']
    teff_h = iso[i_i]['teff']
    Gmag_l = iso[i_i - 1]['Gmag']
    Gmag_h = iso[i_i]['Gmag']
    k_g = (Gmag_h - Gmag_l)/(teff_h - teff_l)
    n_g = Gmag_l - k_g * teff_l
    Gmag_new = k_g * teff + n_g
    # print Gmag_h, Gmag_l, teff_h, teff_l, Gmag_new, teff
    return Gmag_new

    # iso_sub = iso[iso['Gmag'] > iso['Gmag'][np.argmax(iso['teff'])]]
    # idx_iso = np.argsort(np.abs(iso_sub['teff'] - teff))[:2]
    # # get point point between the nearest ones
    # k_g = (iso_sub['Gmag'][idx_iso[0]] - iso_sub['Gmag'][idx_iso[1]]) / (iso_sub['teff'][idx_iso[0]] - iso_sub['teff'][idx_iso[1]])
    # n_g = iso_sub['Gmag'][idx_iso[0]] - k_g * iso_sub['teff'][idx_iso[0]]
    # Gmag_new = k_g * teff + n_g
    # if not np.isfinite(Gmag_new):
    #     print iso_sub['Gmag'][idx_iso]
    #     print iso_sub['teff'][idx_iso]
    #     print teff, k_g, n_g
    # return Gmag_new


def get_mass_from_gmag(iso, gmag_abs=None, gbpmag_abs=None, grpmag_abs=None):
    idx_iso = np.argsort(np.abs(iso['Gmag'] - gmag_abs))[:2]
    # get point point between the nearest ones
    d_frac = (iso['Gmag'][idx_iso[0]] - gmag_abs) / (iso['Gmag'][idx_iso[0]] - iso['Gmag'][idx_iso[1]])
    mass = iso['Mass'][idx_iso[0]] - d_frac * (iso['Mass'][idx_iso[0]] - iso['Mass'][idx_iso[1]])
    return mass


def get_mass_ration(iso, teff1, teff2):
    Gmag_1 = get_gmag_from_teff_MS(iso, teff1)
    Gmag_2 = get_gmag_from_teff_MS(iso, teff2)
    return get_mass_from_gmag(iso, gmag_abs=Gmag_2)/get_mass_from_gmag(iso, gmag_abs=Gmag_1)


def eval_params(p, teff_range):
    t1, t2 = p
    if t2 > t1:
        return False
    if not teff_range[0] < t1 < teff_range[1]:
        return False
    if not teff_range[0] < t2 < teff_range[1]:
        return False
    return True


def lnprob_mag_fit(params, iso, G_mag_ref, G_mag_obj, teff_range):
    if eval_params(params, teff_range):
        teff_1, teff_2 = params
        Gmag_1 = get_gmag_from_teff_MS(iso, teff_1)
        Gmag_2 = get_gmag_from_teff_MS(iso, teff_2)
        Gmag_comb = add_mag(Gmag_1, Gmag_2)
        # Excess magnitudes
        Gmag_comb_excess = G_mag_ref - Gmag_comb
        Gmag_obj_excess = G_mag_ref - G_mag_obj
        # similarity/probability
        # return -1. * np.abs(Gmag_comb_excess - Gmag_obj_excess)
        return -1. * (Gmag_comb_excess - Gmag_obj_excess)**2
    else:
        return -np.inf


def fit_teff_values(iso, G_mag_ref, G_mag_obj, teff_init,
                    nwalkers=50, n_threds=10, n_steps=500):

    teff_range = [np.min(iso['teff']), np.max(iso['teff'])]
    teff_init_range = 150.

    p0 = list([])
    for i_w in range(nwalkers):
        p0_new = teff_init-teff_init_range/2. + np.random.rand(2) * teff_init_range
        p0.append(p0_new)

    # init samples
    sampler = emcee.EnsembleSampler(nwalkers, 2, lnprob_mag_fit, threads=n_threds,
                                    args=(iso, G_mag_ref, G_mag_obj, teff_range))
    p0, lnp, _ = sampler.run_mcmc(p0, n_steps)

    if n_threds > 1:
        sampler.pool.close()
    return sampler, p0, lnp


def get_cannon(teff, logg, feh):
    # get flux from model - initial implementation by Gregor
    # unpack values
    # thetas, vectorizer, fid, sca = cannon_model_data
    # thetas = cannon_model_data.theta
    # fid = cannon_model_data._fiducials
    # sca = cannon_model_data._scales
    # print '--------------------'
    # print thetas, fid, scac

    # run
    sint = thetas[:, 0] * 0.0
    labs = (np.array([teff, logg, feh]) - fid) / sca
    vec = vectorizer(labs)
    for i, j in enumerate(vec):
        sint += thetas[:, i] * j

    return sint


def get_linelist_mask(wvl_values, d_wvl=0., element=None):
    idx_lines_mask = wvl_values < 0.

    if element is None:
        galah_linelist_use = deepcopy(galah_linelist)
    else:
        galah_linelist_use = galah_linelist[galah_linelist['Element'] == element]

    for line in galah_linelist_use:
        idx_lines_mask[np.logical_and(wvl_values >= line['line_start'] - d_wvl, wvl_values <= line['line_end'] + d_wvl)] = True

    return idx_lines_mask


def get_spectra_complete(s_id):
    flx, wvl, sig = get_spectra_dr52(str(s_id), bands=[1, 2, 3], root=dr53_dir, extension=4, read_sigma=True)
    return np.hstack(flx), np.hstack(sig), np.hstack(wvl)


def lnprob_flx_mag_fit(params, iso, idx_wvl_model, flx_obs, flx_obs_s, G_mag_ref, G_mag_obj,
                       teff_range, logg_cannon, feh_cannon):

    if eval_params(params, teff_range):
        teff_1, teff_2 = params
        Gmag_1 = get_gmag_from_teff_MS(iso, teff_1)
        Gmag_2 = get_gmag_from_teff_MS(iso, teff_2)
        Gmag_comb = add_mag(Gmag_1, Gmag_2)

        # Excess magnitudes, fluxes, ratios,
        Gmag_comb_excess = G_mag_ref - Gmag_comb
        Gmag_obj_excess = G_mag_ref - G_mag_obj
        j_ratio = 10**(-0.4*Gmag_2) / 10**(-0.4*Gmag_1)

        # get spectra
        flx_1 = get_cannon(teff_1, logg_cannon, feh_cannon)[idx_wvl_model]
        flx_2 = get_cannon(teff_2, logg_cannon, feh_cannon)[idx_wvl_model]
        # combine them
        flx_model = 1./(1.+j_ratio) * flx_1 + 1./(1.+(1./j_ratio)) * flx_2
        # print teff_1, teff_2, j_ratio, 1./(1.+j_ratio), 1./(1.+(1./j_ratio))

        # plt.plot(flx_obs, label='obs')
        # plt.plot(flx_1, label='flx1')
        # plt.plot(flx_2, label='flx2')
        # plt.plot(flx_model, label='flx_sum')
        # plt.legend()
        # plt.show()
        # plt.close()

        # similarity/probability
        lnprob_mag = -5e3 * np.abs(Gmag_comb_excess - Gmag_obj_excess)  # linear difference
        # lnprob_mag = -250e3 * ((Gmag_comb_excess - Gmag_obj_excess)**2)  # squared distance
        lnprob_flux = -0.5 * (np.sum((flx_obs - flx_model) ** 2 / flx_obs_s**2 + np.log(2*np.pi*flx_obs_s**2)))
        # print lnprob_mag, lnprob_flux, lnprob_mag/lnprob_flux
        # lnprob_mag = 0
        # lnprob_flux = 0
        return lnprob_mag + lnprob_flux
    else:
        return -np.inf


def fit_spectra_and_teff(iso, obj_data, flx, flx_s, wvl, G_mag_ref, G_mag_obj,
                         nwalkers=50, n_threds=10, n_steps_1=100, n_steps_2=300, suffix=''):

    gmag_exc = G_mag_ref - G_mag_obj
    print ' Mag exces:', gmag_exc

    # define used subset of the wvl data
    wvl_model = cannon_model.dispersion
    idx_cannon_wvl_mask = get_linelist_mask(wvl_model, d_wvl=0.)

    # resample read spectra to the same wvl pixels as cannon mask
    wvl_new = wvl_model[idx_cannon_wvl_mask]
    flx_new = spectra_resample(flx, wvl, wvl_new, k=1)
    flx_s_new = spectra_resample(flx_s, wvl, wvl_new, k=1)

    teff_range = [np.min(iso['teff']), np.max(iso['teff'])]

    def _p0_generate(teff_init_range, mean_val):
        p0 = list([])
        for i_w in range(nwalkers):
            t_2, t_1 = mean_val - teff_init_range / 2. + np.random.rand(2) * teff_init_range
            # order randomly determined temperatures according the specified rule t1 >= t2 limiting to only one solution
            p0_new = [max(t_1, t_2), min(t_1, t_2)]
            p0.append(p0_new)
        return p0

    # init samples
    p0_1 = _p0_generate(350, obj_data['Teff_cannon'])
    sampler = emcee.EnsembleSampler(nwalkers, 2, lnprob_flx_mag_fit, threads=n_threds,
                                    args=(iso, idx_cannon_wvl_mask, flx_new, flx_s_new, G_mag_ref, G_mag_obj, teff_range, obj_data['Logg_cannon'], obj_data['Fe_H_cannon']))

    print '  Initial MCMC run - {:.0f} steps'.format(n_steps_1)
    p0, lnp, _ = sampler.run_mcmc(p0_1, n_steps_1)
    # plot lnprobs and walkers
    plot_lnprob(sampler.flatchain, get_distribution_peaks(sampler.flatchain, plot=False), str(obj_data['sobject_id']) + '_corner' + suffix + '_1.png')

    # evaluate results from the initial burn, compute new priors accordingly based on the best lnprobs in the last step
    idx_lnp_sort = np.argsort(lnp)[::-1][:nwalkers/2]  # select only half of the best
    teff_best_median = np.nanmedian(p0[idx_lnp_sort, :], axis=0)
    p0_2 = _p0_generate(75, teff_best_median)
    print '  Intermediate teff:', teff_best_median

    print '  Final MCMC run - {:.0f} steps'.format(n_steps_2)
    # reset the sampler chain and lnprobability array, rerun the problem in final push
    sampler.reset()
    p0, lnp, _ = sampler.run_mcmc(p0_2, n_steps_2)
    # plot lnprobs and walkers
    plot_lnprob(sampler.flatchain, get_distribution_peaks(sampler.flatchain, plot=False), str(obj_data['sobject_id']) + '_corner' + suffix + '_2.png')

    if n_threds > 1:
        sampler.pool.close()
    return sampler, p0, lnp


def get_distribution_peaks(chain_vals, plot_ref='', plot=False):
    n_par = chain_vals.shape[1]
    peak_center = list([])
    for i_p in range(n_par):
        data = chain_vals[:, i_p]
        hist, bins = np.histogram(data, range=(np.percentile(data, .1), np.percentile(data, 99.9)), bins=150)
        d_bin = bins[1]-bins[0]
        bins = bins[:-1] + d_bin/2.
        model = SkewedGaussianModel()
        params = model.make_params(amplitude=np.max(hist), center=np.nanmedian(data), sigma=1, gamma=0)
        result = model.fit(hist, params, x=bins)
        # print result.fit_report()
        hist_peak = bins[np.argmax(result.best_fit)]
        peak_center.append(hist_peak)
        if plot:
            # output results
            plt.plot(bins, hist, c='black')
            plt.axvline(x=hist_peak)
            plt.plot(bins, result.best_fit, c='red')
            plt.savefig(plot_ref+'_chainvals_'+str(i_p)+'.png', dpi=250)
            plt.close()
    return peak_center


def plot_spectra_comparison(obj_data, res_table, flx, wvl, path='plot.png', title=''):
    print ' Ploting spectra comparison:', path
    res_use = res_table[res_table['sobject_id'] == obj_data['sobject_id']]

    # prepare original spectra
    wvl_model = cannon_model.dispersion
    idx_cannon_wvl_mask = get_linelist_mask(wvl_model, d_wvl=0.)

    # resample read spectra to the same wvl pixels as cannon mask
    wvl_new = wvl_model[idx_cannon_wvl_mask]
    flx_new = spectra_resample(flx, wvl, wvl_new, k=1)

    # get spectra as determined byy the cannon model
    flx_cannon = get_cannon(obj_data['Teff_cannon'], obj_data['Logg_cannon'], obj_data['Fe_H_cannon'])[idx_cannon_wvl_mask]
    j_ratio = 10 ** (-0.4 * res_use['Gmag_iso2'][0]) / 10 ** (-0.4 * res_use['Gmag_iso1'][0])
    flx_1 = get_cannon(res_use['teff_1'][0], obj_data['Logg_cannon'], obj_data['Fe_H_cannon'])[idx_cannon_wvl_mask]
    flx_2 = get_cannon(res_use['teff_2'][0], obj_data['Logg_cannon'], obj_data['Fe_H_cannon'])[idx_cannon_wvl_mask]
    flx_model = 1. / (1. + j_ratio) * flx_1 + 1. / (1. + (1. / j_ratio)) * flx_2

    # plot them all and their differences
    plt.figure(1, figsize=(14, 7))
    # [left, bottom, width, height]
    axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
    axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
    axSpectra.plot(flx_new, c='black', lw=1., label='Observed')
    axSpectra.plot(flx_cannon, c='C1', lw=.5, label='Synt cannon')
    axSpectra.plot(flx_model, c='C2', lw=.5, label='Synt composite')
    axSpectra.set(ylim=(0.3, 1.05), xlim=(0, len(flx_new)))
    axSpectra.legend()

    axDiff.axhline(y=0, c='black', lw=1)
    axDiff.plot(flx_new - flx_cannon, c='C1', lw=.5, label='Synt cannon')
    axDiff.plot(flx_new - flx_model, c='C2', lw=.5, label='Synt composite')
    axDiff.set(ylim=(-0.05, 0.05), xlim=(0, len(flx_new)))

    if title is not None:
        axSpectra.set_title(title)
    plt.savefig(path, dpi=300)
    plt.close()


def combine_abs_lines_velocity_space(flx, wvl, element, path='plot.png'):
    print ' Combine lines for element:', element
    lines_use = galah_linelist[galah_linelist['Element'] == element]
    n_lines = len(lines_use)
    if n_lines < 10:
        print '  Not enough lines found'
    wvl_vel_use = np.linspace(-25., 25., 180)  # in km/s

    flx_velspace = np.full((n_lines, len(wvl_vel_use)), np.nan)
    for i_l in range(n_lines):
        wvl_vel = (1. - lines_use[i_l]['line_centre']/wvl) * 299792.458
        # check if there are any pixels in the vicinity of 0 km/s
        n_vicinity = np.sum(np.abs(wvl_vel) < 15)
        if n_vicinity < 10:
            # line not found in spectrum
            continue
        flx_new_vel = spectra_resample(flx, wvl_vel, wvl_vel_use, k=1)
        flx_velspace[i_l, :] = flx_new_vel

        # add individual line to the plot
        if i_l == 1:
            p_label = 'All'
        else:
            p_label = ''
        plt.plot(wvl_vel_use, flx_new_vel, label=p_label, lw=0.5, alpha=0.2, c='black')

    # compute mean, median and std of collected lines
    plt.plot(wvl_vel_use, np.nanmedian(flx_velspace, axis=0), label='Median', lw=1)
    plt.plot(wvl_vel_use, np.nanmean(flx_velspace, axis=0), label='Mean', lw=1)
    plt.title(element)
    plt.xlabel('Velocity [km/s]')
    plt.ylabel('Flux')
    plt.legend()
    plt.savefig(path, dpi=250)
    plt.close()