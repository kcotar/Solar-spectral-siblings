from scipy import mgrid
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import varconvolve as varcon

from solar_siblings_functions import *

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize, move_to_dir


def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)


def get_spectra_subset(data, range):
    idx = np.logical_and(data[:, 0] < range[1], data[:, 0] > range[0])
    return data[idx, 1], data[idx, 0]


# EMCEE fiting procedure
def lnprior(params):
    a2 = 0
    a3 = 0
    off, a1, a2, a3 = params
    if -0.5 < off < 0.5 and -1. < a1 < 1. and -1. < a2 < 1. and -1. < a3 < 1.:
        return 0.0
    return -np.inf


def lnprob(params, offset_wvl_range, flx, flx_ref, wvl, idx_compare):
    lp = lnprior(params)
    if not np.isfinite(lp):
        return -np.inf
    flx_new = lnlike(params, offset_wvl_range, flx, wvl)
    flx_diff = flx_new - flx_ref
    return -1. * np.sum(np.abs(flx_diff)[idx_compare])
    # flx_diff_fit = spectra_normalize(wvl - np.mean(wvl), flx_diff, steps=10, sigma_low=2., sigma_high=2.,
    #                                  order=12, func='poly', return_fit=True)
    # return -1.*np.sum(np.abs(flx_diff - flx_diff_fit)[idx_compare])


def lnlike(params, offset_wvl_range, flx, wvl):
    off, a1, a2, a3 = params
    wvl_m = (offset_wvl_range[1] + offset_wvl_range[0]) / 2.  # mean value
    wvl_r = (offset_wvl_range[1] - offset_wvl_range[0])  # range of values
    wvl_d = wvl[1] - wvl[0]
    wvl_new = (wvl - wvl_m) / wvl_r
    flx_off = off + a1*wvl_new * a2*wvl_new**2 + a3*wvl_new**3
    flx_new = spectra_normalize(wvl - np.mean(wvl), flx - flx_off, steps=st, sigma_low=sl,
                                sigma_high=sh, order=ord, func='poly')
    return flx_new


data_dir = '/home/klemen/data4_mount/'
solar_data_dir = data_dir+'Solar_data_dr53/'
galah_linelist = Table.read(data_dir+'GALAH_Cannon_linelist_newer.csv')

# reference spectra
solar_ref = pd.read_csv(solar_data_dir + 'solar_spectra.txt', header=None, delimiter=' ', na_values='nan').values
solar_ref_conv = pd.read_csv(solar_data_dir + 'solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values

# Galah spectrum
galah_ext = 0
twilight_spectrum_file = 'twilight_spectrum_galah_ext'+str(galah_ext)+'_dateall_p9.txt'
solar_galah = pd.read_csv(solar_data_dir + twilight_spectrum_file, header=None, delimiter=' ', na_values='nan').values

# # convolve solar spectrum - different, modified and visually modified values
# conv_min = [0.065, 0.070, 0.075, 0.080]
# conv_max = [0.165, 0.165, 0.165, 0.165]
# plt.plot(solar_galah[:, 0], solar_galah[:, 1]-0.02, label='galah', c='black', lw=2)
# for i_c in range(len(conv_min)):
#     kernel_widths = np.linspace(conv_min[i_c], conv_max[i_c], solar_ref.shape[0])
#     solar_ref_varcon = varcon.varconvolve(solar_ref[:, 0], solar_ref[:, 1], kernel, kernel_widths)
#     plt.plot(solar_ref[:, 0], solar_ref_varcon, label=str(i_c), lw=1)
# plt.legend()
# plt.show()
# plt.close()

# # export convolved spectra
# kernel_widths = np.linspace(0.070, 0.165, solar_ref.shape[0])
# solar_ref[:, 1] = varcon.varconvolve(solar_ref[:, 0], solar_ref[:, 1], kernel, kernel_widths)
# txt = open(solar_data_dir + 'solar_spectra_conv.txt', 'w')
# for i_l in range(solar_ref.shape[0]):
#     txt.write(str(solar_ref[i_l, 0])+' '+str(solar_ref[i_l, 1])+'\n')
# txt.close()

min_wvl = list([4705, 5640, 6475, 7690])
max_wvl = list([4905, 5885, 6750, 7900])

sl = 2.
sh = 3.
st = 11
ord = 1

min_wvl_offset = np.array([4715, 5665, 6485, 7760])
max_wvl_offset = np.array([4890, 5865, 6725, 7840])

# offset analysis
perform_analysis = False
fit_for_best = False
fit_on_all_pixels = True
# final offset correction and output of spectra
final_output = True


# final values of offsets
if galah_ext == 0 and final_output:
    param_fitted = [[0.0413, -0.2059, 0.2227, -0.0448],
                    [0.0235, 0.6426, 0.2755, -0.2195],
                    [0.0281, -0.6960, -0.0896, -0.1212],
                    [-0.0377, -0.4947, 0.4748, 0.2465]]
if galah_ext == 4 and final_output:
    # not defined yet for the newest equation
    pass


def get_spectrum_with_offset(flx, wvl, amp, amp2, off, offset_wvl_range):
    d_wvl = offset_wvl_range[1] - offset_wvl_range[0]
    wvl_new = (wvl-(offset_wvl_range[1]-offset_wvl_range[0])/2.)/d_wvl
    y_off_perwvl = wvl_new**2 * amp2/1e4 + wvl_new * amp + off
    flx_new = spectra_normalize(wvl - np.mean(wvl), flx - y_off_perwvl, steps=st, sigma_low=sl,
                                sigma_high=sh, order=ord, func='poly')
    return flx_new


move_to_dir('Twilight_offset_determine_ext'+str(galah_ext)+'_allpx')

if fit_for_best:
    for i_b in range(4):
        print 'Fit for band:', i_b
        wvl_range = (min_wvl[i_b], max_wvl[i_b])
        flux_offset_wvl_range = (min_wvl_offset[i_b], max_wvl_offset[i_b])
        flx_galah, wvl_galah = get_spectra_subset(solar_galah, wvl_range)
        flx_ref, wvl_ref = get_spectra_subset(solar_ref_conv, wvl_range)

        flx_galah = spectra_normalize(wvl_galah-np.mean(wvl_galah), flx_galah, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
        flx_ref = spectra_normalize(wvl_ref-np.mean(wvl_ref), flx_ref, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')

        idx_sim = np.logical_and(wvl_galah >= min_wvl_offset[i_b], wvl_galah <= max_wvl_offset[i_b])
        idx_sim = np.logical_and(idx_sim, flx_galah < 0.93)
        if not fit_on_all_pixels:
            idx_sim = np.logical_and(idx_sim, get_linelist_mask(wvl_galah, d_wvl=1.))

        ndim, nwalkers = 4, 70
        param_max = np.array([0.25, 0.75, 0.75, 0.75][:ndim])
        labels_use = ['off', 'a1', 'a2', 'a3'][:ndim]
        pos = [np.random.rand(ndim)*2.*param_max - param_max for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=None,
                                        args=(flux_offset_wvl_range, flx_galah, flx_ref, wvl_galah, idx_sim))
        lnprob_lim = -150
        print 'Running emcee - 1'
        fit_res, fit_prob, _ = sampler.run_mcmc(pos, 160)

        # filter out bad walkers and reset them
        idx_bad_w = fit_prob < lnprob_lim
        fit_mean_1 = list(np.median(np.array(fit_res)[np.where(~idx_bad_w)[0], :], axis=0))
        if np.sum(idx_bad_w) > 0:
            for i_b_w in np.where(idx_bad_w)[0]:
                fit_res[i_b_w] = fit_mean_1

        # print 'Running emcee - 2'
        # sampler.reset()
        # fit_res, fit_prob, _ = sampler.run_mcmc(fit_res, 1000)

        print 'Finished emcee'

        print 'Plotting walker probabilities'
        sampler_chain_vals = sampler.flatchain
        kernel_fit = np.median(sampler_chain_vals, axis=0)
        walkers_prob = sampler.lnprobability
        for i_w in range(walkers_prob.shape[0]):
            plt.plot(walkers_prob[i_w, :])
        walkers_prob = walkers_prob.flatten()  # without this correction numpy
        walkers_prob = walkers_prob[np.isfinite(walkers_prob)]  # percentile may return incorrect -inf value
        plt.ylim((np.percentile(walkers_prob, 1), np.percentile(walkers_prob, 99)))
        plt.savefig('b'+str(i_b) + '_lnprob.png', dpi=400)
        plt.close()
        c_fig = corner.corner(sampler.flatchain, truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                              labels=labels_use, bins=30)
        c_fig.savefig('b'+str(i_b) + '_corner.png', dpi=200)
        plt.close(c_fig)
        plt.close()

        kernel_fit_max = fit_res[np.nanargmax(fit_prob)]  # max probable walker in the last step
        print '  Run - m_ok: ', fit_mean_1
        print '  Run - m_all:', kernel_fit
        print '  Run - p_max: ', kernel_fit_max

        # plot few of them
        for fit_res_sel in fit_res[np.argsort(fit_prob)[-3:]]:
            plot_p = str('b{:.0f}_o_{:0.3f}_a1_{:0.3f}_a2_{:0.3f}_a3_{:0.3f}_mcmc.png'.format(i_b, fit_res_sel[0], fit_res_sel[1], fit_res_sel[2], fit_res_sel[3]))
            flx_galah_new = lnlike(fit_res_sel, flux_offset_wvl_range, flx_galah, wvl_galah)
            # plt.plot(flx_galah_new - flx_galah)
            # plt.savefig(plot_p+'_diff.png', dpi=200)
            # plt.close()
            flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_ref - flx_galah_new,
                                             steps=st, sigma_low=2., sigma_high=2., order=12, func='poly',
                                             return_fit=True)
            plot_spectra_with_difference(flx_ref, flx_galah_new, wvl_ref, x_range=wvl_range, linelist=galah_linelist,
                                         title='', path=plot_p, diff_func=flx_diff_fit, flux3=flx_galah)

        # params = Parameters()
        # params.add('off', value=1e-1, min=-5.0, max=5.0, vary=True, brute_step=0.001)
        # params.add('amp', value=1e-1, min=-1.0, max=1.0, vary=True, brute_step=0.001)
        # params.add('amp2', value=1e-1, min=-1.0, max=1.0, vary=True, brute_step=0.001)
        #
        #
        # def tel_scale_func(params, eval=False):
        #     flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, params['amp'], params['amp2'], params['off'], flux_offset_wvl_range)
        #     flx_diff = flx_ref - flx_galah_new
        #     flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_diff,
        #                                      steps=st, sigma_low=2., sigma_high=2., order=12, func='poly',
        #                                      return_fit=True)
        #     if eval:
        #         return 1.
        #     else:
        #         print np.sum(np.abs(flx_diff - flx_diff_fit)[idx_sim])
        #         return np.abs(flx_diff - flx_diff_fit)[idx_sim]
        #         # return (flx_diff - flx_diff_fit)[idx_sim]**2
        #         # return (np.abs(flx_diff)[idx_sim])
        #
        # minner = Minimizer(tel_scale_func, params)
        # result = minner.minimize()#method='brute')
        # res_param = result.params
        # # report_fit(result)
        #
        # results_best = [[res_param['amp'].value, res_param['amp2'].value, res_param['off'].value, flux_offset_wvl_range[0], flux_offset_wvl_range[1]]]
        # print results_best
        #
        # # plot them
        # for res in results_best:
        #     plot_p = str('b{:.0f}_a_{:0.3f}_a2_{:0.3f}_o_{:0.3f}_lmfit_allabs.png'.format(i_b, res[0], res[1], res[2]))
        #     flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, res[0], res[1], res[2], (res[3], res[4]))
        #     flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_ref - flx_galah_new,
        #                                      steps=st, sigma_low=2., sigma_high=2., order=12, func='poly', return_fit=True)
        #     plot_spectra_with_difference(flx_ref, flx_galah_new, wvl_ref, x_range=wvl_range, linelist=galah_linelist,
        #                                  title='', path=plot_p, diff_func=flx_diff_fit)

if perform_analysis:
    # TODO: not yet updated for new offset removal functions
    for i_b in range(4):
        wvl_range = (min_wvl[i_b], max_wvl[i_b])
        flux_offset_wvl_range = (min_wvl_offset[i_b], max_wvl_offset[i_b])
        flx_galah, wvl_galah = get_spectra_subset(solar_galah, wvl_range)
        flx_ref, wvl_ref = get_spectra_subset(solar_ref_conv, wvl_range)

        flx_galah = spectra_normalize(wvl_galah-np.mean(wvl_galah), flx_galah, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
        flx_ref = spectra_normalize(wvl_ref-np.mean(wvl_ref), flx_ref, steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')

        results = []
        idx_sim = np.logical_and(wvl_galah >= min_wvl_offset[i_b], wvl_galah <= max_wvl_offset[i_b])
        idx_sim = np.logical_and(idx_sim, get_linelist_mask(wvl_galah, d_wvl=0.1))
        for flux_offset_amp in np.arange(-0.15, 0.15, 0.01):
            for flux_offset in np.arange(-0.15, 0.15, 0.01):
                print flux_offset, flux_offset_amp
                flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, flux_offset_amp, flux_offset_amp2, flux_offset, flux_offset_wvl_range)
                flx_diff = flx_ref-flx_galah_new
                flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_diff,
                                                 steps=st, sigma_low=2., sigma_high=2., order=12, func='poly', return_fit=True)
                plot_sim = np.nansum(np.abs(flx_diff - flx_diff_fit)[idx_sim])
                results.append([flux_offset, flux_offset_amp, min_wvl_offset[i_b], max_wvl_offset[i_b], plot_sim])
        results = np.vstack(results)
        print 'Best for band', i_b
        idx_best = np.argsort(results[:, 4])
        results_best = results[idx_best[:5], :]
        print results_best

        # plot them
        for res in results_best:
            plot_p = str('b{:.0f}_o{:0.3f}_a{:0.3f}.png'.format(i_b, res[0], res[1]))
            flx_galah_new = get_spectrum_with_offset(flx_galah, wvl_galah, res[0], res[1], res[3], (res[4], res[5]))
            flx_diff_fit = spectra_normalize(wvl_galah - np.mean(wvl_galah), flx_ref - flx_galah_new,
                                             steps=st, sigma_low=2., sigma_high=2., order=12, func='poly', return_fit=True)
            plot_spectra_with_difference(flx_ref, flx_galah_new, wvl_ref, x_range=wvl_range, linelist=galah_linelist,
                                         title=str(res[4]), path=plot_p, diff_func=flx_diff_fit)

if final_output:
    # output new spectra as a csv format for later use

    flx_galah_offset = np.full(solar_galah.shape[0], np.nan)
    wvl_galah_offset = solar_galah[:, 0]

    for i_b in range(4):
        wvl_range = (min_wvl[i_b], max_wvl[i_b])
        flx_galah, wvl_galah = get_spectra_subset(solar_galah, wvl_range)
        flux_offset_wvl_range = (min_wvl_offset[i_b], max_wvl_offset[i_b])

        flx_galah_new = lnlike(param_fitted[i_b], flux_offset_wvl_range, flx_galah, wvl_galah)
        idx_use = np.in1d(wvl_galah_offset, wvl_galah)
        flx_galah_offset[idx_use] = flx_galah_new

    # output to file
    out_file = solar_data_dir + twilight_spectrum_file[:-4]+'_offset.txt'
    txt = open(out_file, 'w')
    for i_l in range(len(wvl_galah_offset)):
        txt.write(str(wvl_galah_offset[i_l]) + ' ' + str(flx_galah_offset[i_l]) + '\n')
    txt.close()
