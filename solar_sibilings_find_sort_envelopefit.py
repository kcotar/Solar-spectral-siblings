import sys
from os import chdir, getcwd
from solar_siblings_functions import *
from scipy.stats import kde


def metrices_to_investigate(all_cols):
    metrices = all_cols[1:]  # remove first sobject_id col
    remove_metrices = np.array(['C', 'Ce', 'Co', 'Eu', 'K', 'La', 'Li', 'Mo', 'Nd', 'Ru', 'Sm', 'Sr', 'Zr',
                                'correlation', 'snr_spectrum', 'sqeuclidean', 'minkowski', 'sum'])
    metrices = [s for s in metrices if np.sum(remove_metrices == s) <= 0]
    metrices = [s for s in metrices if '_std' not in s]
    metrices = [s for s in metrices if '_min' not in s]
    metrices = [s for s in metrices if '_max' not in s]
    # metrices = [s for s in metrices if 'px_' not in s]
    return metrices


def fill_results_dictionary(res_dict, key, values):
    if key not in res_dict:
        res_dict[key] = list([])
    res_dict[key].append(values)
    return res_dict


# read reference solar data
suffix = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'
solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix)

# read Galah guess and/or cannon parameters
galah_params = Table.read(galah_data_input+'sobject_iraf_53_reduced_20180327.fits')
# cannon_params = Table.read(galah_data_input+'sobject_iraf_iDR2_180325_cannon.fits')
# cannon_params = cannon_params[cannon_params['flag_cannon'] >= 0]  # no flagging at this point
# cannon_params = cannon_params.filled()

results_dir = '/shared/data-camelot/cotar/_Multiples_binaries_results_iDR3/'
# define directory with simulations of metrics SNR functions
# distance/similarity measurements

argv = sys.argv
if len(argv) >= 2:
    perc_best = int(argv[1])
else:
    perc_best = 5

if len(argv) >= 3:
    idx_run_params = int(argv[2])
    params_str = ['5100_4.55_0.00', '5200_4.53_0.00', '5300_4.51_0.00', '5400_4.48_0.00',
                  '5500_4.46_0.00', '5600_4.43_0.00', '5700_4.40_0.00', '5800_4.37_0.00',
                  '5900_4.34_0.00', '6000_4.30_0.00'][idx_run_params]
    chdir(results_dir + 'Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_G20180327_C181221_withH_refpar_' + params_str)
    print params_str
else:
    chdir(results_dir + 'Distances_Step1_p0_SNRsamples0_ext0_oklinesonly_G20180327_C181221_withHwings')

combined_bands = False

evaluate_bands = list([1,2,3,4])
plot_flux_offsets = [0., 0.05, 0.1, 0.15, 0.2]  # [0., 0.04, 0.08, 0.12, 0.16, 0.2]
snr_multi = 1.
PLOT_RESULTS = True

y_lim_plot = {'braycurtis':0.025, 'canberra':0.025, 'chebyshev':0.25, 'chi2':0.03, 'cityblock':0.045,
              'cosine':0.02, 'minkowski':0.00015, 'sqeuclidean':0.0005, 'wchebyshev':0.3, 'euclidean':0.04}
v_lim_params = {'Teff_cannon':(5400, 5850), 'Fe_H_cannon':(-0.3, 0.15), 'Logg_cannon':(3.9, 4.5)}
final_selected_objects = {}

for i_b in evaluate_bands:
    if combined_bands:
        sim_suffix = ''.join([str(sb) for sb in evaluate_bands])
        print 'Bands', sim_suffix
        sim_res = Table.read('solar_similarity_b'+sim_suffix+'.fits')
        plot_suffix = sim_suffix
    else:
        print 'Band', i_b
        sim_res = Table.read('solar_similarity_b{:.0f}.fits'.format(i_b))
        plot_suffix = str(i_b)
    params_joined = join(sim_res, galah_params, keys='sobject_id', join_type='left')

    d_snr = 2
    w_snr = 4
    snr_col = 'snr_spectrum'
    snr_multi = 1.
    plot_suffix += '_{:02.0f}'.format(perc_best)

    for metric in ['canberra']: #metrices_to_investigate(sim_res.colnames):
        print ' Metric:', metric
        metric_values = params_joined[metric]
        snr_values = params_joined[snr_col] * snr_multi

        min_snr, max_snr = np.nanpercentile(snr_values, [0.5, 99.])
        print ' Min/max snr:', min_snr, max_snr
        eval_at_snrs = np.arange(min_snr, max_snr, d_snr)

        eval_at_snrs_maxdis = list([])
        for eval_snr in eval_at_snrs:
            idx_bin = np.logical_and(snr_values >= eval_snr - w_snr/2.,
                                     snr_values < eval_snr + w_snr/2.)
            eval_at_snrs_maxdis.append(np.percentile(metric_values[idx_bin], perc_best))
        eval_at_snrs_maxdis = np.array(eval_at_snrs_maxdis)

        # fit a function to those points
        pow_multi = 1.
        sim_med = np.median(eval_at_snrs_maxdis)
        if sim_med > 100:
            pow_amp = 250.
        else:
            pow_amp = 1.
        f_init = models.Shift(offset=0.) | (models.Const1D(amplitude=pow_multi) * models.PowerLaw1D(amplitude=pow_amp, x_0=1.0, alpha=1.0) + models.Linear1D(slope=0, intercept=sim_med))
        fitter = fitting.LevMarLSQFitter()
        gg_fit = fitter(f_init, eval_at_snrs, eval_at_snrs_maxdis)

        max_metric_value = gg_fit(snr_values)
        idx_selected_sobjects = metric_values < max_metric_value
        print ' Selected:', np.sum(idx_selected_sobjects)
        final_selected_objects = fill_results_dictionary(final_selected_objects, metric,
                                                         np.int64(list(params_joined['sobject_id'][idx_selected_sobjects].data)))

        plt.scatter(snr_values, metric_values, s=2, lw=0, alpha=0.25, color='black')
        plt.scatter(eval_at_snrs, eval_at_snrs_maxdis, s=10, lw=0, alpha=1, color='C3')
        plt.plot(eval_at_snrs, eval_at_snrs_maxdis, lw=2, alpha=1, color='C3')
        plt.plot(eval_at_snrs, gg_fit(eval_at_snrs), lw=2, alpha=1, color='C0')
        x_lim = (10, np.nanpercentile(snr_values, 99.))
        plt.ylim(0, np.nanpercentile(metric_values, 90))
        plt.title(metric + ' band:' + plot_suffix)
        plt.xlim(x_lim)
        plt.xlabel(snr_col)
        plt.ylabel(metric)
        # plt.show()
        plt.savefig(metric + '_b' + plot_suffix + '_envelope.png', dpi=450)
        plt.close()

txt_out_selection = 'final_selection_{:02.0f}_envelope.txt'.format(perc_best)
txt = open(txt_out_selection, 'w')
seleted_per_metric = list([])
for metric_key in final_selected_objects.keys():
    print 'Selecting objects for metric:', metric_key
    # print final_selected_objects[metric_key]
    selected = np.hstack(final_selected_objects[metric_key])
    uniq_id, repeats_id = np.unique(selected, return_counts=True)
    selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)]
    # selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)-1]  # allow 1 miss detection (bad wvl solution)
    seleted_per_metric.append(selected_uniq)
    # print selected_uniq
    # txt.write(metric_key + ' ('+str(len(selected_uniq))+'):\n')
    n_max = 100
    n_unq = len(selected_uniq)
    print '  ', n_unq
    if n_unq > 5000:
        # txt.write('......................\n')
        pass
    else:
        if n_unq > n_max:
            # write out in multiple lines
            for i_l in range(int(np.ceil(1.*n_unq/n_max))):
                if i_l > 0:
                    txt.write(',')
                txt.write(','.join([str(su) for su in selected_uniq[n_max*i_l:n_max*(i_l+1)]])+'\n')
        else:
            txt.write(','.join([str(su) for su in selected_uniq])+'\n')
    txt.write('\n')
txt.close()