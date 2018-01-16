from os import chdir
from solar_siblings_functions import *


def get_snr_metric_functions(data_dir, band, flux):
    csv_file = 'metrices_snr_function_b{:.0f}_flux{:.2f}.csv'.format(band, flux)
    csv_path = data_dir + csv_file
    if os.path.isfile(csv_path):
        return Table.read(csv_path)
    else:
        return Table()


def metrices_to_investigate(all_cols):
    metrices = all_cols[1:]  # remove first sobject_id col
    remove_metrices = np.array(['C', 'Ce', 'Co', 'Eu', 'K', 'La', 'Li', 'Mo', 'Nd', 'Ru', 'Sm', 'Sr', 'Zr', 'EW', 'correlation'])
    metrices = [s for s in metrices if np.sum(remove_metrices == s) <= 0]
    metrices = [s for s in metrices if '_std' not in s]
    return metrices


def fill_results_dictionary(res_dict, key, values):
    if key not in res_dict:
        res_dict[key] = list([])
    res_dict[key].append(values)
    return res_dict


# read reference solar data
suffix = '_ext0_2_offset'
solar_input_dir = galah_data_input+'Solar_data/'
solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix)

# read Galah guess and/or cannon parameters
galah_params = Table.read(galah_data_input+'sobject_iraf_52_reduced_20171111.fits')

# define directory with simulations of metrics SNR functions
snr_functions_dir = os.getcwd() + '/' + 'Distances_SNR-functions_multioffset_allbands_subsample' + '/'

# distance/similarity measurements
chdir('Distances_Step1_p0_SNRsamples0')
evaluate_bands = list([1, 2, 4])
plot_flux_offsets = [0., 0.04, 0.08, 0.12, 0.16, 0.2]

final_selected_objects = {}

for i_b in evaluate_bands:
    print 'Band', i_b
    sim_res = Table.read('solar_similarity_b{:.0f}.fits'.format(i_b))
    params_joined = join(sim_res, galah_params, keys='sobject_id')

    for metric in metrices_to_investigate(sim_res.colnames):
        print ' Metric:', metric
        metric_values = params_joined[metric]
        snr_col = 'snr_c' + str(i_b) + '_guess'
        snr_values = params_joined[snr_col]
        snr_range = np.linspace(np.min(snr_values), np.max(snr_values), 200)
        plt.errorbar(snr_values, metric_values, yerr=params_joined[metric + '_std'], fmt='o', ms=1, elinewidth=0.3,
                     alpha=0.3)
        plt.ylim(0, np.nanpercentile(metric_values, 92))
        plt.title(metric + ' band:' + str(i_b))
        plt.xlim(0, np.nanpercentile(snr_values, 99.8))
        plt.xlabel(snr_col)
        plt.ylabel(metric)
        # add SNR function for observed metric
        for f_o in plot_flux_offsets:
            snr_metrices_functions = get_snr_metric_functions(snr_functions_dir, i_b, f_o)
            plt.plot(snr_range, metric_by_snr(snr_metrices_functions, metric, snr_range), lw=1, c='red', label='{:.2f}'.format(f_o))
        # plt.show()
        plt.savefig(metric + '_b' + str(i_b) + '.png', dpi=450)
        plt.close()

        # choose objects to be considered in the next step
        max_metric_value = metric_by_snr(get_snr_metric_functions(snr_functions_dir, i_b, 0.1),
                                         metric, params_joined[snr_col])
        idx_selected_sobjects = params_joined[metric] < max_metric_value
        final_selected_objects = fill_results_dictionary(final_selected_objects, metric,
                                                         np.int64(list(params_joined['sobject_id'][idx_selected_sobjects].data)))

txt_out_selection = 'final_selection.txt'
txt = open(txt_out_selection, 'w')
seleted_per_metric = list([])
for metric_key in final_selected_objects.keys():
    print 'Selecting objects for metric:', metric_key
    # print final_selected_objects[metric_key]
    selected = np.hstack(final_selected_objects[metric_key])
    uniq_id, repeats_id = np.unique(selected, return_counts=True)
    selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)]
    seleted_per_metric.append(selected_uniq)
    print selected_uniq
    txt.write(metric_key + ':\n')
    txt.write(','.join([str(su) for su in selected_uniq]))
    txt.write('\n\n')
# objects selected by all metrices
uniq_id, repeats_id = np.unique(np.hstack(seleted_per_metric), return_counts=True)
selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)]
txt.write('All:\n')
txt.write(','.join([str(su) for su in selected_uniq]))
txt.close()

'''
params_joined = join(galah_params, sim_res, keys='sobject_id')
for i_b in [2]:
    for metric in sim_metrics:
        metric_values = params_joined[metric]
        if np.sum(np.isfinite(metric_values)) <= 0:
            continue
        snr_col = 'snr_c'+str(i_b)+'_guess'
        snr_values = params_joined[snr_col]
        snr_range = np.linspace(np.min(snr_values), np.max(snr_values), 200)
        # plt.scatter(snr_values, metric_values, lw=0, s=2)
        plt.errorbar(snr_values, metric_values, yerr=params_joined[metric+'_std'], fmt='o', ms=1, elinewidth=0.3, alpha=0.3)
        # y_perc = np.percentile(metric_values, 500./len(metric_values)*100)
        # plt.axhline(y=y_perc, c='black')
        plt.ylim(0, np.nanpercentile(metric_values, 92))
        plt.title(metric+' band:'+str(i_b))
        plt.xlim(0, np.nanpercentile(snr_values, 99.8))
        plt.xlabel(snr_col)
        plt.ylabel(metric)
        # add SNR function for observed metric
        plt.plot(snr_range, metric_by_snr(snr_metrices_functions1, metric, snr_range), lw=1, c='red', label='0.00')
        plt.plot(snr_range, metric_by_snr(snr_metrices_functions6, metric, snr_range), lw=1, c='red', label='0.10')
        # plt.show()
        plt.savefig(metric+'_b'+str(i_b)+'.png', dpi=450)
        plt.close()

        if '_std' not in metric:
            # compute new values of similarities
            new_metric_value = metric_values - metric_by_snr(snr_metrices_functions1, metric, snr_values)
            plt.scatter(snr_values, new_metric_value, lw=0, s=2, alpha=0.3, label='0.00')
            plt.xlabel(snr_col)
            plt.ylabel(metric)
            plt.ylim(0, np.nanpercentile(new_metric_value, 90))
            plt.xlim(0, np.nanpercentile(snr_values, 99.8))
            plt.savefig(metric + '_b' + str(i_b) + '_SNRnorm.png', dpi=300)
            plt.close()

sim_metrics = [s for s in sim_metrics if '_std' not in s]
# remove metrices
sim_metrics = [s for s in sim_metrics if np.sum(remove_metrices == s) <= 0]
print sim_metrics

for metric in sim_metrics:

    print '-----------------------------'
    print metric
    if not np.isfinite(sim_res[metric]).all():
        continue

    metric_value_snr = metric_by_snr(snr_metrices_functions4, metric, params_joined[snr_col])
    idx_lower = params_joined[metric] < metric_value_snr
    s_lower = params_joined['sobject_id'][idx_lower]
    print ','.join([str(a) for a in s_lower])

    idx_sort = np.argsort(sim_res[metric]-metric_value_snr)
    sim_res[metric][idx_sort] = scores

    print '-- Best'
    s_best = sim_res[idx_sort]['sobject_id'][:120]
    print ','.join([str(a) for a in s_best])
    # plot_spectra(s_best, dr52_dir, galah_params, solar_flx, solar_wvl,
    #              galah_linelist=galah_linelist, save_path=metric + '_b2_spectra_best.png',
    #              band=2, ext=0, y_range=(0.7, 1.05), x_range=(5680, 5710))
    print '-- Worst'
    s_worst = sim_res[idx_sort]['sobject_id'][-120:]
    print ','.join([str(a) for a in s_worst])
    # plot_spectra(s_worst, dr52_dir, galah_params, solar_flx, solar_wvl,
    #              galah_linelist=galah_linelist, save_path=metric + '_b2_spectra_worst.png',
    #              band=2, ext=0, y_range=(0.7, 1.05), x_range=(5680, 5710))

# final
sim_res['final_score'] = np.nanmean(sim_res[sim_metrics].to_pandas().values, axis=1)
print '-- Final best'
idx_sort = np.argsort(sim_res['final_score'])
s_best = sim_res[idx_sort]['sobject_id'][:120]
print ','.join([str(a) for a in s_best])
print '-- Final worst'
s_worst = sim_res[idx_sort]['sobject_id'][-120:]
print ','.join([str(a) for a in s_worst])
# plot_spectra(s_best, dr52_dir, galah_params, solar_flx, solar_wvl,
#                  galah_linelist=galah_linelist, save_path='all_b2_spectra_best.png',
#                  band=2, ext=0, y_range=(0.7, 1.05), x_range=(5680, 5710))

# plt.scatter(params_joined['chebyshev'], params_joined['minkowski'])
# plt.show()
# plt.close()

print sim_res[list(np.hstack(('sobject_id',sim_metrics,'final_score')))]
'''
