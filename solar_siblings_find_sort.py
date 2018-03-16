from os import chdir
from solar_siblings_functions import *
from scipy.stats import kde


def get_snr_metric_functions(data_dir, band, flux):
    csv_file = 'metrices_snr_function_b{:.0f}_flux{:.2f}.csv'.format(band, flux)
    csv_path = data_dir + csv_file
    if os.path.isfile(csv_path):
        return Table.read(csv_path)
    else:
        return Table()


def metrices_to_investigate(all_cols):
    metrices = all_cols[1:]  # remove first sobject_id col
    remove_metrices = np.array(['C', 'Ce', 'Co', 'Eu', 'K', 'La', 'Li', 'Mo', 'Nd', 'Ru', 'Sm', 'Sr', 'Zr', 'EW', 'correlation', 'snr_spectrum'])
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
galah_params = Table.read(galah_data_input+'sobject_iraf_53_reduced_20180222.fits')
cannon_params = Table.read(galah_data_input+'sobject_iraf_iDR2_180108_cannon.fits')
cannon_params = cannon_params[cannon_params['flag_cannon'] == 0]

# define directory with simulations of metrics SNR functions
snr_functions_dir = os.getcwd() + '/' + 'Distances_SNR_models_subsample_guesslike-alllines_gauss_oklinesonly' + '/'

# distance/similarity measurements

chdir('Distances_Step1_p0_SNRsamples0_ext0_oklinesonly')
combined_bands = False

# chdir('Distances_Step2_p0_SNRsamples500_ext0_oklinesonly_origsamp_no-offset-fit')
# combined_bands = True

evaluate_bands = list([1, 2, 3, 4])
plot_flux_offsets = [0., 0.05, 0.1, 0.15, 0.2]  # [0., 0.04, 0.08, 0.12, 0.16, 0.2]
snr_multi = 1.
PLOT_RESULTS = True

y_lim_plot = {'braycurtis':0.025, 'canberra':0.025, 'chebyshev':0.25, 'chi2':0.03, 'cityblock':0.045,
              'cosine':0.02, 'minkowski':0.00015, 'sqeuclidean':0.0005, 'wchebyshev':0.3, 'euclidean':0.04}
v_lim_params = {'Teff_cannon':(5400, 5850), 'Feh_cannon':(-0.3, 0.15), 'Logg_cannon':(3.9, 4.5)}
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
    params_joined = join(sim_res, galah_params, keys='sobject_id')
    params_joined = join(sim_res, cannon_params, keys='sobject_id', join_type='left')

    # # colourful plot by parameters
    # for c_col in ['Teff_cannon', 'Feh_cannon', 'Logg_cannon']:
    #     print ' Param plot: '+c_col
    #     c_data = params_joined[c_col]
    #     plt.scatter(params_joined['snr_spectrum'], params_joined['chi2'],
    #                 lw=0, s=0.75, c=c_data, cmap='viridis',
    #                 vmin=v_lim_params[c_col][0], vmax=v_lim_params[c_col][1])
    #     plt.ylim(0, np.nanpercentile(params_joined['chi2'], 90))
    #     plt.xlim(10, np.nanpercentile(params_joined['snr_spectrum'], 99.))
    #     plt.colorbar()
    #     plt.savefig('cannon_chi2_b' + plot_suffix + '_'+c_col+'.png', dpi=450)
    #     plt.close()

    for metric in metrices_to_investigate(sim_res.colnames):
        print ' Metric:', metric
        metric_values = params_joined[metric]
        snr_col = 'snr_spectrum'
        snr_values = params_joined[snr_col] * snr_multi
        snr_range = np.linspace(10, np.nanpercentile(snr_values, 98.), 600)

        if PLOT_RESULTS:
            if 'median' in metric:
                y_lim = (np.nanpercentile(metric_values, 0.5), np.nanpercentile(metric_values, 99.5))
                plt.scatter(snr_values, metric_values, s=2, lw=0, alpha=0.25, color='black')
            else:
                # check if has min max error bat
                y_lim = (0, np.nanpercentile(metric_values, 90))
                plt.errorbar(snr_values, metric_values, yerr=params_joined[metric + '_std'], fmt='.', ms=2,
                             elinewidth=0.3, alpha=0.25, color='black', markeredgewidth=0)
                # plt.errorbar(snr_values, metric_values, yerr=[params_joined[metric + '_min'], params_joined[metric + '_max']], fmt='.', ms=2,
                #              elinewidth=0.3, alpha=0.25, color='black', markeredgewidth=0)

            # # show flagged spectra in this scatter plot
            # if combined_bands:
            #     idx_red_flag = params_joined['red_flag'] > 0
            # else:
            #     idx_red_flag = np.bitwise_and(params_joined['red_flag'], 2**(i_b-1)) > 0
            #     if i_b >= 3:
            #         # also search for molecfit problems
            #         idx_red_flag = np.logical_or(idx_red_flag,
            #                                      np.bitwise_and(params_joined['red_flag'], 2**(i_b - 3 + 4)) > 0)
            # # add them to the plot
            # if np.sum(idx_red_flag) > 0:
            #     plt.scatter(snr_values[idx_red_flag], metric_values[idx_red_flag], s=2, lw=0, alpha=1, color='C1')

            x_lim = (10, np.nanpercentile(snr_values, 99.))
            plt.ylim(y_lim)
            plt.title(metric + ' band:' + plot_suffix)
            plt.xlim(x_lim)
            plt.xlabel(snr_col)
            plt.ylabel(metric)

            # add SNR function for observed metric
            if 'median' not in metric and not combined_bands:
                for f_o in plot_flux_offsets:
                    snr_metrices_functions = get_snr_metric_functions(snr_functions_dir, i_b, f_o)
                    plt.plot(snr_range, metric_by_snr(snr_metrices_functions, metric, snr_range), lw=1, c='C2', label='{:.2f}'.format(f_o))
            # plt.show()
            plt.savefig(metric + '_b' + plot_suffix + '_g-all.png', dpi=450)
            plt.close()

        # choose objects to be considered in the next step
        if 'px_' not in metric:
            sel_limit = 0.1
            max_metric_value = metric_by_snr(get_snr_metric_functions(snr_functions_dir, i_b, sel_limit),
                                             metric, params_joined[snr_col] * snr_multi)
            idx_selected_sobjects = params_joined[metric] < max_metric_value
            final_selected_objects = fill_results_dictionary(final_selected_objects, metric,
                                                             np.int64(list(params_joined['sobject_id'][idx_selected_sobjects].data)))

    if combined_bands:
        # can be executed only once
        break

txt_out_selection = 'final_selection_{:.2f}.txt'.format(sel_limit)
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
    txt.write(metric_key + ' ('+str(len(selected_uniq))+'):\n')
    if len(selected_uniq) > 1000:
        txt.write('...........')
    else:
        txt.write(','.join([str(su) for su in selected_uniq]))
    txt.write('\n\n')
# objects selected by all metrices
uniq_id, repeats_id = np.unique(np.hstack(seleted_per_metric), return_counts=True)
selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)]
txt.write('All:\n')
txt.write(','.join([str(su) for su in selected_uniq]))
txt.close()
