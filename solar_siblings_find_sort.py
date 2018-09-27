from os import chdir, getcwd
from solar_siblings_functions import *
from scipy.stats import kde


def get_snr_metric_functions(data_dir, band, flux):
    csv_file = 'metrices_snr_function_b{:.0f}_flux{:.2f}.csv'.format(band, flux)
    csv_path = data_dir + csv_file
    if os.path.isfile(csv_path):
        return Table.read(csv_path)
    else:
        return Table()


def get_snr_metric_raw(data_dir, band, flux, band_suffix=None):
    if band_suffix is not None:
        csv_file = 'solar_similarity_narrow_b'+band_suffix+'_flux{:.2f}.fits'.format(flux)
    else:
        csv_file = 'solar_similarity_narrow_b{:.0f}_flux{:.2f}.fits'.format(band, flux)
    csv_path = data_dir + csv_file
    if os.path.isfile(csv_path):
        return Table.read(csv_path)
    else:
        return Table()


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
cannon_params = Table.read(galah_data_input+'sobject_iraf_iDR2_180325_cannon.fits')
# cannon_params = cannon_params[cannon_params['flag_cannon'] >= 0]  # no flagging at this point
# cannon_params = cannon_params.filled()

# define directory with simulations of metrics SNR functions
# distance/similarity measurements
# snr_functions_dir = getcwd() + '/' + 'Distances_SNR_models_guesslike-alllines_gauss_oklinesonly' + '/'
# chdir('Distances_Step1_p0_SNRsamples0_ext0_oklinesonly_G20180327_C180325')
combined_bands = False

params_str = '6000_4.18_0.00'
snr_functions_dir = galah_data_input + 'Distances_SNR_models_'+params_str+'_guesslike-alllines_gauss_oklinesonly' + '/'
chdir(galah_data_input + 'Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_G20180327_C180325_refpar_'+params_str)

# snr_functions_dir = os.getcwd() + '/' + 'Distances_SNR_models_guesslike-alllines_gauss_oklinesonly_origsamp_nonoise' + '/'
# chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_comb')
# combined_bands = True

evaluate_bands = list([1,2,3,4])
plot_flux_offsets = [0., 0.05, 0.1, 0.15, 0.2]  # [0., 0.04, 0.08, 0.12, 0.16, 0.2]
snr_multi = 1.
PLOT_RESULTS = False

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
    params_joined = join(sim_res, galah_params, keys='sobject_id')
    params_joined = join(sim_res, cannon_params, keys='sobject_id', join_type='left')

    # # plot coorelations between different metrices
    # plt.scatter(sim_res['canberra'], sim_res['correlation'], s=2, lw=0, c='C2', alpha=0.25, label='correlation')
    # plt.scatter(sim_res['canberra'], sim_res['euclidean'], s=2, lw=0, c='C1', alpha=0.25, label='euclidean')
    # plt.scatter(sim_res['canberra'], sim_res['braycurtis'], s=2, lw=0, c='C3', alpha=0.25, label='braycurtis')
    # plt.scatter(sim_res['canberra'], sim_res['chi2'], s=2, lw=0, c='C4', alpha=0.25, label='chi2')
    # plt.xlim(0., 0.06)
    # plt.ylim(0., 0.12)
    # plt.xlabel('Canberra similarity metric')
    # plt.ylabel('Other similarity metrics')
    # plt.legend(markerscale=6)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig('sim_combine_b' + plot_suffix + '.png', dpi=250)
    # plt.close()
    # continue

    # # colourful plot by parameters
    # c_data_flag = params_joined['cannon_flag']
    # for c_col in ['Teff_cannon', 'Fe_H_cannon', 'Logg_cannon']:
    #     print ' Param plot: '+c_col
    #     c_data = params_joined[c_col]
    #     c_data[c_data_flag > 0] = np.nan  # remove/flag bad parameters
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
        snr_range = np.linspace(10, np.nanpercentile(snr_values, 99.5), 600)

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
            if 'median' not in metric:
                for f_o in plot_flux_offsets:
                    if not combined_bands:
                        snr_metrices_functions = get_snr_metric_functions(snr_functions_dir, i_b, f_o)
                        plt.plot(snr_range, metric_by_snr(snr_metrices_functions, metric, snr_range), lw=1, c='C2', label='{:.2f}'.format(f_o))
                        snr_metrices_raw = get_snr_metric_raw(snr_functions_dir, i_b, f_o)
                        # plt.scatter(snr_metrices_raw['snr_guesslike'], snr_metrices_raw[metric], s=5, color='red', lw=0)
                        n_e = 3
                        plt.errorbar(snr_metrices_raw['snr_guesslike'][::n_e], snr_metrices_raw[metric][::n_e],
                                     yerr=snr_metrices_raw[metric + '_std'][::n_e],
                                     c='C3', alpha=0.75, errorevery=2,
                                     capsize=0, elinewidth=0.75, linewidth=0.5, fmt='o', ms=1.5,
                                     label='Simulations')
                    else:
                        snr_metrices_raw = get_snr_metric_raw(snr_functions_dir, None, f_o, band_suffix=sim_suffix)
                        plt.axhline(snr_metrices_raw[metric][0], label=str(f_o))
                        print snr_metrices_raw[metric][0]

            # plt.show()
            # plt.legend()
            plt.savefig(metric + '_b' + plot_suffix + '_g-all.png', dpi=450)
            plt.close()

        # choose objects to be considered in the next step
        if 'px_' not in metric and 'median' not in metric:
            sel_limit = 0.10
            if not combined_bands:
                max_metric_value = metric_by_snr(get_snr_metric_functions(snr_functions_dir, i_b, sel_limit),
                                                 metric, params_joined[snr_col] * snr_multi)
                metric_2 = 'EW'
                max_metric_value_ew = metric_by_snr(get_snr_metric_functions(snr_functions_dir, i_b, sel_limit),
                                                    metric_2, params_joined[snr_col] * snr_multi)

                idx_selected_sobjects = params_joined[metric] < max_metric_value
                print np.sum(idx_selected_sobjects)
                idx_selected_sobjects = np.logical_and(idx_selected_sobjects, params_joined[metric_2] < max_metric_value_ew)
                print np.sum(idx_selected_sobjects)
                # idx_selected_sobjects = np.logical_and(idx_selected_sobjects,
                #                                        np.abs((params_joined['px_under']-params_joined['px_over']))/((params_joined['px_under']+params_joined['px_over'])) < 0.25)
                # print np.sum(idx_selected_sobjects)

                final_selected_objects = fill_results_dictionary(final_selected_objects, metric,
                                                                 np.int64(list(params_joined['sobject_id'][idx_selected_sobjects].data)))
            else:
                pass



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
    # selected_uniq = uniq_id[repeats_id >= len(evaluate_bands)-1]  # allow 1 miss detection (bad wvl solution)
    seleted_per_metric.append(selected_uniq)
    # print selected_uniq
    txt.write(metric_key + ' ('+str(len(selected_uniq))+'):\n')
    n_max = 150
    n_unq = len(selected_uniq)
    print '  ', n_unq
    if n_unq > 8000:
        txt.write('......................\n')
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
# objects selected by all metrices
uniq_id, repeats_id = np.unique(np.hstack(seleted_per_metric), return_counts=True)
n_eval_metrices = len(final_selected_objects.keys())
selected_uniq = uniq_id[repeats_id >= n_eval_metrices]
selected_uniq_1 = uniq_id[repeats_id >= n_eval_metrices-1]
selected_uniq_2 = uniq_id[repeats_id >= n_eval_metrices-2]
txt.write('All:\n')
txt.write(','.join([str(su) for su in selected_uniq])+'\n\n')
txt.write('All-1:\n')
txt.write(','.join([str(su) for su in selected_uniq_1])+'\n\n')
txt.write('All-2:\n')
txt.write(','.join([str(su) for su in selected_uniq_2])+'\n\n')
txt.close()
