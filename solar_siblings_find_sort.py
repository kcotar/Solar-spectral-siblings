from os import chdir
from solar_siblings_functions import *

chdir('Distances_SNRadded_weighted-p0_test-SNR-norm-ndata-b2-snrnone')
sim_res = Table.read('solar_similarity_narrow.fits')
snr_metrices_functions1 = Table.read('metrices_snr_function_flux0.00.csv')
snr_metrices_functions2 = Table.read('metrices_snr_function_flux0.02.csv')
snr_metrices_functions3 = Table.read('metrices_snr_function_flux0.04.csv')
snr_metrices_functions4 = Table.read('metrices_snr_function_flux0.06.csv')
snr_metrices_functions5 = Table.read('metrices_snr_function_flux0.08.csv')
snr_metrices_functions6 = Table.read('metrices_snr_function_flux0.10.csv')

galah_data_input = '/home/klemen/data4_mount/'
galah_params = Table.read(galah_data_input+'sobject_iraf_52_reduced_20171111.fits')
galah_linelist = Table.read(galah_data_input + 'GALAH_Cannon_linelist_newer.csv')

suffix = '_ext0_2_offset'
solar_input_dir = galah_data_input+'Solar_data/'
solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix)

remove_metrices = np.array(['C', 'Ce', 'Co', 'Eu', 'K', 'La', 'Li', 'Mo', 'Nd', 'Ru', 'Sm', 'Sr', 'Zr', 'EW', 'correlation'])
sim_metrics = sim_res.colnames[1:]
sim_metrics = [s for s in sim_metrics if '_std' not in s]
scores = np.arange(len(sim_res))+1

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

