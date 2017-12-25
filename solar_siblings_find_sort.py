from os import chdir
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np

chdir('Distances_SNRadded_weighted-p1')
sim_res = Table.read('solar_similarity_narrow.fits')

galah_data_input = '/home/klemen/data4_mount/'
galah_params = Table.read(galah_data_input+'sobject_iraf_52_reduced_20171111.fits')

sim_metrics = sim_res.colnames[1:]
scores = np.arange(len(sim_res))+1

params_joined = join(galah_params, sim_res, keys='sobject_id')
for i_b in [2,3]:
    for metric in sim_metrics:
        metric_values = params_joined[metric]
        if np.sum(np.isfinite(metric_values)) <= 0:
            continue
        snr_col = 'snr_c'+str(i_b)+'_iraf'
        snr_values = params_joined[snr_col]
        plt.scatter(snr_values, metric_values, lw=0, s=2)
        y_perc = np.percentile(metric_values, 500./len(metric_values)*100)
        plt.ylim(0, np.nanpercentile(metric_values, 92))
        plt.title(metric+' band:'+str(i_b))
        plt.axhline(y=y_perc, c='black')
        plt.xlim(0, np.nanpercentile(snr_values, 99.8))
        plt.xlabel(snr_col)
        plt.ylabel(metric)
        plt.show()
        # plt.savefig(metric+'_b'+str(i_b)+'.png', dpi=300)
        plt.close()

for metric in sim_metrics:
    print '-----------------------------'
    print metric
    if not np.isfinite(sim_res[metric]).all():
        continue
    idx_sort = np.argsort(sim_res[metric])
    sim_res[metric][idx_sort] = scores

    print '-- Best'
    s_best = sim_res[idx_sort]['sobject_id'][:100]
    print ','.join([str(a) for a in s_best])
    print '-- Worst'
    s_worst = sim_res[idx_sort]['sobject_id'][-100:]
    print ','.join([str(a) for a in s_worst])

# final
sim_res['final_score'] = np.nanmean(sim_res[sim_metrics].to_pandas().values, axis=1)
print '-- Final best'
idx_sort = np.argsort(sim_res['final_score'])
s_best = sim_res[idx_sort]['sobject_id'][:100]
print ','.join([str(a) for a in s_best])
print '-- Final worst'
s_worst = sim_res[idx_sort]['sobject_id'][-100:]
print ','.join([str(a) for a in s_worst])

# plt.scatter(params_joined['chebyshev'], params_joined['minkowski'])
# plt.show()
# plt.close()

print sim_res

