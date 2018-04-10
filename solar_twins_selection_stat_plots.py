import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join
from os import chdir

galah_data_input = '/home/klemen/data4_mount/'

# data-table settings
data_date = '20180327'
cannon_param_file = 'sobject_iraf_iDR2_180325_cannon.fits'

cannon_data = Table.read(galah_data_input+cannon_param_file)

chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_comb')
gp_res = Table.read('solar_similarity_b1234_gp.csv')

# predetermined objects
solar_like_sobjects = gp_res['sobject_id']

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]
cannon_data = join(cannon_data, gp_res, join_type='left')

# make histogram plots for parameters
plot_col = ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']
solar_val = [5771, 4.44, 0.02]
idx_valid_param = cannon_data['flag_cannon'] == 0
for i_p in range(len(plot_col)):
    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plot_vals_median = np.median(plot_vals)
    plt.hist(plot_vals, bins=75, color='black', alpha=0.33)
    plt.axvline(x=solar_val[i_p], color='black', ls='--', lw=2.5)
    plt.axvline(x=plot_vals_median, color='red', ls='--', lw=2.5)
    plt.title('All: {:.0f}    Unflagged: {:.0f}    Median: {:.2f}    Difference: {:.2f}'.format(len(solar_like_sobjects), len(plot_vals), plot_vals_median, plot_vals_median-solar_val[i_p]))
    plt.tight_layout()
    # plt.show()
    plt.savefig('solar_twins_like_'+plot_col[i_p]+'.png', dpi=300)
    plt.close()

    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plt.scatter(plot_vals, cannon_data['canberra'][idx_valid_param], s=4, lw=0, c='black')
    plt.tight_layout()
    # plt.show()
    plt.savefig('solar_twins_like_' + plot_col[i_p] + '_canberra.png', dpi=300)
    plt.close()



# make histogram plots for parameters
plot_col = [col for col in cannon_data.colnames if '_abund_cannon' in col and 'e_' not in col and 'flag_' not in col and len(col.split('_')) < 4]
solar_val = np.zeros(len(plot_col))
for i_p in range(len(plot_col)):
    idx_valid_param = cannon_data['flag_'+plot_col[i_p]] == 0
    print plot_col[i_p]+': '
    if 'La' in plot_col[i_p]:
        idx_Li = cannon_data[plot_col[i_p]] > 0.1
        sob_Li = cannon_data[np.logical_and(idx_Li, idx_valid_param)]['sobject_id']
        print ','.join([str(s) for s in sob_Li])

    # create subsets by similarity levels
    for sim_val in np.arange(0.002, 0.007, 0.001):
        idx_sim = np.logical_and(gp_res['canberra'] > sim_val, gp_res['canberra'] <= sim_val+0.001)
        if np.sum(idx_sim) == 0:
            continue
        cannon_data_sub = cannon_data[np.in1d(cannon_data['sobject_id'], gp_res['sobject_id'][idx_sim])]

    idx_valid_param = cannon_data['flag_'+plot_col[i_p]] == 0
    if np.sum(idx_valid_param) < 3:
        continue
    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plot_vals_median = np.median(plot_vals)
    plot_vals_std = 2.*np.std(plot_vals)
    # idx_out = np.abs(plot_vals - plot_vals_median) > plot_vals_std
    idx_out = np.abs(plot_vals - plot_vals_median) < 0.05
    # if np.sum(idx_out) > 0:
        # print ','.join([str(so) for so in cannon_data['sobject_id'][idx_valid_param][idx_out]])
    plt.hist(plot_vals, bins=100, range=(-1.2, 1.2), color='black', alpha=0.33)
    # plt.axvline(x=solar_val[i_p], color='black', ls='--', lw=2.5)
    plt.axvline(x=plot_vals_median, color='red', ls='--', lw=1.5, alpha=0.8)
    plt.axvline(x=plot_vals_median-plot_vals_std, color='red', ls='--', lw=1.5, alpha=0.33)
    plt.axvline(x=plot_vals_median+plot_vals_std, color='red', ls='--', lw=1.5, alpha=0.33)

    plt.title('All: {:.0f}    Unflagged: {:.0f}    Median: {:.2f}    Difference: {:.2f}'.format(len(solar_like_sobjects), len(plot_vals), plot_vals_median, plot_vals_median-solar_val[i_p]))
    # plt.tight_layout()
    # plt.show()
    # plt.savefig('GP_solar_twins_like_'+plot_col[i_p]+'_'+str(sim_val)+'.png', dpi=300)
    plt.savefig('solar_twins_like_'+plot_col[i_p]+'.png', dpi=300)
    plt.close()
