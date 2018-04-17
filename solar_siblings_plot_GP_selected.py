from os import chdir
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

bands = [1, 2, 3, 4]
b_suffx = ''.join([str(b) for b in bands])
input_dir = 'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp__norm-obs_lnlike-obs'
# input_dir = '/home/klemen/data4_mount/Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325__1_norm-obs/'
chdir(input_dir)

sim_file_root = 'solar_similarity_b'+b_suffx+'_gp'
sim_data = Table.read(sim_file_root+'.csv')
gp_data = Table.read('GP_fit_res.txt', format='ascii.csv')

suffix = '_ext0_dateall_offset'
solar_input_dir = '/home/klemen/data4_mount/Solar_data_dr53/'
solar_data = pd.read_csv(solar_input_dir + 'twilight_spectrum_galah'+suffix+'.txt', header=None, delimiter=' ', na_values='nan').values
solar_wvl = solar_data[:, 0]
solar_flx = solar_data[:, 1]

# get reference solar spectrum

plot_spectra = sim_data['sobject_id']
# plot_spectra = []

for s_b in bands:
    wvl = pd.read_csv('gp_median_wvl_b' + str(s_b) + '.csv', header=None, sep=',').values[0]
    flx = pd.read_csv('gp_median_flx_b' + str(s_b) + '.csv', header=None, sep=',').values
    flx_multi_all = list([])
    for s_id in plot_spectra:
        idx = np.where(sim_data['sobject_id'] == s_id)[0]
        flx_multi = gp_data[gp_data['sobject_id'] == s_id]['cont_norm_b' + str(s_b)][0]
        flx_vals = flx[idx, :][0]
        plt.plot(wvl, flx_vals, lw=1, alpha=0.75)#, c='blue')
        # plt.plot(wvl, flx_vals / flx_multi, lw=1, alpha=0.2, c='green')
        flx_multi_all.append(flx_multi)
    idx_sol_use = np.logical_and(solar_wvl >= wvl[0], solar_wvl <= wvl[-1])
    plt.plot(solar_wvl[idx_sol_use], solar_flx[idx_sol_use], lw=2, c='black')
    plt.show()
    plt.close()

    plt.hist(flx_multi_all, bins=75)
    plt.show()
    plt.close()
