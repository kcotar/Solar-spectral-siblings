from os import chdir
from solar_siblings_functions import *

suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'

# read Galah guess and/or cannon parameters
galah_params = Table.read(galah_data_input+'sobject_iraf_53_reduced_20180214.fits')

# distance/similarity measurements
chdir('Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_2')
evaluate_bands = list([1,2,3,4])

final_selected_objects = {}

for i_b in evaluate_bands:
    print 'Band', i_b

    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=1)  # every_nth_solar_pixel[i_b - 1])

    sim_res = pd.read_csv('solar_spectral_diff_b{:.0f}.csv'.format(i_b), header=None, sep=',')
    s_ids_diff = sim_res[0][1:].values
    sim_res = sim_res.values
    wvl_pos = sim_res[0, 1:]
    # find discontinues in wvl data
    insert_nans = list([])
    for i_p in range(len(wvl_pos)-1):
        if wvl_pos[i_p+1] - wvl_pos[i_p] > 0.1:
            insert_nans.append(i_p+1)
    wvl_pos = np.insert(wvl_pos, insert_nans, np.full(len(insert_nans), np.nan))

    # get snr values for observed sobject_ids
    idx_match = np.in1d(galah_params['sobject_id'], s_ids_diff)
    s_ids_snr = galah_params['snr_c'+str(i_b)+'_guess'][idx_match]

    # determine snr ranges to be evaluated
    snr_ranges = np.arange(10, 200, 10)

    # iterate over all possible snr ranges
    for i_snr in range(len(snr_ranges)-1):
        print ' SNR:', snr_ranges[i_snr]
        idx_snr_match = np.logical_and(s_ids_snr >= snr_ranges[i_snr], s_ids_snr < snr_ranges[i_snr+1])
        idx_snr_match = np.where(idx_snr_match)[0] + 1  # +1 as first row was omitted when reading sobject ids from csv
        n_per_snr = len(idx_snr_match)

        if n_per_snr < 10:
            print '  Low per SNR'
            continue

        median_diff = np.nanmedian(sim_res[idx_snr_match, 1:], axis=0)  # remove first col (id info)
        std_diff = np.nanstd(sim_res[idx_snr_match, 1:], axis=0)  # remove first col (id info)

        median_diff = np.insert(median_diff, insert_nans, np.full(len(insert_nans), np.nan))
        std_diff = np.insert(std_diff, insert_nans, np.full(len(insert_nans), np.nan))

        idx_band_mask = get_band_mask(solar_wvl, i_b)
        idx_lines_mask = get_linelist_mask(solar_wvl)

        x_range = (solar_wvl[idx_band_mask][0]-2, solar_wvl[idx_band_mask][-1]+2)
        plt.figure(1, figsize=(12, 7))
        axSpectra = plt.axes([0.05, 0.3, 0.92, 0.65])
        axDiff = plt.axes([0.05, 0.05, 0.92, 0.20])
        axSpectra.plot(solar_wvl[idx_band_mask], solar_flx[idx_band_mask], c='black', lw=0.5)
        axDiff.axhline(y=0, c='black', lw=0.5)
        # axDiff.errorbar(wvl_pos, median_diff, yerr=std_diff, c='blue', lw=0.5)
        axDiff.plot(wvl_pos, median_diff, c='blue', lw=0.5)
        axSpectra.set(ylim=(0.3, 1.1), xlim=x_range)
        axDiff.set(ylim=(-0.05, 0.05), xlim=x_range)
        axSpectra.set_title('Number of objects: {:.0f}'.format(n_per_snr))
        #plt.show()
        plt.savefig('diff_b' + str(i_b) + '_snr{:03.0f}'.format(snr_ranges[i_snr]) + '.png', dpi=350)
        plt.close()
