from os import chdir
from solar_siblings_functions import *

suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'

# read Galah guess and/or cannon parameters
galah_params = Table.read(galah_data_input+'sobject_iraf_53_reduced_20180214.fits')

# distance/similarity measurements
chdir('Distances_Step1_p0_SNRsamples0_ext0')
evaluate_bands = [1] #list([1, 2, 3, 4])

final_selected_objects = {}

for i_b in evaluate_bands:
    print 'Band', i_b

    solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=1)  # every_nth_solar_pixel[i_b - 1])

    sim_res = np.loadtxt('solar_spectral_diff_b{:.0f}.csv'.format(i_b), delimiter=',')
    wvl_pos = sim_res[0, 1:]
    # find discontinues in wvl data
    insert_nans = list([])
    for i_p in range(len(wvl_pos)-1):
        if wvl_pos[i_p+1] - wvl_pos[i_p] > 0.1:
            insert_nans.append(i_p+1)

    median_diff = np.nanmedian(sim_res[1:, 1:], axis=0)  # remove first col (id info)
    std_diff = np.nanstd(sim_res[1:, 1:], axis=0)  # remove first col (id info)

    wvl_pos = np.insert(wvl_pos, insert_nans, np.full(len(insert_nans), np.nan))
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
    axSpectra.set_title('')
    plt.show()
    # plt.savefig(str(s_obj) + '_' + str(galah_object['snr_c2_guess'].data[0]) + bands_suffix + '.png', dpi=350)
    plt.close()
