from Mass_distribution_binaries_functions import *


# data-table settings
data_date = '20180327'
cannon_param_file = 'sobject_iraf_iDR2_180325_cannon.fits'
cannon_data = Table.read(galah_data_input+cannon_param_file)

sim_res_subdir = 'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb'
gp_res = Table.read(sim_res_subdir+'/solar_similarity_b1234_gp.csv')
gaia_data = Table.read(galah_data_input+'sobject_iraf_53_gaia.fits')

chdir('Binary_candidates_spectra_fit_MS')

# predetermined objects
ref_like_objects = gp_res['sobject_id']

cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], ref_like_objects)]
cannon_data = join(cannon_data, gaia_data, join_type='left', keys='sobject_id')
print len(cannon_data)

# Absolute magnitude of the reference star
Teff_ref = np.median(cannon_data['Teff_cannon'])
Feh_ref = np.median(cannon_data['Fe_H_cannon'])
Logg_ref = np.median(cannon_data['Logg_cannon'])
print Teff_ref
# predetermined values
G_abs_ref = 4.679
star_age = 4.5e9
star_mh = 0

# Absolute magnitude of observed stars
Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
cannon_data['g_mean_mag_abs'] = Gmag_twin_abs

# objects to analyze - possible binary stars
Gmag_thr = 4.2
obj_analyze = cannon_data[Gmag_twin_abs < Gmag_thr]

iso_data = select_isochrone(isochrones_data, star_mh, star_age)
# iso_data_2 = select_isochrone(isochrones_data, star_mh, 4.5e9)
# plt.plot(iso_data['G_BPmag']-iso_data['G_RPmag'], iso_data['Gmag'])
# plt.plot(iso_data['teff'], iso_data['Gmag'])
# plt.axhline(y=iso_data['Gmag'][np.argmax(iso_data['teff'])])
# plt.plot(iso_data_2['teff'], iso_data_2['Gmag'])
# plt.show()
# plt.close()

print 'Mag Sun - teff cannon:', get_gmag_from_teff_MS(iso_data, Teff_ref)
print 'Mag Sun - teff real:', get_gmag_from_teff_MS(iso_data, 5777)
print 'Mag Sun - UVB:', G_abs_ref

print '==========='

res_final_table = Table(names=['sobject_id', 'teff_1', 'teff_2', 'q', 'Gmag_iso1', 'Gmag_iso2', 'Gmag_iso_comb'],
                        dtype=['int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64'])

suffix = '_flx_mag'
res_fits = 'teff_fit_res'+suffix+'.fits'
if path.isfile(res_fits):
    print 'Reading previous results'
    res_final_table = Table.read(res_fits)

for obj in obj_analyze:
    print 'Object:', obj['sobject_id']
    # TODO: cannon flag check
    G_abs_cannon = get_gmag_from_teff_MS(iso_data, obj['Teff_cannon'])

    # load needed data
    obj_s_id = obj['sobject_id']
    flx, flx_s, wvl = get_spectra_complete(obj_s_id)

    # check if results were already obtained for the selected sobject_id
    if np.sum(res_final_table['sobject_id'] == obj_s_id) == 0:
        # VERSION 2
        print 'Suffix:', suffix
        time_s = time()
        sampler, fit_res, fit_prob = fit_spectra_and_teff(iso_data, obj, flx, flx_s, wvl,
                                                          G_abs_ref, obj['g_mean_mag_abs'],
                                                          nwalkers=100, n_threds=5, n_steps_1=150, n_steps_2=500, suffix=suffix)
        time_min = (time() - time_s)/60.
        print 'Total fit time: {:.1f} min'.format(time_min)

        # # VERSION 1
        # # use isochrones to determine Teff1 and Teff2 of the investigated star
        # sampler, fit_res, fit_prob = fit_teff_values(iso_data, G_abs_cannon, obj['g_mean_mag_abs'], Teff_ref,
        #                                              nwalkers=100, n_threds=50, n_steps=50)

        sampler_chain_vals = sampler.flatchain
        kernel_fit = get_distribution_peaks(sampler_chain_vals, plot_ref=str(obj_s_id), plot=False)
        # kernel_fit = np.median(sampler_chain_vals, axis=0)

        gmag_iso_1 = get_gmag_from_teff_MS(iso_data, kernel_fit[0])
        gmag_iso_2 = get_gmag_from_teff_MS(iso_data, kernel_fit[1])
        gmag_iso_final = add_mag(gmag_iso_1, gmag_iso_2)
        mass_ratio = get_mass_ration(iso_data, kernel_fit[0], kernel_fit[1])

        # plot lnprob and walkers
        plot_walkers(sampler.lnprobability, str(obj_s_id) + '_lnprob' + suffix + '_2.png',
                     title='Mag obs: {:.3f}   Mag final: {:.3f}   T1:{:.1f}   T2:{:.1f}'.format(obj['g_mean_mag_abs'], gmag_iso_final, kernel_fit[0], kernel_fit[1]))

        print ' Teff final:', kernel_fit
        print ' Mag final:', gmag_iso_1, gmag_iso_2, gmag_iso_final
        print ' Mag start:', obj['g_mean_mag_abs']
        print ' Mass ratio:', mass_ratio

        res_final_table.add_row([obj['sobject_id'], kernel_fit[0], kernel_fit[1], mass_ratio, gmag_iso_1, gmag_iso_2, gmag_iso_final])

        # output complete final results
        res_final_table.write(res_fits, overwrite=True)

    # PLOT RESULTS
    # presume that results for this object already exist, read and use them

    # output spectral results
    plot_spectra_comparison(obj, res_final_table, flx, wvl, path=str(obj_s_id) + '_spectra_comp.png')
    # output combination of absorption lines for one element
    combine_abs_lines_velocity_space(flx, wvl, 'Fe', path=str(obj_s_id) + '_lines_Fe.png')
    combine_abs_lines_velocity_space(flx, wvl, 'Ti', path=str(obj_s_id) + '_lines_Ti.png')
    combine_abs_lines_velocity_space(flx, wvl, 'V', path=str(obj_s_id) + '_lines_V.png')

    # END of processing for this objects
    print ''

# FINAL RESULTS
plt.hist(res_final_table['q'], range=(0., 1.), bins=50)
plt.title('mass ratio q histogram')
plt.savefig('q_hist'+suffix+'.png', dpi=250)
plt.close()
