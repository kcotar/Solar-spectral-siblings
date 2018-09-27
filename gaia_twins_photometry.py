import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, unique
from os import chdir
from scipy.interpolate import spline
import astropy.coordinates as coord
import astropy.units as un

galah_data_input = '/data4/cotar/'
isochrones_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_Solar/'
evoltracks_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_evolutionary_track/'

# data-table settings
data_date = '20180327'
cannon_param_file = 'sobject_iraf_iDR2_180325_cannon.fits'

cannon_data = Table.read(galah_data_input+cannon_param_file)
apass_data = Table.read(galah_data_input+'photometry/apass_dr53_20180327.csv')['sobject_id','Vmag','e_Vmag','Bmag','e_Bmag', 'gpmag','e_gpmag','u_e_gpmag','rpmag','e_rpmag','u_e_rpmag','ipmag']
wise_data = Table.read(galah_data_input+'photometry/wise_dr52_20171111.csv')['sobject_id','W1mag','W2mag','W3mag','W4mag','e_W1mag','e_W2mag','e_W3mag','e_W4mag']
tmass_data = Table.read(galah_data_input+'photometry/2mass_dr52_20171111.csv')['sobject_id','Jmag','Hmag','Kmag', 'e_Jmag','e_Hmag','e_Kmag']
wise_data = unique(wise_data, keys='sobject_id', keep='first')
tmass_data = unique(tmass_data, keys='sobject_id', keep='first')
apass_data = unique(apass_data, keys='sobject_id', keep='first')

# chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb')
# gp_res = Table.read('solar_similarity_b1234_gp.csv')
# print len(gp_res)
# # predetermined objects
# idx_run_params = 4.2
# solar_like_sobjects = gp_res['sobject_id']

idx_run_params = 8
params_str = ['5100_4.54_0.00', '5200_4.53_0.00', '5300_4.51_0.00', '5400_4.48_0.00',
              '5500_4.45_0.00', '5600_4.41_0.00', '5700_4.36_0.00', '5800_4.31_0.00',
              '5900_4.25_0.00', '6000_4.18_0.00'][idx_run_params]
multi_mag_thr = [5.7, 5.5, 5.3, 5.1,
                 4.9, 4.6, 4.3, 3.9,
                 3.6, 3.4][idx_run_params]

chdir(galah_data_input + 'Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_G20180327_C180325_refpar_'+params_str)
sel_txt = open('final_selection_0.10.txt', 'r')
solar_like_sobjects = sel_txt.read()
sel_txt.close()
solar_like_sobjects = [np.int64(sid) for sid in solar_like_sobjects.split(',')]
print len(solar_like_sobjects)

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]
# cannon_data = join(cannon_data, gp_res, join_type='left', keys='sobject_id')
# print len(cannon_data)
cannon_data = join(cannon_data, apass_data, join_type='left', keys='sobject_id')
cannon_data = join(cannon_data, wise_data, keys='sobject_id', join_type='left')
cannon_data = join(cannon_data, tmass_data, keys='sobject_id', join_type='left')
print len(cannon_data)


# absolute values in Vega magnitudes
# https://arxiv.org/pdf/1804.07788.pdf
B_s = 5.44
V_s = 4.81
R_s = 4.43
I_s = 4.10

# Gaia relations for Gmag
# https://www.aanda.org/articles/aa/pdf/forth/aa32756-18.pdf
G_1 = -0.01746 + 0.008092 * (V_s - I_s) -0.2810 * (V_s - I_s)**2 + 0.03655 * (V_s - I_s)**3 + V_s
G_2 = -0.02269 + 0.01784 * (V_s - R_s) -1.016 * (V_s - R_s)**2 + 0.2225 * (V_s - R_s)**3 + V_s
G_3 = -0.02907 -0.02385 * (B_s - V_s) -0.2297 * (B_s - V_s)**3 -0.001768 * (B_s - V_s)**3 + V_s
G_s = [G_1, G_2, G_3]
G_mean = np.mean(G_s)
print G_mean

# make positional plots, l/b and galaxy position based on gaia dr2
gaia_data = Table.read(galah_data_input+'sobject_iraf_53_gaia.fits')#['sobject_id','source_id','parallax','parallax_error','phot_bp_mean_mag','phot_rp_mean_mag','phot_g_mean_mag', 'pmra', 'pmdec','radial_velocity']
# gaia_data = gaia_data[gaia_data['angDist'] < 0.5]
gaia_data = unique(gaia_data, keys='sobject_id', keep='first')
# gaia_data['parallax'][gaia_data['parallax_error'] > gaia_data['parallax']] = np.nan
cannon_data = join(cannon_data, gaia_data, keys='sobject_id', join_type='left')
# cannon_data = cannon_data.filled(-1)
cannon_data['parsec'] = 1e3/cannon_data['parallax']
# print cannon_data['sobject_id', 'phot_g_mean_mag', 'parsec', 'parallax', 'parallax_error']

# simulate parallax and Gmag as would be observe by the Gaia spacecraft
parallax_sim = np.linspace(0.5, 15, 120)
Gmag_sim = G_mean + 2.5*np.log10(((1e3/parallax_sim)/10.)**2)

Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2) #- cannon_data['a_g_val']
Vmag_twin_abs = cannon_data['Vmag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
cannon_data['g_mean_mag_abs'] = Gmag_twin_abs


# idx_mark = np.in1d(cannon_data['sobject_id'], [150208003201286,150427004801275,160130006301234,160327004601337,160524006601258,160916001801263,161118004701221,161119002801292,170117003101044,170205005401120,170515003101036,170516002101273,171001001601082,171207003601278,180129003101184])
# idx_mark = np.in1d(cannon_data['sobject_id'], li_high_sid)
# idx_mark = Gmag_twin_abs < multi_mag_thr
idx_mark = cannon_data['parallax_error'] > cannon_data['parallax']*0.1

plt.scatter(cannon_data['phot_g_mean_mag'], cannon_data['Vmag'], lw=0, s=3)
plt.xlabel('Gaia G')
plt.ylabel('APASS V')
plt.savefig('gaia_apass.png', dpi=300)
plt.close()

plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], cannon_data['Bmag']-cannon_data['Vmag'], lw=0, s=3)
plt.xlabel('Gaia Gbp-Grp')
plt.ylabel('APASS B-V')
plt.savefig('gaia_apass_color.png', dpi=300)
plt.close()


plt.plot(Gmag_sim, parallax_sim, lw=1, color='C2', ls='--', label='Relation for the Sun')
plt.plot(Gmag_sim-0.75, parallax_sim, lw=1, color='black', ls='--', label='Sun twins binary')
plt.errorbar(cannon_data['phot_g_mean_mag'], cannon_data['parallax'], yerr=cannon_data['parallax_error'],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0, label='Gaia observations')
plt.errorbar(cannon_data['phot_g_mean_mag'][idx_mark], cannon_data['parallax'][idx_mark], yerr=cannon_data['parallax_error'][idx_mark],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='red', markeredgewidth=0, label='')
plt.xlabel('Gaia G mean magnitude')

plt.ylabel('Gaia parallax')
plt.xlim((9, 14))
plt.ylim((0.5, 13))
plt.grid(color='black', ls='--', alpha=0.2)
plt.legend()
plt.tight_layout()
plt.savefig('mag_parallax_gaia.png', dpi=250)
plt.close()

# load isochrones that will be added tot he plot
iso_data = Table.read(isochrones_dir+'isochrones_all.fits')
iso_data = iso_data[iso_data['Age'] == 4500000000.]
track_data = Table.read(evoltracks_dir+'isochrones_all.fits')
track_data['BP_RP'] = track_data['G_BPmag'] - track_data['G_RPmag']
# track_data = track_data[track_data['Age'] >= 1e9]


binaries = [140709003001184, 150408005301259, 150830004601175, 160326002101077, 160401004401051, 160420003301342,
            160424002101194, 160425001901273, 160426004501395, 160524004201306, 160530003901026, 160530003901054,
            160813002101314, 160813004101136, 161104002801361, 161116003801308, 161117004001161, 170107004201309,
            170112003601298, 170510007301226, 170514003301011, 170614004601055, 170710002701354, 150409002601317,
            150411004101331, 150606005901339, 160327006101355, 160330002601095, 160521004801353, 160530005501077,
            160531005101256, 160602001601307, 161008002501018, 161210004201315, 161211003101387, 161217002601138,
            161219005101228, 170117003101044, 170408004501048, 170514003001180, 170906003601357, 170909002601291,
            171207003601278, 171227003601367]
multiples = [140608002501303, 140808002701338, 140808003701104, 141102002701267, 150408005301116, 150413003601344,
             150413005101096, 160107003101157, 160125004501038, 160325002701048, 160401003901215, 160531006101153,
             161009005901171, 161118002601376, 161212002101397, 161217004101075, 161217004101234, 170121002801292,
             170507007801271, 170508002601001, 170508004801312, 170514002401099, 170514003301001, 170911003101383,
             140608002501303, 140808003701104, 150413003601344, 160125004501038, 160401003901215, 161009005901171,
             161118002601376, 161217004101234, 170507007801271, 170514002401099, 170911003101383]
inconclusive = [161009002601018]

c_m = ['Gmag', 'BP_RP']
# add evolutionary track
for mh_ini in np.unique(track_data['MHini']):
    for M_star in [0.99, 1.00, 1.01]:
        track_mag = list([])
        track_mag2 = list([])

        for i_age in np.unique(track_data['Age']):
            track_data_age = track_data[np.logical_and(track_data['Age'] == i_age, track_data['MHini'] == mh_ini)]
            idx_m = np.argmin(np.abs(track_data_age['Mini']-M_star))
            # track_mag.append([np.interp(np.array([M_star]), track_data_age['Mini'][idx_m-3:idx_m+3], track_data_age[c][idx_m-3:idx_m+3], 3) for c in c_m])
            track_mag2.append([np.polyval(np.polyfit(track_data_age['Mini'][idx_m-3:idx_m+3]-M_star, track_data_age[c][idx_m-3:idx_m+3], 3), np.array([0])) for c in c_m])
        # track_mag = np.array(track_mag)
        track_mag2 = np.array(track_mag2)
        # plt.plot(track_mag[:, 1], track_mag[:, 0])
        plt.plot(track_mag2[:, 1], track_mag2[:, 0], label='[M/H]={:.2f}  Mini={:.2f}'.format(mh_ini, M_star))

    idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
    idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
    idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)

    plt.scatter(cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, label='', c='black')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=8, c='C3', label='Discarded')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=8, c='C0', label='Binary')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=8, c='C1', label='Triple')
    plt.ylabel(r'${\it Gaia}$ G absolute magnitude')
    plt.xlabel(r'G$_{bp}$ - G$_{rp}$')
    plt.xlim(0.7, 1.35)

    # plt.scatter(cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'] - cannon_data['e_bp_min_rp_val'], Gmag_twin_abs, s=4, lw=0, label='', c='black')
    # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'] - cannon_data['e_bp_min_rp_val'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=8, c='C3', label='Discarded')
    # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'] - cannon_data['e_bp_min_rp_val'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=8, c='C0', label='Binary')
    # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'] - cannon_data['e_bp_min_rp_val'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=8, c='C1', label='Triple')
    # plt.ylabel(r'${\it Gaia}$ G absolute magnitude - A$_G$')
    # plt.xlabel(r'G$_{bp}$ - G$_{rp}$ - E(Bp-Rp)')
    # plt.xlim(0.1, 1.35)

    plt.ylim(2, 6)
    plt.axhline(multi_mag_thr, ls='--', c='black', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.legend()
    # plt.show()
    # plt.savefig('mag_hr_gaia_bin-multi_evol_mh{:.2f}_A_E.png'.format(mh_ini), dpi=300)
    plt.savefig('mag_hr_gaia_bin-multi_evol_mh{:.2f}.png'.format(mh_ini), dpi=300)
    plt.close()

raise SystemExit

for c_col in ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']:

    # add isochrone
    for i_feh in np.unique(iso_data['MHini']):
        idx_iso = np.logical_and(iso_data['MHini'] == i_feh, iso_data['Mini'] - iso_data['Mass'] < 0.1)
        plt.plot(iso_data['G_BPmag'][idx_iso] - iso_data['G_RPmag'][idx_iso], iso_data['Gmag'][idx_iso],
                 # lw=0.75, label='[M/H] = {:.1f}'.format(i_feh))
                 lw=0.5, label='', c='black')

    # plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, c=cannon_data[c_col], label='')
    # plt.colorbar()
    # plt.scatter((cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=7, label='', c='C3')
    # plt.ylabel(r'{\it Gaia} G absolute magnitude')
    # plt.xlabel(r'G$_{bp}$ - G$_{rp}$')
    # plt.ylim(2, 6)
    # plt.xlim(0.7, 1.35)

    binaries = [140709003001184,150408005301259,150830004601175,160326002101077,160401004401051,160420003301342,160424002101194,160425001901273,160426004501395,160524004201306,160530003901026,160530003901054,160813002101314,160813004101136,161104002801361,161116003801308,161117004001161,170107004201309,170112003601298,170510007301226,170514003301011,170614004601055,170710002701354,150409002601317,150411004101331,150606005901339,160327006101355,160330002601095,160521004801353,160530005501077,160531005101256,160602001601307,161008002501018,161210004201315,161211003101387,161217002601138,161219005101228,170117003101044,170408004501048,170514003001180,170906003601357,170909002601291,171207003601278,171227003601367]
    multiples = [140608002501303,140808002701338,140808003701104,141102002701267,150408005301116,150413003601344,150413005101096,160107003101157,160125004501038,160325002701048,160401003901215,160531006101153,161009005901171,161118002601376,161212002101397,161217004101075,161217004101234,170121002801292,170507007801271,170508002601001,170508004801312,170514002401099,170514003301001,170911003101383,140608002501303,140808003701104,150413003601344,160125004501038,160401003901215,161009005901171,161118002601376,161217004101234,170507007801271,170514002401099,170911003101383]
    inconclusive = [161009002601018]
    idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
    idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
    idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)
    plt.scatter(cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, label='', c='black')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=8, c='C3', label='Discarded')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=8, c='C0', label='Binary')
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=8, c='C1', label='Triple')
    # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_inc], Gmag_twin_abs[idx_inc], lw=0, s=7, c='C1', label='Inconclusive')
    plt.ylabel(r'${\it Gaia}$ G absolute magnitude')
    plt.xlabel(r'G$_{bp}$ - G$_{rp}$')
    plt.ylim(2, 6)
    plt.xlim(0.7, 1.35)

    plt.axhline(multi_mag_thr, ls='--', c='black', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.legend()
    # plt.savefig('mag_hr_gaia_' + c_col + '.png', dpi=300)
    plt.savefig('mag_hr_gaia_bin-multi.png', dpi=300)
    plt.close()

idx_mag_multiple = Gmag_twin_abs < multi_mag_thr

for g_p in ['astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']:
    x_range = (np.nanpercentile(cannon_data[g_p], 1), np.nanpercentile(cannon_data[g_p], 99))
    plt.hist(cannon_data[~idx_mag_multiple][g_p], range=x_range, bins=250, color='blue', alpha=0.5)
    plt.hist(cannon_data[idx_mag_multiple][g_p], range=x_range, bins=250, color='red', alpha=0.5)
    plt.gca().set_yscale('log')
    plt.savefig('Gp_'+g_p+'.png', dpi=250)
    plt.close()

print ','.join([str(s) for s in cannon_data[idx_mag_multiple]['sobject_id']])


plt.figure(1, figsize=(8, 6))
ax = list([])
ax.append(plt.axes([0.08, 0.27, 0.90, 0.71]))
ax.append(plt.axes([0.08, 0.07, 0.90, 0.2]))

cannon_data = cannon_data[cannon_data['parallax_error'] < cannon_data['parallax']*0.1]
Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
idx_mag_multiple = Gmag_twin_abs < multi_mag_thr

# all photometric data together plot - and their distribution
to_abs_mag = (-2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)).reshape(-1, 1)
c_plot = ['Bmag','Vmag','gpmag','rpmag','ipmag','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','Jmag','Hmag','Kmag','W1mag','W2mag','W3mag','W4mag']
x_pos = np.repeat([np.arange(len(c_plot))], len(cannon_data), axis=0)

phot_data_multi = cannon_data[c_plot][idx_mag_multiple].to_pandas().values + to_abs_mag[idx_mag_multiple]
median_1 = np.nanmedian(phot_data_multi, axis=0)
# prepare data to be ploted as a violin
phot_data_multi_pos = np.arange(len(c_plot))*3
plot_d_list = list([])
for i_c in range(len(c_plot)):
    plot_data = phot_data_multi[:, i_c]
    plot_data = plot_data[np.isfinite(plot_data)]
    plot_d_list.append(plot_data)

violin_parts = ax[0].violinplot(plot_d_list, phot_data_multi_pos, showmeans=False, showextrema=False, showmedians=False, widths=1)
for pc in violin_parts['bodies']:
    pc.set_facecolor('C1')
    pc.set_edgecolor('black')
    pc.set_alpha(0.75)

phot_data_multi = cannon_data[c_plot][~idx_mag_multiple].to_pandas().values + to_abs_mag[~idx_mag_multiple]
median_2 = np.nanmedian(phot_data_multi, axis=0)
# prepare data to be ploted as a violin
phot_data_multi_pos = np.arange(len(c_plot))*3 + 1
plot_d_list = list([])
for i_c in range(len(c_plot)):
    plot_data = phot_data_multi[:, i_c]
    plot_data = plot_data[np.isfinite(plot_data)]
    plot_d_list.append(plot_data)

violin_parts = ax[0].violinplot(plot_d_list, phot_data_multi_pos, showmeans=False, showextrema=False, showmedians=False, widths=1)
for pc in violin_parts['bodies']:
    pc.set_facecolor('C2')
    pc.set_edgecolor('black')
    pc.set_alpha(0.75)

import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red')
ax[0].legend([mpatches.Patch(color='C1'), mpatches.Patch(color='C2')],
           ['Multiple', 'Single'], loc=2)
ax[0].set(ylim=(0, 6.5), xlim=(-1, 44), ylabel='Magnitude',
          xticks=np.arange(len(c_plot))*3 + 0.5, xticklabels=['' for i in range(len(c_col))])
ax[0].invert_yaxis()
ax[0].grid(ls='--', color='black', alpha=0.2)

ax[1].scatter(np.arange(len(c_plot))*3 + 0.5, median_2-median_1, lw=0, s=10, c='black')
ax[1].set(ylim=(0.8, 1.0), xlim=(-1, 44), ylabel='Difference', xlabel='Photometric filter',
          xticks=np.arange(len(c_plot)) * 3 + 0.5,
          xticklabels=['B', 'V', "g'", "r'", "i'", 'G', r'G$_{BP}$', r'G$_{RP}$', 'J', 'H', 'K', 'W1', 'W2', 'W3', 'W4'],
          yticks=[0.85, 0.95])
ax[1].grid(ls='--', color='black', alpha=0.2)

plt.savefig('multi_mag_plot.png', dpi=300)
plt.close()
