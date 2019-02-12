import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.table import Table, join, unique
from os import chdir
import astropy.coordinates as coord
import astropy.units as un
from dustmaps.bayestar import BayestarWebQuery
from gaia_twins_photometry_functions import fit_MS_lin_line

gmag_MS_ebv = [6.46440852913791, 6.406135638759127, 6.326063057758769, 6.270662313171322, 6.186901128845153, 6.112655265836556, 6.047964331920713, 5.975597998849081, 5.864775024621654, 5.752935773085529, 5.70615132715316, 5.607932103569, 5.510605841548255, 5.394819317583499, 5.320540582224229, 5.206141954573217, 5.105337368096131, 4.997662315562364, 4.882464233710049, 4.756352703758858, 4.644881458729178, 4.527907502273222, 4.408222001484199, 4.2871798967317405, 4.117467298488372, 3.994238242422548, 3.861777335534732, 3.72487273109238, 3.633423030730339, 3.477971871745545, 3.3525566275598284, 3.1855673220231195, 2.9850365948134816, 2.830982397100216]
bprp_MS_ebv = [1.2908581206431133, 1.2588210724583755, 1.2268114976583142, 1.2073005303122892, 1.1797431892275174, 1.1542746627289153, 1.1258597376937827, 1.0968862630807177, 1.0720300797647617, 1.0410309764107133, 1.02608802792935, 0.9972604104823493, 0.9779415187202556, 0.9411407394522913, 0.942366454124274, 0.9106812550604175, 0.8997248316011284, 0.877800406566954, 0.8593017040995061, 0.8325740629627321, 0.8233155569426165, 0.7927167614889337, 0.7802655762838091, 0.7615338470099067, 0.7389512071183755, 0.7321674777076357, 0.7169088241276089, 0.6948854046524628, 0.6878793860732788, 0.6748835060196896, 0.63182029448188, 0.6369998924908571, 0.6304756782387, 0.6046636421838443]

plt.rcParams['font.size'] = 15

galah_data_input = '/shared/ebla/cotar/'
results_data_input = '/shared/data-camelot/cotar/_Multiples_binaries_results_iDR3/'

isochrones_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_Solar/'
evoltracks_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_evolutionary_track/'

# data-table settings
# cannon_param_file = 'galah_cannon_gaia_photometry_20180327_ebv-corr.fits'
# cannon_version = 'Cannon iDR2'
# suffix = '_ebv_c2'

cannon_param_file = 'galah_cannon_DR3_gaia_photometry_20181221_ebv-corr.fits'
cannon_version = 'SME iDR3'
suffix = '_ebv_c3'

cannon_data = Table.read(galah_data_input+cannon_param_file)

argv = sys.argv
if len(argv) >= 2:
    envel_val = int(argv[1])
else:
    envel_val = 10

if len(argv) >= 3:
    idx_run_params = int(argv[2])
    params_str = ['5100_4.55_0.00', '5200_4.53_0.00', '5300_4.51_0.00', '5400_4.48_0.00',
                  '5500_4.46_0.00', '5600_4.43_0.00', '5700_4.40_0.00', '5800_4.37_0.00',
                  '5900_4.34_0.00', '6000_4.30_0.00'][idx_run_params]
    multi_mag_thr = [5.7, 5.5, 5.3, 5.2,
                     4.9, 4.6, 4.3, 3.9,
                     3.6, 3.4][idx_run_params]
    suffix += '_' + params_str.split('_')[0]
    chdir(
        results_data_input + 'Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_G20180327_C181221_withH_refpar_' + params_str)
else:
    chdir(results_data_input + 'Distances_Step1_p0_SNRsamples0_ext0_oklinesonly_G20180327_C181221_withH')
    multi_mag_thr = 4.2
    suffix += ''


envel_str = '{:02.0f}'.format(envel_val)
sel_txt = open('final_selection_' + envel_str + '_envelope.txt', 'r')
solar_like_sobjects = sel_txt.read()
sel_txt.close()
solar_like_sobjects = [np.int64(sid) for sid in solar_like_sobjects.split(',')]
suffix += '_' + envel_str
print suffix

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]

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
# print G_mean, np.std(G_s), G_s
G_bp_sun = -0.05204 + 0.4830 * (V_s - I_s) - 0.2001 * (V_s - I_s)**2 + 0.02186 * (V_s - I_s)**3 + V_s
G_rp_sun = 0.0002428 - 0.8675 * (V_s - I_s) - 0.02866 * (V_s - I_s)**2 + V_s

sun_mag_list = np.array([5.44, 4.81, 5.23, 4.53, 4.19, G_mean, G_bp_sun, G_rp_sun, 3.67, 3.32, 3.27, 3.26, 3.28, 3.26, 3.27])

# simulate parallax and Gmag as would be observe by the Gaia spacecraft
parallax_sim = np.linspace(0.5, 15, 120)
Gmag_sim = G_mean + 2.5*np.log10(((1e3/parallax_sim)/10.)**2)

# Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((cannon_data['r_est'])/10.)**2)
Gmag_bp_rp = cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag']
# idx_mark = cannon_data['parallax_error'] > cannon_data['parallax']*0.2
idx_mark = cannon_data['ruwe'] > 1.4
# idx_mark = ~np.isfinite(cannon_data['parallax'])

# fit simple linear MS ---model to the median distribution of magnitudes
bp_rp_eval_complete = np.arange(0., 2., 0.1)
bin_bp_rb = 0.25
gg_fit = fit_MS_lin_line(Gmag_twin_abs, Gmag_bp_rp, path='MS_binary_line' + suffix + '.png', d_above=bin_bp_rb)

idx_above = Gmag_twin_abs <= (gg_fit(Gmag_bp_rp) - bin_bp_rb)
print 'Total:', len(solar_like_sobjects), 'above:', np.sum(idx_above)
print ','.join([str(s_sid) for s_sid in cannon_data['sobject_id'][idx_above]])

plt.figure(figsize=(7, 5.5))
plt.plot(Gmag_sim, parallax_sim, lw=1, color='black', ls='--', label='Relation for the Sun', alpha=0.8)
plt.plot(Gmag_sim-0.75, parallax_sim, lw=1, color='black', ls='-.', label='Sun twin binary', alpha=0.8)
plt.plot(Gmag_sim-1.2, parallax_sim, lw=1, color='black', ls=':', label='Sun twin triple', alpha=0.8)
plt.errorbar(cannon_data['phot_g_mean_mag'][~idx_mark], cannon_data['parallax'][~idx_mark], yerr=cannon_data['parallax_error'][~idx_mark],
             fmt='.', ms=5, elinewidth=1, alpha=0.9, color='C0', markeredgewidth=0, label='Gaia observations')
plt.errorbar(cannon_data['phot_g_mean_mag'][idx_mark], cannon_data['parallax'][idx_mark], yerr=cannon_data['parallax_error'][idx_mark],
             fmt='.', ms=3, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0, label='', marker='X')
plt.xlabel(u'${\it Gaia}$ G mean magnitude')
if 'ebv' in suffix:
    plt.xlabel(u'Extinction corrected ${\it Gaia}$ G mean magnitude')
plt.ylabel(u'${\it Gaia}$ parallax [mas]')
plt.xlim((9.5, 14))
plt.ylim((0.5, 11))
plt.xticks([10,11,12,13], ['10','11','12','13'])
plt.grid(color='black', ls='--', alpha=0.2)
plt.legend()
plt.tight_layout()
plt.savefig('mag_parallax_gaia'+suffix+'.png', dpi=250)
plt.close()

# load isochrones that will be added tot he plot
iso_data = Table.read(isochrones_dir+'isochrones_all.fits')
iso_data = iso_data[iso_data['Age'] == 4500000000.]
track_data = Table.read(evoltracks_dir+'isochrones_all.fits')
track_data['BP_RP'] = track_data['G_BPmag'] - track_data['G_RPmag']
track_data = track_data[(track_data['Mini'] - track_data['Mass'])/track_data['Mini'] < 0.3]

# >1
single_plus = [140711002901268,140814003301261,150408005301259,150706005901329,150828004701312,150830005101144,160420005301119,160531001601046,160816004201351,160817003101363,170515003601030,171102003301124,171227001601269,140807005001072,150108002201366,150330002601242,150411006101106,150705005401377,150828003701040,151230003201236,160107002601011,160326002101077,160327004101343,160327004601337,160328004201051,160401004401168,160403003601392,160420004301224,160426004501395,160524004201306,160530005001359,160813002101314,160817003101276,161007002801208,161007003301032,170114004101216,170407002601374,170414004101039,170515003101035,170603003101048,170614004601061,170713004601335,170805003601381,171001001601097,171102004501327,171106001901092]
# 2
binaries = [140806002301027,150408005301116,150409005601390,150411004101331,150703005101329,150830005101021,151111002101059,160327006101355,160330002101179,160330102801144,160331004301386,160424004201264,160531005101256,161118004701364,161217004101075,161219005101228,170108004601305,170115001601273,170117003101267,170408004501048,170408005501159,170517001801023,170601003601121,170711005101185,170801002801082,170801004601371,170829001901315,170909002601291,140413004401324,170112002601348,170418002701078,170905003101295,171001001601037]
# >2
binaries_plus = [140805003601343,140811004501343,151230003201396,160330002601095,160402006101178,160425003101385,161006003901387,161008002501018,161217002601138,170416004801356,170508004801312,170514003001180,170711004501329,170712004201047]
# 3
multiples = [140314004401277,140805003101303,150408004101169,151111002101116,160401004401123,170514002401099]
# i
inconclusive = [150330002601306]


c_m = ['Gmag', 'BP_RP']
i_ls = ['--', '-.', ':']
# add evolutionary track
for mh_ini in [0.0]: #np.unique(track_data['MHini']):
    print 'mh_ini:', mh_ini
    plt.figure(figsize=(7, 5.5))

    for i_a, M_star in enumerate([0.97, 1.00, 1.03]):
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
        str_mini = u'{:.2f}'.format(M_star)
        plt.plot(track_mag2[:, 1], track_mag2[:, 0], label=u'M$_{ini}$ = '+str_mini+u'M$_{\odot}$', ls=i_ls[i_a], c='black', alpha=0.7)

    idx_sin_p = np.in1d(cannon_data['sobject_id'], single_plus)
    idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
    idx_bin_p = np.in1d(cannon_data['sobject_id'], binaries_plus)
    idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
    idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)

    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[~idx_mark], Gmag_twin_abs[~idx_mark], s=7, lw=0, label='', c='C0', alpha=0.8)
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], s=10, lw=0, label='RUWE > 1.4', c='black', marker='X', alpha=0.8)

    plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete), c="C2", alpha=0.8, label='MS fit')
    plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete) - bin_bp_rb, ls='--', c="C2", alpha=0.8, label='Binaries limit')

    if 'ebv' in suffix:
        plt.ylabel(u'M$_G$ - A$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$ - A')
        plt.xlim(0.68, 1.)
        plt.xticks([0.7, 0.8, 0.9, 1.0], ['0.7', '0.8', '0.9', '1.0'])
        plt.ylim(2.2, 5.5)
    else:
        plt.ylabel(u'M$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$')
        plt.xlim(0.69, 1.1)
        plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1], ['0.7', '0.8', '0.9', '1.0', '1.1'])
        plt.ylim(2.2, 5.5)

    # plt.axhline(multi_mag_thr, ls='--', c='C2', alpha=0.8)
    plt.gca().invert_yaxis()
    plt.legend()
    # plt.legend(loc=2)
    plt.grid(ls='--', alpha=0.3)
    plt.tight_layout()
    # plt.show()
    plt.savefig('mag_hr_gaia_bin-multi_evol_mh{:.2f}_nores'.format(mh_ini)+suffix+'.png', dpi=300)
    plt.close()

l_meh = ['--', '-.', ':']
v_meh = [-0.49969851157094075, 0, 0.49994017852310874]
# print np.unique(iso_data['MHini'])
# add isochrones
plt.figure(figsize=(7, 5.5))
for i_meh, meh in enumerate(v_meh):
    idx_iso = np.logical_and(iso_data['MHini'] == meh, iso_data['Mini'] - iso_data['Mass'] < 0.1)
    plt.plot(iso_data['G_BPmag'][idx_iso] - iso_data['G_RPmag'][idx_iso], iso_data['Gmag'][idx_iso],
             lw=1.5, label='[M/H] = {:4.1f}'.format(meh), ls=l_meh[i_meh], color='black', alpha=0.7)

idx_sin_p = np.in1d(cannon_data['sobject_id'], single_plus)
idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
idx_bin_p = np.in1d(cannon_data['sobject_id'], binaries_plus)
idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)
idx_norm = (idx_sin_p + idx_bin + idx_bin_p + idx_mul) == 0

plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_norm], Gmag_twin_abs[idx_norm], s=8, lw=0, label='', c='black', alpha=1)
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=14, c='C2', label='    3')
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin_p], Gmag_twin_abs[idx_bin_p], lw=0, s=14, c='C3', label=u'$\geq$2')
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=14, c='C0', label='    2')
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_sin_p], Gmag_twin_abs[idx_sin_p], lw=0, s=14, c='C1', label=u'$\geq$1')

# plt.plot(bprp_MS_ebv, gmag_MS_ebv, c='C0', label='')

if 'ebv' in suffix:
    plt.ylabel(u'M$_G$ - A$_G$')
    plt.xlabel(u'G$_{BP}$ - G$_{RP}$ - A')
    plt.xlim(0.68, 1.)
    plt.xticks([0.7, 0.8, 0.9, 1.0], ['0.7', '0.8', '0.9', '1.0'])
    plt.ylim(2.2, 5.5)
else:
    plt.ylabel(u'M$_G$')
    plt.xlabel(u'G$_{BP}$ - G$_{RP}$')
    plt.xlim(0.69, 1.1)
    plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1], ['0.7', '0.8', '0.9', '1.0', '1.1'])
    plt.ylim(2.2, 5.5)

# plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete) - bin_bp_rb, ls='--', c="C2", alpha=0.8,label='Binaries limit')
plt.gca().invert_yaxis()
plt.legend()
# plt.legend(loc=2)
plt.grid(ls='--', alpha=0.3)
plt.tight_layout()
# plt.show()
plt.savefig('mag_hr_gaia_bin-multi_iso_4.5Gyr_res'.format(mh_ini)+suffix+'.png', dpi=300)
plt.close()

# raise SystemExit

# idx_mark = cannon_data['flag_cannon'] != 0
idx_mark = cannon_data['flag_cannon'] < 0
for c_col in ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon', 'Vsini_cannon']:
    plt.figure(figsize=(7, 5.5))
    plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0,
                vmin=np.nanpercentile(cannon_data[c_col][~idx_mark], 1), vmax=np.nanpercentile(cannon_data[c_col][~idx_mark], 99),
                c=cannon_data[c_col], label='')
    plt.colorbar()
    plt.scatter((cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=5, label='Flaged', c='C3')
    if 'ebv' in suffix:
        plt.ylabel(u'M$_G$ - A$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$ - A')
        plt.xlim(0.68, 1.)
        plt.xticks([0.7, 0.8, 0.9, 1.0], ['0.7', '0.8', '0.9', '1.0'])
        plt.ylim(2.0, 6.5)
    else:
        plt.ylabel(u'M$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$')
        plt.xlim(0.69, 1.1)
        plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1], ['0.7', '0.8', '0.9', '1.0', '1.1'])
        plt.ylim(2.2, 5.5)

    plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete) - bin_bp_rb, ls='--', c='black', alpha=0.3, label='Binaries limit')
    plt.gca().invert_yaxis()
    plt.title(cannon_version + ': ' + c_col)
    plt.tight_layout()
    plt.legend()
    plt.savefig('mag_hr_gaia_' + c_col +suffix+'.png', dpi=250)
    # plt.savefig('mag_hr_gaia_bin-multi.png', dpi=250)
    plt.close()

    plt.figure(figsize=(7, 3.5))
    c_range = np.nanpercentile(cannon_data[c_col][~idx_mark], [0.5, 99.5])
    plt.hist(cannon_data[c_col][~idx_above], range=c_range, bins=22, alpha=0.3, label='Single candidates', color='C2')
    plt.hist(cannon_data[c_col][~idx_above], range=c_range, bins=22, alpha=1, label='', color='C2', histtype='step')
    plt.hist(cannon_data[c_col][idx_above], range=c_range, bins=22, alpha=0.3, label='Multiple candidates', color='C1')
    plt.hist(cannon_data[c_col][idx_above], range=c_range, bins=22, alpha=1, label='', color='C1', histtype='step')
    plt.axvline(np.nanmedian(cannon_data[c_col][~idx_above]), label='', color='C2', ls='--')
    plt.axvline(np.nanmedian(cannon_data[c_col][idx_above]), label='', color='C1', ls='--')
    plt.xlabel(u'$v \sin i$ [km s$^{-1}$]')
    plt.ylabel('Number of candidates')
    # plt.title(c_col)
    plt.grid(ls='--', alpha=0.2, color='black')
    plt.tight_layout()
    plt.legend()
    plt.savefig('hist_' + c_col + suffix + '.png', dpi=250)
    plt.close()

# raise SystemExit

# for g_p in ['astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']:
#     x_range = (np.nanpercentile(cannon_data[g_p], 1), np.nanpercentile(cannon_data[g_p], 99))
#     plt.hist(cannon_data[~idx_mag_multiple][g_p], range=x_range, bins=250, color='blue', alpha=0.5)
#     plt.hist(cannon_data[idx_mag_multiple][g_p], range=x_range, bins=250, color='red', alpha=0.5)
#     plt.gca().set_yscale('log')
#     plt.savefig('Gp_'+g_p+'.png', dpi=250)
#     plt.close()

# do a MS fit on selected objects
idx_mag_multiple = Gmag_twin_abs <= (gg_fit(Gmag_bp_rp) - bin_bp_rb)

print 'N candidates:', np.sum(idx_mag_multiple)
print ','.join([str(sid) for sid in cannon_data['sobject_id'][idx_mag_multiple]])

plt.figure(figsize=(7, 6))
# fig, ax = plt.subplots(1, 2, sharex=True)
ax = list([])
ax.append(plt.axes([0.12, 0.28, 0.87, 0.70]))
ax.append(plt.axes([0.12, 0.09, 0.87, 0.19]))

# all photometric data together plot - and their distribution
to_abs_mag = (-2.5*np.log10(((cannon_data['r_est'])/10.)**2)).reshape(-1, 1)
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
    pc.set_alpha(0.5)

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
    pc.set_alpha(0.5)

stars = ax[0].scatter(phot_data_multi_pos, sun_mag_list, lw=0, s=100, c='black', marker='*', label='Sun', alpha=0.8)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red')
ax[0].legend([mpatches.Patch(color='C2', alpha=0.5), mpatches.Patch(color='C1', alpha=0.5), stars],
           ['Single candidates', 'Multiple candidates', 'Sun'], loc=2)
ax[0].set(ylim=(0, 6.5), xlim=(-1, 44), ylabel=u'Absolute magnitude\n ',
          xticks=np.arange(len(c_plot))*3 + 0.5, xticklabels=['' for i in range(len(c_plot))])
if 'ebv' in suffix:
    ax[0].set(ylabel=u'Absolute magnitude - extinction\n ')
ax[0].invert_yaxis()
ax[0].grid(ls='--', color='black', alpha=0.2)

print 'Med diffs:', median_2-median_1
ax[1].scatter(np.arange(len(c_plot))*3 + 0.5, median_2-median_1, lw=0, s=10, c='black')
ax[1].set(ylim=(0.45, 0.65), xlim=(-1, 44), ylabel=u'Difference', xlabel=u'Photometric filter',
          xticks=np.arange(len(c_plot)) * 3 + 0.5,
          xticklabels=['B', 'V', "g'", "r'", "i'", 'G', r'G$_{BP}$', r'G$_{RP}$', 'J', 'H', 'K', 'W1', 'W2', 'W3', 'W4'],
          yticks=[0.5, 0.6])
ax[1].grid(ls='--', color='black', alpha=0.2)

plt.savefig('multi_mag_plot'+suffix+'.png', dpi=250)
plt.close()
