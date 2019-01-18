import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, unique
from os import chdir
from scipy.interpolate import spline
import astropy.coordinates as coord
import astropy.units as un
from dustmaps.bayestar import BayestarWebQuery


gmag_MS_ebv = [6.46440852913791, 6.406135638759127, 6.326063057758769, 6.270662313171322, 6.186901128845153, 6.112655265836556, 6.047964331920713, 5.975597998849081, 5.864775024621654, 5.752935773085529, 5.70615132715316, 5.607932103569, 5.510605841548255, 5.394819317583499, 5.320540582224229, 5.206141954573217, 5.105337368096131, 4.997662315562364, 4.882464233710049, 4.756352703758858, 4.644881458729178, 4.527907502273222, 4.408222001484199, 4.2871798967317405, 4.117467298488372, 3.994238242422548, 3.861777335534732, 3.72487273109238, 3.633423030730339, 3.477971871745545, 3.3525566275598284, 3.1855673220231195, 2.9850365948134816, 2.830982397100216]
bprp_MS_ebv = [1.2908581206431133, 1.2588210724583755, 1.2268114976583142, 1.2073005303122892, 1.1797431892275174, 1.1542746627289153, 1.1258597376937827, 1.0968862630807177, 1.0720300797647617, 1.0410309764107133, 1.02608802792935, 0.9972604104823493, 0.9779415187202556, 0.9411407394522913, 0.942366454124274, 0.9106812550604175, 0.8997248316011284, 0.877800406566954, 0.8593017040995061, 0.8325740629627321, 0.8233155569426165, 0.7927167614889337, 0.7802655762838091, 0.7615338470099067, 0.7389512071183755, 0.7321674777076357, 0.7169088241276089, 0.6948854046524628, 0.6878793860732788, 0.6748835060196896, 0.63182029448188, 0.6369998924908571, 0.6304756782387, 0.6046636421838443]

plt.rcParams['font.size'] = 15

galah_data_input = '/shared/ebla/cotar/'
isochrones_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_Solar/'
evoltracks_dir = galah_data_input + 'isochrones/padova_Gaia_DR2_evolutionary_track/'

# data-table settings
data_date = '20180327'
cannon_param_file = 'galah_cannon_gaia_photometry_20180327_ebv-corr.fits'
cannon_data = Table.read(galah_data_input+cannon_param_file)

suffix = '_ebv'
# suffix = ''

chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb')
gp_res = Table.read('solar_similarity_b1234_gp.csv')
print len(gp_res)
# predetermined objects
multi_mag_thr = 4.2
solar_like_sobjects = gp_res['sobject_id']

idx_run_params = 9
params_str = ['5100_4.54_0.00', '5200_4.53_0.00', '5300_4.51_0.00', '5400_4.48_0.00',
              '5500_4.45_0.00', '5600_4.41_0.00', '5700_4.36_0.00', '5800_4.31_0.00',
              '5900_4.25_0.00', '6000_4.18_0.00'][idx_run_params]
multi_mag_thr = [5.7, 5.5, 5.3, 5.2,
                 4.9, 4.6, 4.3, 3.9,
                 3.6, 3.4][idx_run_params]

chdir(galah_data_input + 'Distances_Step1_p0_SNRsamples0_ext4_oklinesonly_G20180327_C180325_refpar_'+params_str)
sel_txt = open('final_selection_0.10.txt', 'r')
solar_like_sobjects = sel_txt.read()
sel_txt.close()
solar_like_sobjects = [np.int64(sid) for sid in solar_like_sobjects.split(',')]
suffix += '_' + params_str.split('_')[0]
print suffix

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]

# print len(cannon_data)

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
print 'Total:', len(solar_like_sobjects), 'above:', np.sum(Gmag_twin_abs < multi_mag_thr)
print cannon_data['sobject_id'][np.logical_and(Gmag_twin_abs < multi_mag_thr, cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'] < 0.6)]

# idx_mark = cannon_data['parallax_error'] > cannon_data['parallax']*0.2
idx_mark = cannon_data['ruwe'] > 1.4
# idx_mark = ~np.isfinite(cannon_data['parallax'])

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
plt.xlim((9, 14))
plt.ylim((0.5, 13))
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

# 2 and >2
binaries = [140412001701123,140608001401220,140708005301189,140805002101116,140806003501357,140807004501202,140813003801242,141102002401103,141231003501047,150107002201216,150108001001207,150211002201182,150211005201034,150409004601049,150429001601375,150603002801033,150603004301022,150703001601061,150706005401044,150718004401157,150827004001317,150828004701328,150828005201087,150830005101263,151110004201125,151229005001270,151229005501152,160124002601146,160130005801230,160325003701176,160328004201105,160330102801363,160331002201006,160331002201156,160331004801074,160331005301007,160401003401328,160402004601164,160424003101103,160426004001205,160519002601136,160519003601178,160519006001177,160523005801165,160527002101137,160529005901105,160530002201108,160530003901205,160530003901360,160530005501086,160602001601163,160812003101256,160813004601070,160815003601073,160816004701025,160919003001218,160919003501225,161009004301038,161011003401158,161013001601360,161013003201388,161013004901167,161107003901105,161109002601057,161109004401022,161116004401044,161118004001126,161210004201159,161213003101091,161213003101355,161217005101052,161217005101125,161218002101211,161218002601113,170112001601131,170114003101192,170115005201161,170206003701057,170219002601044,170407005201028,170407005201070,170408003501020,170412002901220,170413003601201,170413005101386,170414004601355,170506003401087,170507006201086,170507011701283,170508004301278,170509003201249,170513004001177,170517001801075,170602002701083,170602002701315,170603008301199,170604006101059,170604006101388,170711005801066,170712004201086,170713005601065,170801004001110,170802001801229,170906004101129,170906004601078,170910001801075,170910001801104,170910002601107,170910002601108,170910006101013,170910006101377,170911003601133,170911005301354,170911005801345,171029003801283,171031003301069,171205003101140,171206005101299,171230004101299,171230004101309,171230004601085,180101002101111,180101002101231,180101003101335,180129003101072,180129003601229,180129004401240]
# 3
multiples = [141103003601188,141231003001363,150204001601334,150401003601336,150401004101146,150429001601116,150602003301082,150606002901183,150606002901197,150703005601301,160326000101229,160327005101097,160330001601096,160420005301099,160422002501087,160422003001334,160519006501256,160520002601147,160524003601122,160524004201322,160530003901281,160812002601124,160817002601186,160919005101149,161008003001047,161116003801393,161118004701207,170105003101066,170114001601106,170216002801161,170407002601147,170410003901086,170413002601061,170415001501219,170514003301376,170601003101388,170602003201116,170713001601370,170713002101086,170911002601045,170911005801328,170912002901191,171101003001066,171205004601179,171227005801095,180101001601210]
# ?
inconclusive = [140611001601132,150210004201221,150411006601018,160325003701204,160327004101204,160524006101292,160723002601396,161107001601077,161217002101146,170113001601146,170412003401337,170508001601282,170515005101004,170711002001099,171207002701147]

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

    idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
    idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
    idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)

    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[~idx_mark], Gmag_twin_abs[~idx_mark], s=7, lw=0, label='', c='C0', alpha=0.8)
    plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], s=10, lw=0, label='', c='black', marker='X', alpha=0.8)
    if 'ebv' in suffix:
        plt.ylabel(u'M$_G$ - A$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$ - A')
        plt.xlim(0.60, 1.0)
        plt.xticks([0.7, 0.8, 0.9, 1.0], ['0.7', '0.8', '0.9', '1.0'])
        plt.ylim(1.8, 5.0)
    else:
        plt.ylabel(u'M$_G$')
        plt.xlabel(u'G$_{BP}$ - G$_{RP}$')
        plt.xlim(0.69, 1.1)
        plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1], ['0.7', '0.8', '0.9', '1.0', '1.1'])
        plt.ylim(2.2, 5.5)

    plt.axhline(multi_mag_thr, ls='--', c='C2', alpha=0.8)
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

idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)
idx_norm = (idx_bin + idx_mul + idx_inc) == 0

plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_norm], Gmag_twin_abs[idx_norm], s=8, lw=0, label='', c='black', alpha=1)
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=11, c='C3', label='Triple')
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=11, c='C2', label='Binary')
plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_inc], Gmag_twin_abs[idx_inc], lw=0, s=8, c='C0', label='')

plt.plot(bprp_MS_ebv, gmag_MS_ebv, c='C0', label='')

if 'ebv' in suffix:
    plt.ylabel(u'M$_G$ - A$_G$')
    plt.xlabel(u'G$_{BP}$ - G$_{RP}$ - A')
    plt.xlim(0.60, 1.0)
    plt.xticks([0.7, 0.8, 0.9, 1.0], ['0.7', '0.8', '0.9', '1.0'])
    plt.ylim(1.8, 5.0)
else:
    plt.ylabel(u'M$_G$')
    plt.xlabel(u'G$_{BP}$ - G$_{RP}$')
    plt.xlim(0.69, 1.1)
    plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1], ['0.7', '0.8', '0.9', '1.0', '1.1'])
    plt.ylim(2.2, 5.5)

plt.axhline(multi_mag_thr, ls='--', c='C2', alpha=0.8)
plt.gca().invert_yaxis()
plt.legend()
# plt.legend(loc=2)
plt.grid(ls='--', alpha=0.3)
plt.tight_layout()
# plt.show()
plt.savefig('mag_hr_gaia_bin-multi_iso_4.5Gyr_res'.format(mh_ini)+suffix+'.png', dpi=300)
plt.close()

raise SystemExit

# for c_col in ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']:
#
#     # add isochrone
#     plt.figure(figsize=(7, 5.5))
#     for i_feh in np.unique(iso_data['MHini']):
#         idx_iso = np.logical_and(iso_data['MHini'] == i_feh, iso_data['Mini'] - iso_data['Mass'] < 0.1)
#         plt.plot(iso_data['G_BPmag'][idx_iso] - iso_data['G_RPmag'][idx_iso], iso_data['Gmag'][idx_iso],
#                  lw=0.75, label='[M/H] = {:.1f}'.format(i_feh))
#                  # lw=0.5, label='', c='black')
#
#     plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, c=cannon_data[c_col], label='')
#     plt.colorbar()
#     plt.scatter((cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=7, label='', c='C3')
#     plt.ylabel(u'Absolute ${\it Gaia}$ G magnitude')
#     plt.xlabel(u'Colour index G$_{BP}$ - G$_{RP}$')
#     plt.yticks([6, 5.5, 5, 4.5, 4, 3.5, 3, 2.5, 2], ['6', '', '5', '', '4', '', '3', '', '2'])
#     plt.xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3], ['', '0.8', '', '1.0', '', '1.2', ''])
#     plt.ylim(2, 6)
#     plt.xlim(0.7, 1.35)
#
#     # binaries = [140709003001184,150408005301259,150830004601175,160326002101077,160401004401051,160420003301342,160424002101194,160425001901273,160426004501395,160524004201306,160530003901026,160530003901054,160813002101314,160813004101136,161104002801361,161116003801308,161117004001161,170107004201309,170112003601298,170510007301226,170514003301011,170614004601055,170710002701354,150409002601317,150411004101331,150606005901339,160327006101355,160330002601095,160521004801353,160530005501077,160531005101256,160602001601307,161008002501018,161210004201315,161211003101387,161217002601138,161219005101228,170117003101044,170408004501048,170514003001180,170906003601357,170909002601291,171207003601278,171227003601367]
#     # multiples = [140608002501303,140808002701338,140808003701104,141102002701267,150408005301116,150413003601344,150413005101096,160107003101157,160125004501038,160325002701048,160401003901215,160531006101153,161009005901171,161118002601376,161212002101397,161217004101075,161217004101234,170121002801292,170507007801271,170508002601001,170508004801312,170514002401099,170514003301001,170911003101383,140608002501303,140808003701104,150413003601344,160125004501038,160401003901215,161009005901171,161118002601376,161217004101234,170507007801271,170514002401099,170911003101383]
#     # inconclusive = [161009002601018]
#     # idx_bin = np.in1d(cannon_data['sobject_id'], binaries)
#     # idx_mul = np.in1d(cannon_data['sobject_id'], multiples)
#     # idx_inc = np.in1d(cannon_data['sobject_id'], inconclusive)
#     # plt.figure(figsize=(7, 5.5))
#     # plt.scatter(cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, label='', c='black')
#     # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=8, c='C3', label='Discarded')
#     # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_bin], Gmag_twin_abs[idx_bin], lw=0, s=8, c='C0', label='Binary')
#     # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_mul], Gmag_twin_abs[idx_mul], lw=0, s=8, c='C1', label='Triple')
#     # # plt.scatter((cannon_data['phot_bp_mean_mag'] - cannon_data['phot_rp_mean_mag'])[idx_inc], Gmag_twin_abs[idx_inc], lw=0, s=7, c='C1', label='Inconclusive')
#     # plt.ylabel(u'Absolute ${\it Gaia}$ G magnitude')
#     # plt.xlabel(u'Colour index G$_{BP}$ - G$_{RP}$')
#     # plt.ylim(2, 6)
#     # plt.xlim(0.7, 1.35)
#
#     plt.axhline(multi_mag_thr, ls='--', c='black', alpha=0.3)
#     plt.gca().invert_yaxis()
#     plt.tight_layout()
#     plt.legend()
#     plt.savefig('mag_hr_gaia_' + c_col + '.png', dpi=250)
#     # plt.savefig('mag_hr_gaia_bin-multi.png', dpi=250)
#     plt.close()


# for g_p in ['astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']:
#     x_range = (np.nanpercentile(cannon_data[g_p], 1), np.nanpercentile(cannon_data[g_p], 99))
#     plt.hist(cannon_data[~idx_mag_multiple][g_p], range=x_range, bins=250, color='blue', alpha=0.5)
#     plt.hist(cannon_data[idx_mag_multiple][g_p], range=x_range, bins=250, color='red', alpha=0.5)
#     plt.gca().set_yscale('log')
#     plt.savefig('Gp_'+g_p+'.png', dpi=250)
#     plt.close()
idx_mag_multiple = Gmag_twin_abs < multi_mag_thr
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

stars = ax[0].scatter(phot_data_multi_pos, sun_mag_list, lw=0, s=100, c='black', marker='*', label='Sun', alpha=0.8)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red')
ax[0].legend([mpatches.Patch(color='C1'), mpatches.Patch(color='C2'), stars],
           ['Multiple', 'Single', 'Sun'], loc=2)
ax[0].set(ylim=(0, 6.5), xlim=(-1, 44), ylabel=u'Absolute magnitude\n ',
          xticks=np.arange(len(c_plot))*3 + 0.5, xticklabels=['' for i in range(len(c_plot))])
if 'ebv' in suffix:
    ax[0].set(ylabel=u'Absolute magnitude - extinction\n ')
ax[0].invert_yaxis()
ax[0].grid(ls='--', color='black', alpha=0.2)

ax[1].scatter(np.arange(len(c_plot))*3 + 0.5, median_2-median_1, lw=0, s=10, c='black')
ax[1].set(ylim=(0.8, 1.0), xlim=(-1, 44), ylabel=u'Difference', xlabel=u'Photometric filter',
          xticks=np.arange(len(c_plot)) * 3 + 0.5,
          xticklabels=['B', 'V', "g'", "r'", "i'", 'G', r'G$_{BP}$', r'G$_{RP}$', 'J', 'H', 'K', 'W1', 'W2', 'W3', 'W4'],
          yticks=[0.85, 0.95])
ax[1].grid(ls='--', color='black', alpha=0.2)

plt.savefig('multi_mag_plot'+suffix+'.png', dpi=250)
plt.close()
