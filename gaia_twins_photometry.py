import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, unique
from os import chdir
import astropy.coordinates as coord
import astropy.units as un

galah_data_input = '/home/klemen/data4_mount/'

# data-table settings
data_date = '20180327'
cannon_param_file = 'sobject_iraf_iDR2_180325_cannon.fits'

cannon_data = Table.read(galah_data_input+cannon_param_file)
apass_data = Table.read(galah_data_input+'photometry/apass_dr53_20180327.csv')['sobject_id','Vmag','e_Vmag','Bmag','e_Bmag']
apass_data = unique(apass_data, keys='sobject_id', keep='first')

chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb')
gp_res = Table.read('solar_similarity_b1234_gp.csv')
gp_res_li = Table.read('solar_similarity_b1234_gp_Li.csv')
li_high_sid = gp_res_li[gp_res_li['canberra'] > 0.01]['sobject_id']

# predetermined objects
solar_like_sobjects = gp_res['sobject_id']

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]
cannon_data = join(cannon_data, gp_res, join_type='left', keys='sobject_id')
print len(cannon_data)
cannon_data = join(cannon_data, apass_data, join_type='left', keys='sobject_id')
print len(cannon_data)


# absolute values in Vega magnitudes
B_s = 5.44
V_s = 4.81
R_s = 4.43
I_s = 4.10

# Gaia relations for Gmag
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
gaia_data['parallax'][gaia_data['parallax_error'] > gaia_data['parallax']] = np.nan
cannon_data = join(cannon_data, gaia_data, keys='sobject_id', join_type='left')
# cannon_data = cannon_data.filled(-1)
cannon_data['parsec'] = 1e3/cannon_data['parallax']
# print cannon_data['sobject_id', 'phot_g_mean_mag', 'parsec', 'parallax', 'parallax_error']

# simulate parallax and Gmag as would be observe by the Gaia spacecraft
parallax_sim = np.linspace(1, 14, 100)
Gmag_sim = G_mean + 2.5*np.log10(((1e3/parallax_sim)/10.)**2)

Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
Vmag_twin_abs = cannon_data['Vmag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
cannon_data['g_mean_mag_abs'] = Gmag_twin_abs

# idx_mark = np.in1d(cannon_data['sobject_id'], [150208003201286,150427004801275,160130006301234,160327004601337,160524006601258,160916001801263,161118004701221,161119002801292,170117003101044,170205005401120,170515003101036,170516002101273,171001001601082,171207003601278,180129003101184])
# idx_mark = np.in1d(cannon_data['sobject_id'], li_high_sid)
idx_mark = Gmag_twin_abs < 4.2

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


plt.plot(Gmag_sim, parallax_sim, lw=1, color='C2', ls='--', label='Relation for Sun')
plt.errorbar(cannon_data['phot_g_mean_mag'], cannon_data['parallax'], yerr=cannon_data['parallax_error'],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0, label='Gaia observations')
# plt.scatter(cannon_data['phot_g_mean_mag'][idx_mark], cannon_data['parallax'][idx_mark], lw=0, s=10, label='', c='red')
plt.xlabel('Gaia G mean magnitude')

# plt.errorbar(cannon_data['Vmag'], cannon_data['parallax'], yerr=cannon_data['parallax_error'],
#              fmt='.', ms=5, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0, label='APASS observations')
# plt.scatter(cannon_data['Vmag'][idx_mark], cannon_data['parallax'][idx_mark], lw=0, s=10, label='', c='red')
# plt.xlabel('APASS V magnitude')

plt.ylabel('Gaia parallax')
plt.xlim((9, 14))
plt.ylim((0.5, 14))
plt.grid(color='black', ls='--', alpha=0.2)
plt.legend()
plt.tight_layout()
plt.savefig('mag_parallax_gaia.png', dpi=250)
plt.close()

for c_col in ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']:
    # plt.scatter((cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=7, label='', c='red')
    # plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, c=cannon_data[c_col])
    # plt.ylabel('Gmag absolute')
    # plt.xlabel(r'G$_{bp}$ - G$_{rp}$')

    plt.scatter((cannon_data['Bmag']-cannon_data['Vmag'])[idx_mark], Vmag_twin_abs[idx_mark], lw=0, s=7, label='', c='red')
    plt.scatter(cannon_data['Bmag']-cannon_data['Vmag'], Vmag_twin_abs, s=4, lw=0, c=cannon_data[c_col])
    plt.ylabel('APASS V')
    plt.xlabel(r'APASS B-V')

    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.tight_layout()
    plt.savefig('mag_hr_APASS_'+c_col+'.png', dpi=300)
    plt.close()

idx = np.logical_and(Gmag_twin_abs < 4.2, cannon_data['phot_g_mean_mag'] < 12)
# print cannon_data[Gmag_twin_abs > 4.2]['sobject_id', 'source_id', 'ra_2', 'dec_2', 'phot_bp_rp_excess_factor', 'phot_variable_flag', 'a_g_val', 'g_mean_mag_abs']
print cannon_data[idx]['sobject_id', 'source_id', 'ra_2', 'dec_2', 'phot_bp_rp_excess_factor', 'phot_variable_flag', 'a_g_val', 'g_mean_mag_abs', 'parallax']

print cannon_data[Gmag_twin_abs > 4.2]['sobject_id', 'source_id', 'l', 'b', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']
print cannon_data[Gmag_twin_abs < 4.2]['sobject_id', 'source_id', 'l', 'b', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']

for g_p in ['astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']:
    range = (np.nanpercentile(cannon_data[g_p], 1), np.nanpercentile(cannon_data[g_p], 99))
    plt.hist(cannon_data[Gmag_twin_abs > 4.2][g_p], range=range, bins=250, color='blue', alpha=0.5)
    plt.hist(cannon_data[Gmag_twin_abs < 4.2][g_p], range=range, bins=250, color='red', alpha=0.5)
    plt.gca().set_yscale('log')
    plt.savefig('Gp_'+g_p+'.png', dpi=250)
    plt.close()

print ','.join([str(s) for s in cannon_data[Gmag_twin_abs < 4.2]['sobject_id']])
