import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, unique
from os import chdir
import astropy.coordinates as coord
import astropy.units as un

galah_data_input = '/home/klemen/data4_mount/'
isochrones_dir = '/home/klemen/data4_mount/isochrones/padova_Gaia_DR2_Solar/'

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

chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb')
gp_res = Table.read('solar_similarity_b1234_gp.csv')
print len(gp_res)

# predetermined objects
solar_like_sobjects = gp_res['sobject_id']

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]
cannon_data = join(cannon_data, gp_res, join_type='left', keys='sobject_id')
print len(cannon_data)
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
parallax_sim = np.linspace(1, 14, 100)
Gmag_sim = G_mean + 2.5*np.log10(((1e3/parallax_sim)/10.)**2)

Gmag_twin_abs = cannon_data['phot_g_mean_mag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
Vmag_twin_abs = cannon_data['Vmag'] - 2.5*np.log10(((1e3/cannon_data['parallax'])/10.)**2)
cannon_data['g_mean_mag_abs'] = Gmag_twin_abs

# idx_mark = np.in1d(cannon_data['sobject_id'], [150208003201286,150427004801275,160130006301234,160327004601337,160524006601258,160916001801263,161118004701221,161119002801292,170117003101044,170205005401120,170515003101036,170516002101273,171001001601082,171207003601278,180129003101184])
# idx_mark = np.in1d(cannon_data['sobject_id'], li_high_sid)
# idx_mark = Gmag_twin_abs < 4.2
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
plt.errorbar(cannon_data['phot_g_mean_mag'], cannon_data['parallax'], yerr=cannon_data['parallax_error'],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0, label='Gaia observations')
plt.errorbar(cannon_data['phot_g_mean_mag'][idx_mark], cannon_data['parallax'][idx_mark], yerr=cannon_data['parallax_error'][idx_mark],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='red', markeredgewidth=0, label='')
plt.xlabel('Gaia G mean magnitude')

plt.ylabel('Gaia parallax')
plt.xlim((9, 14))
plt.ylim((0.5, 14))
plt.grid(color='black', ls='--', alpha=0.2)
plt.legend()
plt.tight_layout()
plt.savefig('mag_parallax_gaia.png', dpi=250)
plt.close()

# load isochrones that will be added tot he plot
iso_data = Table.read(isochrones_dir+'isochrones_all.fits')
iso_data = iso_data[iso_data['Age'] == 4500000000.]

for c_col in ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']:
    plt.scatter(cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'], Gmag_twin_abs, s=4, lw=0, c=cannon_data[c_col], label='')
    plt.colorbar()
    plt.scatter((cannon_data['phot_bp_mean_mag']-cannon_data['phot_rp_mean_mag'])[idx_mark], Gmag_twin_abs[idx_mark], lw=0, s=7, label='', c='red')
    plt.ylabel('Gmag absolute')
    plt.xlabel(r'G$_{bp}$ - G$_{rp}$')
    plt.ylim(2, 6)
    plt.xlim(0.7, 1.35)

    # add isochrone
    for i_feh in np.unique(iso_data['MHini']):
        idx_iso = np.logical_and(iso_data['MHini'] == i_feh, iso_data['Mini']-iso_data['Mass'] < 0.1)
        plt.plot(iso_data['G_BPmag'][idx_iso]-iso_data['G_RPmag'][idx_iso], iso_data['Gmag'][idx_iso],
                 lw=0.75, label='[M/H] = {:.1f}'.format(i_feh))

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.legend()
    plt.savefig('mag_hr_gaia_'+c_col+'.png', dpi=300)
    plt.close()

# idx = np.logical_and(Gmag_twin_abs < 4.2, cannon_data['phot_g_mean_mag'] < 12)
# print cannon_data[Gmag_twin_abs > 4.2]['sobject_id', 'source_id', 'ra_2', 'dec_2', 'phot_bp_rp_excess_factor', 'phot_variable_flag', 'a_g_val', 'g_mean_mag_abs']
# print cannon_data[idx]['sobject_id', 'source_id', 'ra_2', 'dec_2', 'phot_bp_rp_excess_factor', 'phot_variable_flag', 'a_g_val', 'g_mean_mag_abs', 'parallax']

idx_mag_multiple = Gmag_twin_abs < 4.2
# print cannon_data[~idx_mag_multiple]['sobject_id', 'source_id', 'l', 'b', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']
# print cannon_data[idx_mag_multiple]['sobject_id', 'source_id', 'l', 'b', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig']

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
idx_mag_multiple = Gmag_twin_abs < 4.2

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
