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
cannon_data = cannon_data[cannon_data['sobject_id']>140301000000000]
ra_dec_all = coord.ICRS(ra=cannon_data['ra']*un.deg, dec=cannon_data['dec']*un.deg)

chdir('Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb')
gp_res = Table.read('solar_similarity_b1234_gp.csv')

# predetermined objects
solar_like_sobjects = gp_res['sobject_id']

# cannon data subsets
cannon_data = cannon_data[np.in1d(cannon_data['sobject_id'], solar_like_sobjects)]
cannon_data = join(cannon_data, gp_res, join_type='left')


def _prepare_hist_data(d, bins, range, norm=True):
    heights, edges = np.histogram(d, bins=bins, range=range)
    width = np.abs(edges[0] - edges[1])
    if norm:
        heights = 1.*heights / np.nanmax(heights)
    return edges[:-1], heights, width


# make histogram plots for parameters
plot_col = ['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']
x_label = ['Teff [K]', 'Logg [dex]', '[Fe/H] [dex]']
solar_val = [5771, 4.44, 0.02]
idx_valid_param = cannon_data['flag_cannon'] == 0
for i_p in range(len(plot_col)):
    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plot_vals_median = np.median(plot_vals)
    plt.hist(plot_vals, bins=75, color='black', alpha=0.33)
    plt.axvline(x=solar_val[i_p], color='black', ls='--', lw=2.5)
    plt.axvline(x=plot_vals_median, color='red', ls='--', lw=2.5)
    plt.title('All: {:.0f}    Unflagged: {:.0f}    Median: {:.2f}    Difference: {:.2f}'.format(len(solar_like_sobjects), len(plot_vals), plot_vals_median, plot_vals_median-solar_val[i_p]))
    plt.xlabel(x_label[i_p])
    plt.ylabel('Number of objects in bin')
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1)
    # plt.show()
    plt.savefig('solar_twins_like_'+plot_col[i_p]+'.png', dpi=300)
    plt.close()

    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plt.scatter(plot_vals, cannon_data['canberra'][idx_valid_param], s=4, lw=0, c='black')
    plt.tight_layout()
    # plt.show()
    plt.savefig('solar_twins_like_' + plot_col[i_p] + '_canberra.png', dpi=300)
    plt.close()

plot_col = [col for col in cannon_data.colnames if '_abund_cannon' in col and 'e_' not in col and 'flag_' not in col and len(col.split('_')) < 4]
plot_col = np.sort(plot_col)
element = [col.split('_')[0] for col in plot_col]

# one plot for
n_x = 5
n_y = 6
fig, ax = plt.subplots(n_y, n_x, figsize=(8, 8), sharex='all', sharey='all')
fig.subplots_adjust(hspace=0, wspace=0, left=0.025, right=0.975, top=0.975, bottom=0.075)

for i_p in range(len(plot_col)):
    x_p = i_p % n_x
    y_p = int(i_p / n_x)
    idx_valid_param = cannon_data['flag_'+plot_col[i_p]] == 0
    plot_vals = cannon_data[plot_col[i_p]][idx_valid_param]
    plot_vals_median = np.median(plot_vals)
    plot_vals_std = 2. * np.std(plot_vals)
    h_edg, h_hei, h_wid = _prepare_hist_data(plot_vals, 100, (-1., 1.))
    # plots
    ax[y_p, x_p].bar(h_edg, h_hei, width=h_wid, color='black', alpha=0.4)
    ax[y_p, x_p].axvline(x=plot_vals_median, color='red', ls='--', lw=1.5, alpha=0.75)
    ax[y_p, x_p].axvline(x=plot_vals_median - plot_vals_std, color='red', ls='--', lw=1.5, alpha=0.25)
    ax[y_p, x_p].axvline(x=plot_vals_median + plot_vals_std, color='red', ls='--', lw=1.5, alpha=0.25)
    # make it nicer
    ax[y_p, x_p].set_yticks([0.25, 0.5, 0.75])
    ax[y_p, x_p].set_yticklabels(['', '', ''])
    ax[y_p, x_p].set_xticks([-0.75, 0., 0.75])
    ax[y_p, x_p].grid(ls='--', alpha=0.15, color='black')
    # write out labels to the plot
    ax[y_p, x_p].text(-0.9, 0.85, element[i_p])  # element name
    n_valid = np.sum(idx_valid_param)
    ax[y_p, x_p].text(-0.9, 0.475, '{:0.0f}'.format(n_valid))  # number of valid cannon observations
    if n_valid > 0:
        ax[y_p, x_p].text(-0.9, 0.1, '{:0.2f}'.format(plot_vals_median - 0.))  # difference
    # add xlabel
    if y_p == n_y-1 and x_p == 2:
        ax[y_p, x_p].set_xlabel('[X/Fe] abundance distribution')  # element name
plt.savefig('abund_all_selection.png', dpi=300)
plt.close()

# raise SystemExit

# make histogram plots for abundances
solar_val = np.zeros(len(plot_col))
for i_p in range(len(plot_col)):
    idx_valid_param = cannon_data['flag_'+plot_col[i_p]] == 0
    print plot_col[i_p]+': '
    if 'Li' in plot_col[i_p]:
        idx_Li = cannon_data[plot_col[i_p]] > 0.1
        sob_Li = cannon_data[np.logical_and(idx_Li, idx_valid_param)]['sobject_id']
        # print ','.join([str(s) for s in sob_Li])

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

# make positional plots, l/b and galaxy position based on gaia dr2
gp_final_selected = [150208003201286,150427004801275,160130006301234,160327004601337,160524006601258,160916001801263,161118004701221,161119002801292,170117003101044,170205005401120,170515003101036,170516002101273,171001001601082,171207003601278,180129003101184]
idx_hig = np.in1d(cannon_data['sobject_id'], gp_final_selected)
gaia_data = Table.read(galah_data_input+'CANNON_180325_Gaia_DR2_xmatch.csv')['angDist','sobject_id','source_id','parallax','parallax_error','phot_g_mean_mag', 'pmra', 'pmdec','radial_velocity']
gaia_data = gaia_data[gaia_data['angDist'] < 0.5]
gaia_data = unique(gaia_data, keys='sobject_id', keep='first')
gaia_data['parallax'][gaia_data['parallax_error']>gaia_data['parallax']] = np.nan
cannon_data = join(cannon_data, gaia_data, keys='sobject_id', join_type='left')
# cannon_data = cannon_data.filled(-1)
cannon_data['parsec'] = 1e3/cannon_data['parallax']
print cannon_data['sobject_id','phot_g_mean_mag','parsec','parallax','parallax_error']
cannon_data['sobject_id','phot_g_mean_mag','parsec','parallax','parallax_error'].write('distance_data.fits', overwrite=True)
cannon_data['sobject_id','phot_g_mean_mag','parsec','parallax','parallax_error'][idx_hig].write('distance_data_sel.fits', overwrite=True)

plt.errorbar(cannon_data['phot_g_mean_mag'], cannon_data['parallax'], yerr=cannon_data['parallax_error'],
             fmt='.', ms=5, elinewidth=1, alpha=0.8, color='black', markeredgewidth=0)
plt.xlabel('Gaia G mean magnitude')
plt.ylabel('Gaia parallax')
plt.tight_layout()
plt.savefig('mag_parallax.png', dpi=250)
plt.close()

ra_dec = coord.ICRS(ra=cannon_data['ra']*un.deg, dec=cannon_data['dec']*un.deg,
                    distance=cannon_data['parsec']*un.pc)

# galactic cartesian
gal_coord = ra_dec.transform_to(coord.Galactocentric).cartesian

print 'Max dist:', np.nanmax(cannon_data['parsec'])
print 'Max dist sub:', np.nanmax(cannon_data['parsec'][idx_hig])

print 'Z height:', np.nanmax(gal_coord.z.value), np.nanmin(gal_coord.z.value)
print 'Z height:', np.nanmax(gal_coord.z.value[idx_hig]), np.nanmin(gal_coord.z.value[idx_hig])

plt.hist(gal_coord.z.value, range=(-700, 700), bins=100)
plt.savefig('z_hist.png', dpi=250)
plt.close()

# l/ coordinates and plot them
print 'l b plot'
l_b = ra_dec.transform_to(coord.Galactic)
l_b_all = ra_dec_all.transform_to(coord.Galactic)
plt.figure()
plt.subplot(111, projection="mollweide")
plt.scatter(np.deg2rad(l_b_all.l.value)-np.pi, np.deg2rad(l_b_all.b.value), lw=0, s=1, c='#D0D0D0')
plt.scatter(np.deg2rad(l_b.l.value)-np.pi, np.deg2rad(l_b.b.value), lw=0, s=5, c='black')
plt.scatter(np.deg2rad(l_b.l.value)[idx_hig]-np.pi, np.deg2rad(l_b.b.value)[idx_hig], lw=0, s=5, c='red')
plt.grid(alpha=0.75, ls='--')
plt.xlabel('Galactic longitude')
plt.ylabel('Galactic latitude')
plt.tight_layout()
plt.savefig('l_b.png', dpi=350)
plt.close()





