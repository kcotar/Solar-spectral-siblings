import numpy as np
from astropy.table import Table, join

galah_data_input = '/home/klemen/data4_mount/'

# data read
gaia_data = Table.read('Gaia_benchmark_data.csv')
cannon_data = Table.read(galah_data_input+'sobject_iraf_iDR2_180325_cannon.fits')
galah_gaia_match = Table.read('Gaia_benchmark_stars_GALAH.csv')['sobject_id', 'std_name']

# apply filters and flags
# cannon_data = cannon_data[cannon_data['flag_cannon'] == 0]
# gaia_data = gaia_data[np.abs(gaia_data['Teff_bench'] - 5777) < 500]
# gaia_data = gaia_data[np.abs(gaia_data['Logg_bench'] - 4.44) < 0.5]
# gaia_data = gaia_data[np.abs(gaia_data['Feh_bench'] - 0.02) < 0.5]

# join sets
galah_gaia_match = join(galah_gaia_match, gaia_data, keys='std_name')
galah_gaia_match = join(galah_gaia_match, cannon_data, keys='sobject_id')
print 'matched stars:', len(galah_gaia_match)

# output result
galah_gaia_match['Teff_diff'] = galah_gaia_match['Teff_cannon'] - galah_gaia_match['Teff_bench']
galah_gaia_match['Logg_diff'] = galah_gaia_match['Logg_cannon'] - galah_gaia_match['Logg_bench']
galah_gaia_match['Feh_diff'] = galah_gaia_match['Fe_H_cannon'] - galah_gaia_match['Feh_bench']

print galah_gaia_match['sobject_id', 'std_name', 'Teff_bench', 'e_Teff_bench', 'Teff_cannon', 'e_Teff_cannon', 'flag_cannon']
print galah_gaia_match['sobject_id', 'std_name', 'Logg_bench', 'e_Logg_bench', 'Logg_cannon', 'e_Logg_cannon', 'flag_cannon']
print galah_gaia_match['sobject_id', 'std_name', 'Feh_bench', 'e_Feh_bench', 'Fe_H_cannon', 'e_Fe_H_cannon', 'flag_cannon']
print galah_gaia_match['sobject_id', 'std_name', 'Teff_diff', 'Logg_diff', 'Feh_diff', 'flag_cannon', 'red_flag']