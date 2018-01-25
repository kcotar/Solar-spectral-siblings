import numpy as np
from astropy.table import Table
from astropy.io import fits

galah_data_input = '/home/klemen/data4_mount/'

zoo_data = Table.read(galah_data_input+'zooniverse_subjects.csv')
# galah_data = Table.read(galah_data_input+'sobject_iraf_52_reduced_20171111.fits')
# galah_data = Table.read('/media/storage/HERMES_REDUCED/dr5.3/sobject_iraf_53.fits')
galah_data = Table.read('/media/storage/HERMES_REDUCED/dr5.2/sobject_iraf_52.fits')
cannon_data = Table.read(galah_data_input+'sobject_iraf_iDR2_171103_cannon.fits')
galah_data = galah_data[galah_data['flag_guess']==0]
cannon_data = cannon_data[cannon_data['flag_guess']==0]

set_pid = 16051
zoo_data = zoo_data[zoo_data['subject_set_id'] == set_pid]

sobject_ids = zoo_data['metadata']
field_ids = [np.int64(s.split(':')[1].split('_')[0][2:-3]+'000') for s in sobject_ids]
sobject_ids = [np.int64(s.split(':')[1].split('_')[0][2:]) for s in sobject_ids]

print sobject_ids
print field_ids
print

# open fits
fits_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
for s_id in sobject_ids:
    str_s_id = str(s_id)
    fits_path = fits_dir + str_s_id[0:6] + '/standard/com/' + str_s_id + '1.fits'
    # print fits_path
    fits_data = fits.open(fits_path, memmap=True)
    print str_s_id, ',', fits_data[0].header.get('STD_NAME')


print
print galah_data[np.in1d(galah_data['sobject_id'], sobject_ids)]['sobject_id','ra', 'dec', 'galah_id','teff_guess','feh_guess', 'logg_guess', 'flag_guess']
print cannon_data[np.in1d(cannon_data['sobject_id'], sobject_ids)]['sobject_id','galah_id', 'ra','dec','Teff_cannon','Feh_cannon','Logg_cannon','flag_cannon']


# for f_id in sobject_ids:
#     print galah_data[(galah_data['sobject_id']  ==f_id) < 400]['sobject_id','ra','dec','teff_guess','feh_guess', 'logg_guess']
#     print cannon_data[(cannon_data['sobject_id'] == f_id) < 400]['sobject_id', 'Teff_cannon','Feh_cannon','Logg_cannon','flag_cannon']