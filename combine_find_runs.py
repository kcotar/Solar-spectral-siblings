import numpy as np
from astropy.table import Table, vstack

root_dir = '/home/klemen/data4_mount/'
sub_dirs = ['Distances_Step2_p0_SNRsamples750_ext0_oklinesonly_origsamp_1',
            'Distances_Step2_p0_SNRsamples750_ext0_oklinesonly_origsamp_2',
            'Distances_Step2_p0_SNRsamples750_ext0_oklinesonly_origsamp_3']

dir_out = '/home/klemen/Solar-spectral-siblings/Distances_Step2_p0_SNRsamples750_ext0_oklinesonly_origsamp_comb/'

sim_data = list([])
gp_data = list([])
for sub_d in sub_dirs:
    print sub_d
    try:
        gp_data.append(Table.read(root_dir + sub_d + '/GP_fit_res.txt', format='ascii.csv'))
        sim_data.append(Table.read(root_dir + sub_d + '/solar_similarity_b1234.fits'))
    except:
        pass

vstack(sim_data).write(dir_out+'solar_similarity_b1234.fits', overwrite=True)
vstack(gp_data).write(dir_out+'GP_fit_res.txt', format='ascii.csv', overwrite=True)