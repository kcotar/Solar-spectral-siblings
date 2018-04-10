import numpy as np
from os import system
from astropy.table import Table, vstack
from solar_siblings_functions import get_used_elements
import pandas as pd

root_dir = '/home/klemen/data4_mount/'
sub_dirs = ['Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_1',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_2',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_3']
dir_out = '/home/klemen/Solar-spectral-siblings/Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_comb/'

root_dir = '/home/klemen/data4_mount/'
sub_dirs = ['Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_1',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_2',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_3',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_4',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_5',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_6']
dir_out = '/home/klemen/Solar-spectral-siblings/Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_interp_comb/'


root_dir = '/home/klemen/Solar-spectral-siblings/'
sub_dirs = ['Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_1',
            'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_2']
dir_out = '/home/klemen/Solar-spectral-siblings/Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb_2/'
system('mkdir '+dir_out)


sim_data = list([])
gp_data = list([])
gp_sim_data = list([])

for sub_d in sub_dirs:
    print sub_d
    try:
        gp_data.append(Table.read(root_dir + sub_d + '/GP_fit_res.txt', format='ascii.csv'))
        gp_sim_data.append(Table.read(root_dir + sub_d + '/solar_similarity_b1234_gp.csv', format='ascii.csv'))
        sim_data.append(Table.read(root_dir + sub_d + '/solar_similarity_b1234.fits'))
    except:
        print 'Something went wrong'
        pass

if len(sim_data) > 0:
    vstack(sim_data).write(dir_out+'solar_similarity_b1234.fits', overwrite=True)
vstack(gp_sim_data).write(dir_out+'solar_similarity_b1234_gp.csv', format='ascii.csv', overwrite=True)
vstack(gp_data).write(dir_out+'GP_fit_res.txt', format='ascii.csv', overwrite=True)

# stack abundance files if they exist
for element in get_used_elements():
    print element
    gp_sim_data = list([])
    for sub_d in sub_dirs:
        print sub_d
        try:
            gp_sim_data.append(Table.read(root_dir + sub_d + '/solar_similarity_b1234_gp_'+element+'.csv', format='ascii.csv'))
        except:
            print 'Something went wrong'
            pass
    vstack(gp_sim_data).write(dir_out + 'solar_similarity_b1234_gp_'+element+'.csv', format='ascii.csv', overwrite=True)

# stack GP created median spectra
for i_b in [1, 2, 3, 4]:
    print i_b
    wvl_data = list([])
    flx_data = list([])
    for sub_d in sub_dirs:
        print sub_d
        try:
            flx_data.append(pd.read_csv(root_dir + sub_d + '/gp_median_flx_b'+str(i_b)+'.csv', header=None, sep=',').values)
        except:
            print 'Something went wrong'
            pass
    flx_data = np.vstack(flx_data)
    np.savetxt(dir_out + 'gp_median_flx_b'+str(i_b)+'.csv', flx_data, delimiter=',', fmt='%0.3f')
