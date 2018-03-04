from os import chdir
from scipy.io import readsav

file = 'data.sav'

data = readsav(file)
s_data = data['data'][0]

chdir('3D_spectra')
txt_rm = open('readme.txt', 'w')
for i_f in range(len(s_data[0])):
    filen_name_out = s_data[0][i_f]
    str_rm = str(filen_name_out)+','+str(s_data[1][i_f])+','+str(s_data[2][i_f])+','+str(s_data[3][i_f])+','+str(s_data[4][i_f])
    print str_rm
    txt_rm.write(str_rm+'\n')
    s_wvl = s_data[5][i_f, :]
    s_flx = s_data[6][i_f, :]
    s_con = s_data[7][i_f, :]
    # output
    txt_out = open(filen_name_out+'.csv', 'w')
    for i_l in range(len(s_wvl)):
        txt_out.write('{:f},{:f},{:f}\n'.format(s_wvl[i_l], s_flx[i_l], s_con[i_l]))
    txt_out.close()
txt_rm.close()
