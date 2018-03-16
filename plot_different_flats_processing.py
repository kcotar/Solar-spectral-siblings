import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('/home/klemen/data4_mount/Solar_data_dr53/')
solar_ref_conv = np.loadtxt('solar_spectra_conv.txt')

plt.plot(solar_ref_conv[:, 0], solar_ref_conv[:, 1], label='Ref', lw=2, c='black')
for p in [5, 7, 9, 11]:
    s_d = np.loadtxt('twilight_spectrum_galah_ext0_dateall_p'+str(p)+'.txt')
    plt.plot(s_d[:, 0], s_d[:, 1], label=str(p))
plt.legend()
plt.show()
plt.close()
