import matplotlib.pyplot as plt
import pandas as pd
import os
from glob import glob
from astropy.table import Table

os.chdir('/home/klemen/data4_mount/Solar_data_dr52/')
galah_linelist = Table.read('/home/klemen/data4_mount/GALAH_Cannon_linelist_newer.csv')

for spectrum_file in glob('twilight_spectrum_galah_ext0_date*.txt'):
    spectrum = pd.read_csv(spectrum_file, header=None, delimiter=' ', na_values='nan').values
    plt.plot(spectrum[:,0], spectrum[:,1])
for line in galah_linelist:
    plt.axvspan(line['line_start'], line['line_end'], lw=0, color='black', alpha=0.2)
plt.ylim((0,1.3))
plt.show()
plt.close()
