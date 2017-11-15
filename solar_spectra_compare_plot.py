from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import varconvolve as varcon
from scipy import mgrid

def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)

data_input_dir = '/home/klemen/GALAH_data/Solar_data/'

solar_ref = pd.read_csv(data_input_dir + 'solar_spectra.txt', header=None, delimiter=' ', na_values='nan').values
solar_g1 = pd.read_csv(data_input_dir + 'b1_solar_galah.txt', header=None, delimiter=' ', na_values='nan').values
solar_g2 = pd.read_csv(data_input_dir + 'b2_solar_galah.txt', header=None, delimiter=' ', na_values='nan').values
solar_g3 = pd.read_csv(data_input_dir + 'b3_solar_galah.txt', header=None, delimiter=' ', na_values='nan').values
solar_g4 = pd.read_csv(data_input_dir + 'b4_solar_galah.txt', header=None, delimiter=' ', na_values='nan').values

x = solar_ref[:, 0]
y = solar_ref[:, 1]
v = np.ones_like(x)
y1 = varcon.varconvolve(x, y, kernel, v*0.08)
y2 = varcon.varconvolve(x, y, kernel, v*0.09)
y3 = varcon.varconvolve(x, y, kernel, v*0.1)
y4 = varcon.varconvolve(x, y, kernel, v*0.11)

plt.plot(x, y, c='black', lw=2)
plt.plot(x, y1, c='red', lw=2)
plt.plot(x, y2, c='green', lw=2)
plt.plot(x, y3, c='magenta', lw=2)
plt.plot(x, y4, c='orange', lw=2)

plt.plot(solar_g1[:, 0], solar_g1[:, 1], c='blue', lw=2)
plt.plot(solar_g2[:, 0], solar_g2[:, 1], c='blue', lw=2)
plt.plot(solar_g3[:, 0], solar_g3[:, 1], c='blue', lw=2)
plt.plot(solar_g4[:, 0], solar_g4[:, 1], c='blue', lw=2)
plt.show()
plt.close()
