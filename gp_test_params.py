from solar_siblings_functions import *


suffix_solar_ref = '_ext0_dateall_offset'
solar_input_dir = galah_data_input+'Solar_data_dr53/'
solar_wvl, solar_flx = get_solar_data(solar_input_dir, suffix_solar_ref, every_nth=8)

n_px = 750
n_sample = 50

solar_wvl = solar_wvl[:n_px]
solar_flx = solar_flx[:n_px]
plt.plot(solar_wvl, solar_flx, c='blue', alpha=1, lw=1.)
solar_flx_class = mean_flux_class(solar_flx)

solar_flx_orig = np.array(solar_flx)
solar_flx = solar_flx - (1. - solar_flx)*0.3
solar_flx += np.random.randn(n_px)*0.03
# solar_flx = spectrum_offset_norm(9.9e-01, solar_flx)

gp = george.GP(get_kernel(np.array([3.05118524e-04,10.93553260e-03,1.92497689e-05,20.26058399e+00])*np.array([1.,1.,1.,1.])), mean=solar_flx_class)
# gp = george.GP(get_kernel([1, 0.0001, 1, 0.0002]))
gp.compute(solar_wvl, yerr=0.03)
gp_noise_pred = gp.sample_conditional(solar_flx, solar_wvl, size=n_sample)
# gp_noise_pred_orig = gp.sample_conditional(solar_flx_orig, solar_wvl, size=n_sample)

plt.plot(solar_wvl, solar_flx, c='red', alpha=1, lw=0.3)
plt.plot(solar_wvl, np.median(gp_noise_pred, axis=0), c='black', alpha=1, lw=0.3)
# plt.plot(solar_wvl, np.median(gp_noise_pred_orig, axis=0), c='green', alpha=1, lw=0.3)
for i_pred in range(n_sample):
    plt.plot(solar_wvl, gp_noise_pred[i_pred, :], c='black', alpha=0.3, lw=0.3)
    # plt.plot(solar_wvl, gp_noise_pred_orig[i_pred, :], c='green', alpha=0.3, lw=0.3)
plt.show()
plt.close()
