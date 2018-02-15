import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, imp
from glob import glob
from astropy.table import Table
from astropy.modeling import models, fitting
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit

imp.load_source('helper_functions', '../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize, spectra_resample


spectra_dir = '/media/storage/HERMES_REDUCED/dr5.3/'
galah_linelist = Table.read('/home/klemen/data4_mount/GALAH_Cannon_linelist_newer.csv')
telluric_spectrum_O2 = np.loadtxt('/home/klemen/data4_mount/telluric_O2_conv.dat')
telluric_spectrum_H2O = np.loadtxt('/home/klemen/data4_mount/telluric_H2O_conv.dat')
solar_ref_conv = pd.read_csv('/home/klemen/data4_mount/Solar_data_dr53/solar_spectra_conv.txt', header=None, delimiter=' ', na_values='nan').values
run_cob_id = Table.read('observations.csv')

# renormalize data to the same continuum level
sl = 2.
sh = 3.
st = 15
ord = 1
# telluric_spectrum_O2[:, 1] = spectra_normalize(telluric_spectrum_O2[:, 0]-np.mean(telluric_spectrum_O2[:, 0]), telluric_spectrum_O2[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
# telluric_spectrum_H2O[:, 1] = spectra_normalize(telluric_spectrum_H2O[:, 0]-np.mean(telluric_spectrum_H2O[:, 0]), telluric_spectrum_H2O[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
# solar_ref_conv[:, 1] = spectra_normalize(solar_ref_conv[:, 0]-np.mean(solar_ref_conv[:, 0]), solar_ref_conv[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
# solar_ref_conv_res = spectra_resample(solar_ref_conv[:, 1], solar_ref_conv[:, 0], telluric_spectrum_H2O[:, 0], k=1)

os.chdir('/home/klemen/data4_mount/Solar_data_dr53/')
for spectrum_file in glob('twilight_spectrum_galah_ext0_date*.txt'):
    # determine date of flat
    cob_id = spectrum_file.split('_')[-1].split('.')[0][-10:]
    print cob_id
    if cob_id == 'dateall':
        continue
    run_ids = np.unique(run_cob_id[run_cob_id['cob_id'] == np.int64(cob_id)]['run_id'])
    print run_ids
    # determine needed diagnostic files are named like 09sep30018.txt
    diag_dir = spectra_dir + cob_id[:6] + '/diagnostics/'
    for ib in [4]:
        tel_spectrum = []
        for run_id in run_ids:
            tel_spectrum = []
            for run_tel_file in glob(diag_dir+'*'+str(ib)+str(run_id)[-4:]+'.txt'):
                data = np.loadtxt(run_tel_file)
                if ib == 4:  # remove portion of the nir band with strong tellurics
                    idx_use = data[:, 0] > 7670
                    data = data[np.where(idx_use)[0], :]
                    tel_spect_norm = data[:, 3]
                else:
                    # normalize spectrum
                    tel_spect_norm = spectra_normalize(data[:, 0]-np.mean(data[:, 0]), data[:, 3], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
                # tel_spectrum.append(tel_spect_norm)
                tel_spectrum = tel_spect_norm  # np.mean(np.vstack(tel_spectrum), axis=0)
                tel_wvl = data[:, 0] / 1.000275  # use refractive index of air to shift tellurics into observed frame

                h2o_spectrum = spectra_resample(telluric_spectrum_H2O[:, 1], telluric_spectrum_H2O[:, 0], tel_wvl, k=1)
                o2_spectrum = spectra_resample(telluric_spectrum_O2[:, 1], telluric_spectrum_O2[:, 0], tel_wvl, k=1)

                # rescale spectrum model to the telluric spectrum as determined by molecfit
                h2o_w = 1.
                o2_w = 1.
                tel_model = 1. - ((1. - h2o_spectrum) * h2o_w + (1. - o2_spectrum) * o2_w)

                params = Parameters()
                params.add('a1', value=h2o_w, min=0., max=15.)
                params.add('a2', value=o2_w, min=0., max=15.)

                def tel_scale_func(params, data, h2o, o2, eval=False):
                    model = 1. - ((1. - h2o) * params['a1'] + (1. - o2) * params['a2'])
                    if eval:
                        return model
                    else:
                        return model - data

                minner = Minimizer(tel_scale_func, params, fcn_args=(tel_spectrum, h2o_spectrum, o2_spectrum))
                result = minner.minimize()
                report_fit(result)

                # write error report

                plt.plot(tel_wvl, tel_scale_func(result.params, 1, h2o_spectrum, o2_spectrum, eval=True), c='green', alpha=0.5, lw=1., label='Telluric model fit')
                plt.plot(tel_wvl, tel_model, c='black', alpha=0.5, lw=1., label='Telluric model')
                plt.plot(tel_wvl, tel_spectrum, c='red', alpha=0.5, lw=1., label='Molecfit normalized')
                plt.legend()
                plt.show()
                plt.close()

        #
        #         plt.plot(data[:, 0], data[:, 1], c='black', alpha=0.3, lw=.3)
        #         plt.plot(data[:, 0], data[:, 3], c='red', alpha=0.3, lw=.3)
        # plt.ylim((0.5, 1.1))
        # plt.show()
        # # plt.savefig('tel_'+cob_id+'_'+str(ib)+'.png', dpi=450)
        # plt.close()

    print ''
    continue
    spectrum = pd.read_csv(spectrum_file, header=None, delimiter=' ', na_values='nan').values
    spectrum[:, 1] = spectra_normalize(spectrum[:, 0] - np.mean(spectrum[:, 0]), spectrum[:, 1], steps=st, sigma_low=sl, sigma_high=sh, order=ord, func='poly')
    plt.plot(spectrum[:, 0], spectrum[:, 1], c='red', lw=1, alpha=0.5)
    # add telluric spectrum
    plt.plot(solar_ref_conv[:, 0], solar_ref_conv[:, 1], c='black', lw=1.5)
    telluric_contr = (1. - telluric_spectrum_H2O[:, 1])*4. + (1. - telluric_spectrum_O2[:, 1])*3.
    plt.plot(telluric_spectrum_H2O[:, 0], solar_ref_conv_res - telluric_contr, lw=1.5, c='blue')
    for line in galah_linelist:
        plt.axvspan(line['line_start'], line['line_end'], lw=0, color='black', alpha=0.2)
    plt.ylim((0, 1.3))
    plt.show()
    plt.close()

