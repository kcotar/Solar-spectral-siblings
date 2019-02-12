import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting


def fit_MS_lin_line(Gmag_twin_abs, Gmag_bp_rp, path=False, d_above=0.2):
    bp_rp_s = 0.01
    bp_rp_w = 0.02
    bp_rp_r = np.nanpercentile(Gmag_bp_rp, [20, 80])
    bp_rp_eval = np.arange(bp_rp_r[0], bp_rp_r[1], bp_rp_s)
    bp_rp_eval_complete = np.arange(0., 2., 0.1)
    idx_fit_use = Gmag_bp_rp > 0
    for i_n in range(5):
        # print 'Median on:', np.sum(idx_fit_use)
        gmag_med = list([])
        for bp_rp in bp_rp_eval:
            gmag_med.append(np.nanmedian(Gmag_twin_abs[np.logical_and(np.logical_and(Gmag_bp_rp > bp_rp-bp_rp_w/2.,
                                                                                     Gmag_bp_rp <= bp_rp+bp_rp_w/2.),
                                                                      idx_fit_use)]))
        gmag_med = np.array(gmag_med)

        f_init = models.Linear1D(slope=2., intercept=2.)
        fitter = fitting.LevMarLSQFitter()
        gg_fit = fitter(f_init, bp_rp_eval, gmag_med)
        idx_fit_use = np.logical_and(idx_fit_use, gg_fit(Gmag_bp_rp)-d_above < Gmag_twin_abs)
        # print gg_fit

        if path is not None:
            plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete), c="C"+str(i_n), alpha=0.5)
            plt.plot(bp_rp_eval_complete, gg_fit(bp_rp_eval_complete) - d_above, ls='--', c="C"+str(i_n), alpha=0.5)

    if path is not None:
        plt.scatter(Gmag_bp_rp, Gmag_twin_abs, lw=0, s=5)
        plt.scatter(bp_rp_eval, gmag_med, lw=0, s=5)
        plt.xlim(0.60, 1.15)
        plt.ylim(6.5, 2.0)
        plt.savefig(path, dpi=250)
        plt.close()

    return gg_fit
