from os import chdir
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np

bands = [1, 2, 3, 4]
b_suffx = ''.join([str(b) for b in bands])
input_dir = 'Distances_Step2_p0_SNRsamples1000_ext0_oklinesonly_origsamp_G20180327_C180325_multiabund_comb'
simulation_dir = '../Distances_SNR_models_guesslike-alllines_gauss_oklinesonly_origsamp_nonoise/'
chdir(input_dir)

sim_file_root = 'solar_similarity_b'+b_suffx+'_gp'
sim_data = Table.read(sim_file_root+'.csv')

for metric in ['canberra', 'cityblock', 'chi2']:
    ic = 0
    for offset in [0.0, 0.04, 0.08, 0.12, 0.16, 0.20]:
        simulated = Table.read(simulation_dir+'solar_similarity_narrow_b'+b_suffx+'_flux{:0.2f}.fits'.format(offset))
        plt.axhline(simulated[metric][0], label='{:0.2f}'.format(offset), alpha=0.75, c='C'+str(ic))
        ic += 1
    plt.scatter(sim_data['snr_spectrum'], sim_data[metric], lw=0, s=4, c='black', label='')
    plt.legend()
    # plt.show()
    plt.savefig(metric+'.png', dpi=300)
    plt.close()

all_best = list([])
sim_c = 'cityblock'
for sim_c in ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'minkowski', 'wchebyshev', 'sqeuclidean', 'euclidean', 'chi2', 'EW']:
    print sim_c
    perc_best = 50
    idx_best = np.full(len(sim_data), True)
    for element in ['Al','Fe','O','Si','Sc','V','Ti','Ba','C','Zn','Zr','Li','Na']:
        sim_el_data = Table.read(sim_file_root + '_' + element + '.csv')
        perc_thr = np.percentile(sim_el_data[sim_c], perc_best)
        idx_best = np.logical_and(idx_best, sim_el_data[sim_c] <= perc_thr)
        print 'Element:', element, ' left:', np.sum(idx_best)
    print sim_data[idx_best]
    s_obj = sim_data[idx_best]['sobject_id']
    all_best.append(s_obj)
    print ','.join([str(s) for s in s_obj])
    print ' - - - - - - - - - - - - - - '
    print ' - - - - - - - - - - - - - - '


all_best = np.hstack(all_best)
u_sobj, c_sobj = np.unique(all_best, return_counts=True)
for i in np.argsort(c_sobj)[::-1]:
    print u_sobj[i], c_sobj[i]
