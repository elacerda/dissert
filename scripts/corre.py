#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
from matplotlib import pyplot as plt
from scipy import stats as st
import PCAlifa as PCA

#plt.rcParams['text.usetex'] = True

fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = [3800, 6850])
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.runPCA()

K = P.K

prop = {
    'arr'   : [ K.at_flux__z, np.log10(K.aZ_flux__z / 0.019), K.A_V, K.v_0, K.v_d ],
    'label' : [ r'$\log\ t\ [yr]$', r'$\log\ Z\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$', r'eigenvector' ],
#    'xmin'  : [ 0,0,0,0,0,0 ],
#    'xmax'  : [ 0,0,0,0,200,0 ],
}

nPC = 6

nCols = len(prop['label'])

f, axArr = plt.subplots(nPC, nCols)
f.set_size_inches(19.2, 10.8)

for i in range(nPC):
    #iPC = i + 4
    iPC = i
    axArr[i, 0].set_ylabel('PC%d' % iPC)
    
    #y = P.tomo_obs__zk[:, iPC]
    y = P.tomo_obs_norm__zk[:, iPC]
    #y = P.tomo_syn_norm__zk[:, iPC]

    for j in range(nCols)[:-1]:
        ax = axArr[i, j]
        ax.set_ylim(y.min(), y.max())

        x = prop['arr'][j]

        rhoPearson, pvalPearson = st.pearsonr(x, y)
        rhoSpearman, pvalSpearman = st.spearmanr(x, y)
        #txt = 'p: %.2f\ns: %.2f' % (rhoPearson, rhoSpearman)
        txt = 's: %.2f' % rhoSpearman
        ax.scatter(x, y, marker = 'o', s = 0.1)
        #, label = 's: %.2f' % rhoSpearman)
        #leg = ax.legend(loc='best', fancybox=True)


        if prop['label'][j] == r'$\sigma_\star\ [km/s]$':
            ax.set_xlim(0, 230)

        textbox = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        
        ax.text(0.72, 0.15, txt,
                fontsize = 10, 
                transform = ax.transAxes,
                verticalalignment = 'top',
#                weight = 'bold',
                bbox = textbox)
#
#                horizontalalignment = 'right',
#                verticalalignment = 'center',
#                multialignment = 'right',
        
        plt.setp(ax.get_yticklabels(), visible = False)

    #axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_obs__lk[:, iPC], 'k-')
    axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_obs_norm__lk[:, iPC], 'k-')
    #axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_syn_norm__lk[:, iPC], 'k-')
    plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

f.subplots_adjust(hspace = 0.0)
f.subplots_adjust(wspace = 0.1)

plt.setp([a.get_xticklabels() for a in f.axes], rotation = 90)
plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

for i in range(nCols):
    axArr[0, i].set_title(prop['label'][i])

#f.suptitle(u'Correlações $F_{obs}$', fontsize=16)
#f.savefig('%s/%s-correl-f_obs-PCvsPhys.pdf' % (output_dir, K.califaID))

f.suptitle(u'Correlações $F_{obs}$ norm.', fontsize=16)
f.savefig('%s/%s-correl-f_obs_norm-PCvsPhys.pdf' % (output_dir, K.califaID))

#f.suptitle(u'Correlações $F_{syn}$ norm.', fontsize=16)
#f.savefig('%s/%s-correl-f_syn_norm-PCvsPhys.pdf' % (output_dir, K.califaID))
