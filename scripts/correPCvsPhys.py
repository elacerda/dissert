#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import PCAlifa as PCA
import argparse as ap

def parser_args():
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--fitsfile', '-f',
                        help = 'The file must be named KXXXX*.fits',
                        metavar = 'PyCASSO FITS FILE',
                        type = str,
                        default = None)
    parser.add_argument('--lc', '-l',
                        help = 'Lambda constrains',
                        metavar = 'LAMBDA',
                        type = int,
                        nargs = 2,
                        default = [3800, 6850])
    parser.add_argument('--outputimgsuffix', '-o',
                        help = 'Suffix of image file. Sometimes denote the image type. (Ex.: image.png)',
                        type = str,
                        default = 'png')
    parser.add_argument('--outputdir', '-d',
                        help = 'Image output directory',
                        metavar = 'DIR',
                        type = str,
                        default = 'png')

    return parser.parse_args()

###############################################################################
###############################################################################

args = parser_args()

print('Output directory: %s' % args.outputdir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_obs()
P.tomograms_obs()
P.PCA_obs_norm()
P.tomograms_obs_norm()
P.PCA_syn_norm()
P.tomograms_syn_norm()

K = P.K

#xxx

def plot_corr_axes(x, y, ax):
    rhoPearson, pvalPearson = st.pearsonr(x, y)
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = 's: %.2f' % rhoSpearman

    ax.scatter(x, y, marker = 'o', s = 0.1)

    textbox = dict(boxstyle='round', facecolor='wheat', alpha=0.)
    
    ax.text(0.76, 0.15, txt,
            fontsize = 8, 
            transform = ax.transAxes,
            verticalalignment = 'top',
            bbox = textbox)
    
    plt.setp(ax.get_yticklabels(), visible = False)

def plot_evec_ax(l, evec, ax, *kwargs):
    ax.plot(l, evec, *kwargs)
    ax.axhline(y = 0, c = 'r')
    plt.setp(ax.get_yticklabels(), visible = False)

prop = {
    'arr'   : [ K.at_flux__z, np.log10(K.aZ_flux__z / 0.019), K.A_V, K.v_0, K.v_d ],
    'label' : [ r'$\langle \log\ t \rangle_L\ [yr]$', r'$\log\ \langle Z \rangle_L\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$', r'eigenvector' ],
}

nPCs = 6

nCols = len(prop['label'])

f, axArr = plt.subplots(nPCs, nCols)
f.set_size_inches(18, 10)

for i in range(nPCs):
    iPC = i
    nPC = i + 1
    axArr[i, 0].set_ylabel('PC%d' % nPC)
    
    if iPC == 2:
        plot_evec_ax(P.l_obs, P.eigVec_obs_norm__lk[:, iPC], axArr[i, nCols - 1], 'k-')
        y = P.tomo_obs_norm__zk[:, iPC]
    else:
        plot_evec_ax(P.l_obs, -1. * P.eigVec_obs_norm__lk[:, iPC], axArr[i, nCols - 1], 'k-')
        y = -1 * P.tomo_obs_norm__zk[:, iPC]

    for j in range(nCols)[:-1]:
        ax = axArr[i, j]
        ax.set_ylim(y.min(), y.max())

        x = prop['arr'][j]

        plot_corr_axes(x, y, ax)

        if prop['label'][j] == r'$\sigma_\star\ [km/s]$':
            ax.set_xlim(0, 230)

f.subplots_adjust(hspace = 0.05)
f.subplots_adjust(wspace = 0.1)

plt.setp([a.get_xticklabels() for a in f.axes], rotation = 90)
plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

for i in range(nCols):
    axArr[0, i].set_title(prop['label'][i])

f.suptitle(u'Correlações $f_{obs}$', fontsize=16)
f.savefig('%s/%s-correl-f_obs_norm-PCvsPhys.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))

f, axArr = plt.subplots(nPCs, nCols)
f.set_size_inches(18, 10)

for i in range(nPCs):
    iPC = i
    nPC = i + 1
    axArr[i, 0].set_ylabel('PC%d' % nPC)

    if i > 1:
        y = P.tomo_syn_norm__zk[:, iPC]
        plot_evec_ax(P.l_obs, P.eigVec_syn_norm__lk[:, iPC], axArr[i, nCols - 1], 'k-')
    else:
        y = -1. * P.tomo_syn_norm__zk[:, iPC]
        plot_evec_ax(P.l_obs, -1. * P.eigVec_syn_norm__lk[:, iPC], axArr[i, nCols - 1], 'k-')
    

    for j in range(nCols)[:-1]:
        ax = axArr[i, j]
        ax.set_ylim(y.min(), y.max())

        x = prop['arr'][j]

        plot_corr_axes(x, y, ax)

        if prop['label'][j] == r'$\sigma_\star\ [km/s]$':
            ax.set_xlim(0, 230)


f.subplots_adjust(hspace = 0.05)
f.subplots_adjust(wspace = 0.1)

plt.setp([a.get_xticklabels() for a in f.axes], rotation = 90)
plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

for i in range(nCols):
    axArr[0, i].set_title(prop['label'][i])

f.suptitle(u'Correlações $f_{syn}$', fontsize=16)
f.savefig('%s/%s-correl-f_syn_norm-PCvsPhys.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))

f, axArr = plt.subplots(nPCs, nCols)
f.set_size_inches(18, 10)

for i in range(nPCs):
    iPC = i
    nPC = i + 1
    axArr[i, 0].set_ylabel('PC%d' % nPC)

    if i == 2:
        y = P.tomo_obs__zk[:, iPC]
        plot_evec_ax(P.l_obs, P.eigVec_obs__lk[:, iPC], axArr[i, nCols - 1], 'k-')
    else:
        y = -1. * P.tomo_obs__zk[:, iPC]
        plot_evec_ax(P.l_obs, -1. * P.eigVec_obs__lk[:, iPC], axArr[i, nCols - 1], 'k-')

    for j in range(nCols)[:-1]:
        ax = axArr[i, j]
        ax.set_ylim(y.min(), y.max())

        x = prop['arr'][j]

        plot_corr_axes(x, y, ax)

        if prop['label'][j] == r'$\sigma_\star\ [km/s]$':
            ax.set_xlim(0, 230)


f.subplots_adjust(hspace = 0.05)
f.subplots_adjust(wspace = 0.1)

plt.setp([a.get_xticklabels() for a in f.axes], rotation = 90)
plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

for i in range(nCols):
    axArr[0, i].set_title(prop['label'][i])

f.suptitle(u'Correlações $F_{obs}$', fontsize=16)
f.savefig('%s/%s-correl-f_obs-PCvsPhys.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))

