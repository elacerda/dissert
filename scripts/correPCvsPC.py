#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import PCAlifa as PCA
import argparse as ap

###############################################################################
###############################################################################

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
    parser.add_argument('--outputdir', '-o',
                        help = 'Image output directory',
                        metavar = 'DIR',
                        type = str,
                        default = 'png')

    return parser.parse_args()

def corr_PCvsPC(tomo__zk, nPCs, title, fnamepref):
    prop = {
        'arr'   : [ P.K.at_flux__z, np.log10(P.K.aZ_flux__z / 0.019), P.K.A_V, P.K.v_0, P.K.v_d ],
        'label' : [ r'$\log\ t\ [yr]$', r'$\log\ Z\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$' ],
        'fname' : [ 'at_flux', 'aZ_flux', 'AV', 'v0', 'vd' ],
        'cmap'  : [ 'hot_r', 'hot_r', 'hot_r', 'spectral', 'spectral' ],
    }

    for p_i, p in enumerate(prop['arr']):
        z = p
        nRows = nPCs - 1
        nCols = nPCs - 1
        f, axArr = plt.subplots(nRows, nCols)

        for ax in f.axes:
            ax.set_axis_off()

        fig_width_pt = 1080.
        inches_per_pt = 1.0 / 72.27
        golden_mean = (5 ** 0.5 - 1.0) / 2.0
        fig_width = fig_width_pt * inches_per_pt
        fig_height = fig_width * golden_mean

        f.set_size_inches(fig_width, fig_height)
        f.set_dpi(100)

        for i in range(nRows):
            y = tomo__zk[:, i]

            for j in range(i, nCols):
                ax = axArr[i, j]
                ax.set_axis_on()

                plt.setp(ax.get_xticklabels(), visible = False)
                plt.setp(ax.get_yticklabels(), visible = False)

                x = tomo__zk[:, j + 1]

                ax.set_xlim(x.min(), x.max())
                ax.set_ylim(y.min(), y.max())

                if (prop['fname'])[p_i] == 'vd' :
                    prc = np.percentile(p, 98.)
                    im = ax.scatter(x, y, c = z, edgecolor = 'None', marker = 'o', s = 0.1, cmap = prop['cmap'][p_i], vmax = prc)
                else:
                    im = ax.scatter(x, y, c = z, edgecolor = 'None', marker = 'o', s = 0.1, cmap = prop['cmap'][p_i])

                plt.setp(ax.get_yticklabels(), visible = False)

            axArr[i, i].set_ylabel('PC%d' % (i + 1))
            plt.setp(axArr[i, i].get_xticklabels(), visible = True, rotation = 45)
            plt.setp(axArr[i, i].get_yticklabels(), visible = True, rotation = 45)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.0)
        f.subplots_adjust(right = 0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = f.colorbar(im, cax = cbar_ax)
        cb.set_label(prop['label'][p_i])

        for i in range(nCols):
            axArr[0, i].set_title('PC%d' % (i + 2))

        plt.suptitle(r'%s' % title)
        f.savefig('%scorre_PCxPC_%s.%s' % (fnamepref, prop['fname'][p_i], args.outputimgsuffix))

###############################################################################
###############################################################################

args = parser_args()
print('Output directory: %s' % args.outputdir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

nPCs = 6

P.PCA_obs()
P.tomograms_obs()
corr_PCvsPC(P.tomo_obs__zk, nPCs, r'$F_{obs}$.', '%s/%s-f_obs-' % (args.outputdir, P.K.califaID))

P.PCA_obs_norm()
P.tomograms_obs_norm()
corr_PCvsPC(P.tomo_obs_norm__zk, nPCs, r'$f_{obs}$', '%s/%s-f_obs_norm-' % (args.outputdir, P.K.califaID))

P.PCA_syn_norm()
P.tomograms_syn_norm()
corr_PCvsPC(P.tomo_syn_norm__zk, nPCs, r'$f_{syn}$', '%s/%s-f_syn_norm-' % (args.outputdir, P.K.califaID))
