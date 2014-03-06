#!/usr/bin/python
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
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

args = parser_args()

print('Output directory: %s' % args.outputdir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_obs_norm()
P.tomograms_obs_norm()

K = P.K

#xxx

nPCs = 5
hfig = nPCs * 5

f = plt.figure(figsize = (15, hfig))
gs = gridspec.GridSpec(nPCs, 2, width_ratios = [5, 8])

col = 0

for i in range(nPCs):
    tn = i + 1
    gsi = 2 * i

    if i == 2:
        tomo_pc = P.tomo_obs_norm__kyx[i, :, :]
        pc = P.eigVec_obs_norm__lk[:, i]
    else:
        tomo_pc = -1. * P.tomo_obs_norm__kyx[i, :, :]
        pc = -1. * P.eigVec_obs_norm__lk[:, i]

    eval = P.eigVal_obs_norm__k[i]
    evals = P.eigVal_obs_norm__k
    f_obs_mean = P.ms_obs_norm__l
    l = P.l_obs

    ax1 = plt.subplot(gs[gsi])
    ax2 = plt.subplot(gs[gsi + 1])
    ax1.set_title(r'tomograma $%02i$' % tn)

    im = ax1.imshow(tomo_pc, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
    f.colorbar(ax = ax1, mappable = im, use_gridspec = True)
    eval_perc = 100. * eval / evals.sum()

    if i == 0:
        ax2.set_title(r'PCA com $f_{obs}$. - var: $%.2f\ \%%$' % eval_perc)
    else:
        ax2.set_title(r'var: $%.2f\ \%%$' % eval_perc)

    ax2.plot(l, pc, 'k')
    ax2.set_ylabel(r'$PC %02i$' % tn)
#    plt.setp(ax2.get_xticklabels(), rotation = 45)
    ax2.set_ylabel(r'$\lambda$ [\AA]')
    ax2.grid()
    ax3 = ax2.twinx()
    ax3.plot(l, f_obs_mean, color = '0.65')
    ax3.set_ylabel(r'Espectro medio')
    ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

    for xmin in ax2.xaxis.get_minorticklocs():
        ax2.axvline(x = xmin, ls = ':', c = 'grey')

f.tight_layout()
f.savefig('%s/%s-tomo-obs-norm.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
