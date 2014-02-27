#!/usr/bin/python
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

output_fmt = 'pdf'

fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = [3800, 6850])
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_syn_norm()
P.tomograms_syn_norm()

K = P.K

nPCs = 5
hfig = nPCs * 5
#hfig = 15

f = plt.figure(figsize = (15, hfig))
gs = gridspec.GridSpec(nPCs, 2, width_ratios = [5, 8])
#f = plt.figure(figsize = (30, hfig))
#gs = gridspec.GridSpec(3, 4, width_ratios = [5, 8, 5, 8])

for i in range(nPCs):
    tn = i + 1
    gsi = 2 * i

    if i > 1:
        tomo_pc = P.tomo_syn_norm__kyx[i, :, :]
        pc = P.eigVec_syn_norm__lk[:, i]
    else:
        tomo_pc = -1. * P.tomo_syn_norm__kyx[i, :, :]
        pc = -1. * P.eigVec_syn_norm__lk[:, i]

    eval = P.eigVal_syn_norm__k[i]
    evals = P.eigVal_syn_norm__k
    f_syn_mean = P.ms_syn_norm__l
    l = P.l_syn

    ax1 = plt.subplot(gs[gsi])
    ax2 = plt.subplot(gs[gsi + 1])
    ax1.set_title(r'tomograma $%02i$' % tn)

    im = ax1.imshow(tomo_pc, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap='hot_r')
    f.colorbar(ax = ax1, mappable = im, use_gridspec = True)
    eval_perc = 100. * eval / evals.sum()

    if i == 0:
        ax2.set_title(r'PCA com $f_{syn}$. - var: $%.2f\ \%%$' % eval_perc)
    else:
        ax2.set_title(r'var: $%.2f\ \%%$' % eval_perc)

    ax2.plot(l, pc, 'k')
    ax2.set_ylabel(r'$PC %02i$' % tn)
    plt.setp(ax2.get_xticklabels(), rotation = 45)
    ax2.grid()
    ax3 = ax2.twinx()
    ax3.plot(l, f_syn_mean, color = '0.65')
    ax3.set_ylabel(r'Espectro medio')
    ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

    for xmin in ax2.xaxis.get_minorticklocs():
        ax2.axvline(x = xmin, ls = ':', c = 'grey')

f.tight_layout()
f.savefig('%s/%s-tomo-syn-norm.%s' % (output_dir, K.califaID, output_fmt))
