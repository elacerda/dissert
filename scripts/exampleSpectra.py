#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt
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

K = P.K

#xxx
nSpec = 5
mask = P.maskQFlag & P.maskEmLines & P.maskLambdaConstrains
mask2d = (mask[..., np.newaxis] & np.ones_like(K.f_obs, dtype = np.bool))
mf_obs__lz = np.ma.masked_array(K.f_obs, mask = ~mask2d)
mf_obs_sigma__l = mf_obs__lz.std(axis = 1)
mf_obs_norm__lz = np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d)
mf_obs_norm_sigma__l = mf_obs_norm__lz.std(axis = 1)

scinot = 1.e-16

f, axArr = plt.subplots(nSpec)
f.set_size_inches((10, 2 * nSpec))

for i, zone in enumerate(np.asarray([0, 1./5, 2./5, 3./5, 4./5]) * K.N_zone):
    z_i = np.int(zone)

    ax = axArr[i]
    ax.plot(K.l_obs, mf_obs__lz[:, z_i] / scinot, 'k-')
    ax.set_ylabel(r'Zona %i' % z_i)
    ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
    ax.set_xlim(args.lc)
    ax.grid()

    for xmin in ax.xaxis.get_minorticklocs():
        ax.axvline(x = xmin, ls = ':', c = 'grey')

    plt.setp(ax.get_xticklabels(), visible = False)

f.suptitle(u'Exemplos de espectro da galáxia %s' % K.galaxyName)
ax = f.add_axes( [0., 0., 1, 1] )
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.text(.02, 0.5, r'$F_{obs}\ [10^{-16} erg/s/cm^2/\AA]$', rotation='vertical', horizontalalignment='center', verticalalignment='center')
plt.setp(axArr[nSpec - 1].get_xticklabels(), visible = True, rotation = 45)
f.savefig('%s/%s-exampleSpectra.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))

#########################################################################################

f, axArr = plt.subplots(1,2)
f.set_size_inches((10, 3))
ax1 = axArr[0]
ax2 = axArr[1]

plus1sigma = mf_obs__lz.mean(axis = 1) + 1. * mf_obs_sigma__l
plus2sigma = mf_obs__lz.mean(axis = 1) + 2. * mf_obs_sigma__l

plus2sigma_norm = mf_obs_norm__lz.mean(axis = 1) + 2. * mf_obs_norm_sigma__l
minus2sigma_norm = mf_obs_norm__lz.mean(axis = 1) - 2. * mf_obs_norm_sigma__l

p5, p50, p95 = np.percentile(mf_obs__lz, [5, 50, 95], axis = 1)
p5[~mask] = np.nan
p50[~mask] = np.nan
p95[~mask] = np.nan
p5n, p50n, p95n = np.percentile(mf_obs_norm__lz, [5, 50, 95], axis = 1)
p5n[~mask] = np.nan
p50n[~mask] = np.nan
p95n[~mask] = np.nan
mfo_max__l = mf_obs__lz.max(axis = 1)
mfo_min__l = mf_obs__lz.min(axis = 1)

mfo_max_norm__l = mf_obs_norm__lz.max(axis = 1)
mfo_min_norm__l = mf_obs_norm__lz.min(axis = 1)

#ax1.plot(K.l_obs, p5 / scinot, 'k-', lw = 0.5)
#ax1.plot(K.l_obs, p50 / scinot, 'k-', lw = 0.5)
ax1.plot(K.l_obs, p95 / scinot, 'k-', lw = 0.5)
ax1.fill_between(K.l_obs, mfo_max__l / scinot, mfo_min__l / scinot, edgecolor='gray', facecolor='lightgray')
ax1.set_ylabel(r'$F_{obs}\ [10^{-16} erg/s/cm^2/\AA]$')
ax1.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax1.grid()

#ax2.plot(K.l_obs, p5n, 'k-', lw = 0.5)
#ax2.plot(K.l_obs, p50n, 'k-', lw = 0.5)
ax2.plot(K.l_obs, p95n, 'k-', lw = 0.5)
ax2.fill_between(K.l_obs, mfo_max_norm__l, mfo_min_norm__l, edgecolor='gray', facecolor='lightgray')
ax2.set_ylabel(r'$f_{obs}$')
ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax2.grid()

#f.subplots_adjust(left=0.07, bottom=0.1, top=0.95, wspace=0.2, hspace=0)
f.savefig('%s/%s-exampleSpectraFill.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
