#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import PCAlifa as PCA
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
                        default = '../figuras')

    return parser.parse_args()

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

args = parser_args()

print('Output directory: %s' % args.outputdir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_obs()
P.PCA_obs_norm()
P.PCA_syn_norm()

K = P.K

maxPCs = 15

f = plt.figure(figsize = (8, 5), dpi = 100)

PCs = np.linspace(1, maxPCs, maxPCs)

eigval_obs = 100. * P.eigVal_obs__k / P.eigVal_obs__k.sum()
eigval_obs_norm = 100. * P.eigVal_obs_norm__k / P.eigVal_obs_norm__k.sum()
eigval_syn_norm = 100. * P.eigVal_syn_norm__k / P.eigVal_syn_norm__k.sum()

plt.plot(PCs, eigval_obs[:maxPCs], 'k+--', label = '$F_{obs}$')
plt.plot(PCs, eigval_obs_norm[:maxPCs], 'k^-', label = '$f_{obs}$')
plt.plot(PCs, eigval_syn_norm[:maxPCs], 'k*-', label = '$f_{syn}$')
plt.legend()
plt.ylim([0, eigval_obs_norm[1] * 1.1])
plt.xlim([1, maxPCs])
plt.xticks(range(1,maxPCs + 1))
plt.title(r'%s - %s' % (K.galaxyName, K.califaID))
plt.xlabel(r'PC')
plt.ylabel(r'Var. [$\%%$]')
plt.grid()
f.savefig('%s/%s-screetest.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))

sqN_z = np.sqrt(P.K.N_zone)
a = sqN_z / len(P.l_obs)
ed__l = P.f_obs_norm__zl.std(axis = 0) * sqN_z
edk_mean = np.sqrt(np.abs(P.eigVal_obs_norm__k)).sum() * a
cd__l = np.abs(P.f_obs_norm__zl - P.f_obs_norm__zl.mean(axis = 0)).max(axis = 0)
k90_obs = ((np.cumsum(P.eigVal_obs_norm__k) / P.eigVal_obs_norm__k.sum()) <= 0.9).sum() + 1
k90_syn = ((np.cumsum(P.eigVal_syn_norm__k) / P.eigVal_syn_norm__k.sum()) <= 0.9).sum() + 1

#print '%s & %.5f & %.5f & %.5f & %.5f & %.5f & %d' % (K.califaID, a, ed__l.mean(), edk_mean, cd__l.mean(), (100. * P.eigVal_obs_norm__k[:6] / P.eigVal_obs_norm__k.sum())[:6].sum(), k90)
print '%s & %.5f & %.5f & %d & %d & %.5f & %.5f' % (K.califaID, 
                                                    (100. * P.eigVal_obs_norm__k[:6] / P.eigVal_obs_norm__k.sum())[:6].sum(), 
                                                    (100. * P.eigVal_syn_norm__k[:6] / P.eigVal_syn_norm__k.sum())[:6].sum(), 
                                                    k90_obs, k90_syn, 
                                                    P.eigVal_obs_norm__k.sum(), P.eigVal_syn_norm__k.sum(),
                                                   )

f, axArr = plt.subplots(3, 1)
f.set_size_inches((7, 9))

ax = axArr[0]
ed_mean = ed__l.mean()
ax.plot(P.l_obs, ed__l)
ax.axhline(y = ed_mean, color = 'k', ls='--')
ax.set_xlabel(r'$\lambda$')
ax.set_ylabel(r'Distancia Euclidiana')
ax.grid()
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax.text(P.l_obs[1], ed_mean, r'%.2f' % ed_mean,
        fontsize = 12,
        backgroundcolor='w',
        verticalalignment = 'center',
        horizontalalignment = 'center')

for xmin in ax.xaxis.get_minorticklocs():
    ax.axvline(x = xmin, ls = ':', c = 'grey')

ax = axArr[1]
ax.plot(P.l_obs, cd__l)
ax.axhline(y = cd__l.mean(), color = 'k', ls='--')
ax.set_xlabel(r'$\lambda$')
ax.set_ylabel(r'Distancia de Chebyshev')
ax.grid()
md_mean = cd__l.mean()
t = P.l_obs[0] + 0.1 * (P.l_obs[-1] - P.l_obs[0])
ax.text(P.l_obs[1], md_mean, r'%.2f' % md_mean,
        fontsize = 12,
        backgroundcolor='w',
        verticalalignment = 'center',
        horizontalalignment = 'center')

ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

for xmin in ax.xaxis.get_minorticklocs():
    ax.axvline(x = xmin, ls = ':', c = 'grey')

ax = axArr[2]
subpos = [0.7, 0.5, 0.26, 0.45]
ax.plot(P.eigVal_obs_norm__k)
ax.set_xlim([0, len(P.l_obs)])
ax.set_yscale('log')
ax.set_xlabel(r'autoespectro')
ax.set_ylabel(r'$\log\ \Lambda_k$')
subax = add_subplot_axes(ax, subpos)
subax.plot(P.eigVal_obs_norm__k, 'k+-')
subax.set_ylim([P.eigVal_obs_norm__k[0:maxPCs].min(), P.eigVal_obs_norm__k[0:maxPCs].max()])
subax.set_xlim([0, maxPCs + 1])
subax.set_yscale('log')
f.suptitle('%s' % K.califaID)
f.savefig('%s/%s-dist.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
