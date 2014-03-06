#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import matplotlib.gridspec as gridspec
import argparse as ap
import PCAlifa as PCA

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
    parser.add_argument('--galaxyimgfile', '-g',
                        help = 'The image of the galaxy',
                        metavar = 'FILE',
                        type = str,
                        default = None)
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

output_fmt = args.outputimgsuffix
output_dir = args.outputdir
lc = args.lc
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

prop = {
    'arr'   : [ K.at_flux__yx, np.log10(K.aZ_flux__yx / 0.019), K.A_V__yx, K.v_0__yx, K.v_d__yx ],
    'label' : [ r'$\langle \log\ t rangle_L\ [yr]$', r'$\log\ \langle Z \rangle_L\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$' ],
    'name'  : [ 'at_flux', 'aZ_flux', 'AV', 'v0', 'vd' ]
}

f, axArr = plt.subplots(3, 3)
f.set_size_inches((15, 13))

for ax in f.axes:
    ax.set_axis_off()

galimg = plt.imread(args.galaxyimgfile)

ax = axArr[0,1]
ax.set_axis_on()
plt.setp(ax.get_xticklabels(), visible = False)
plt.setp(ax.get_yticklabels(), visible = False)

ax.imshow(galimg)

ax = axArr[1,0]
ax.set_axis_on()
fobs_norm__yx = K.zoneToYX(K.fobs_norm / 1.e-16, extensive = False)
ax.set_title(r'$\log\ F_{\lambda 5635}\ [10^{-16} erg/s/cm^2/\AA]$')
im = ax.imshow(np.log10(fobs_norm__yx), origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1,1]
ax.set_axis_on()
p_i = 0
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1,2]
ax.set_axis_on()
p_i = 1
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2,0]
ax.set_axis_on()
p_i = 2
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2,1]
ax.set_axis_on()
p_i = 3
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2,2]
ax.set_axis_on()
p_i = 4
ax.set_title(prop['label'][p_i])
prc = np.percentile(K.v_d, 98.)
print prc
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r', vmin = 0, vmax = prc)
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

plt.suptitle(r'%s - %s' % (K.galaxyName, K.califaID))
f.savefig('%s/%s-apresent.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
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

#xxx

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
#xxx

nPCs = 5
hfig = nPCs * 5

f = plt.figure(figsize = (15, hfig))
gs = gridspec.GridSpec(nPCs, 2, width_ratios = [5, 8])

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
#    plt.setp(ax2.get_xticklabels(), rotation = 45)
    ax2.set_ylabel(r'$\lambda$ [\AA]')
    ax2.grid()
    ax3 = ax2.twinx()
    ax3.plot(l, f_syn_mean, color = '0.65')
    ax3.set_ylabel(r'Espectro medio')
    ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

    for xmin in ax2.xaxis.get_minorticklocs():
        ax2.axvline(x = xmin, ls = ':', c = 'grey')

f.tight_layout()
f.savefig('%s/%s-tomo-syn-norm.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
