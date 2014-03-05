#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt


fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

#lc = [3800, 6850]
lc = [3840, 6840]

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

K = P.K

#xxx
output_fmt = 'pdf'
nSpec = 5
mask = P.maskQFlag & P.maskEmLines & P.maskLambdaConstrains

f, axArr = plt.subplots(nSpec)
f.set_size_inches((10, 2 * nSpec))

for i, zone in enumerate(np.asarray([0, 1./5, 2./5, 3./5, 4./5]) * K.N_zone):
    z_i = np.int(zone)

    ax = axArr[i]
    ax.plot(K.l_obs, np.ma.masked_array(K.f_obs[:, z_i] / 1.e-16, mask = ~mask), 'k-')
    ax.set_ylabel(r'Zona %i' % z_i)
    ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
    ax.set_xlim(lc)
    ax.grid()

    for xmin in ax.xaxis.get_minorticklocs():
        ax.axvline(x = xmin, ls = ':', c = 'grey')

    plt.setp(ax.get_xticklabels(), visible = False)

f.suptitle(u'Exemplos de espectro da galáxia %s' % K.galaxyName)
ax = f.add_axes( [0., 0., 1, 1] )
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.text( 
    .02, 0.5, r'$F_{obs}\ [10^{-16} erg/s/cm^2/\AA]$', rotation='vertical',
    horizontalalignment='center', verticalalignment='center'
)

#f.subplots_adjust(left=0.07, bottom=0.1, top=0.95, wspace=0.2, hspace=0.1)

#f.subplots_adjust(hspace = 0.0)
plt.setp(axArr[nSpec - 1].get_xticklabels(), visible = True, rotation = 45)
f.savefig('%s/%s-exampleSpectra.%s' % (output_dir, K.califaID, output_fmt))

f, axArr = plt.subplots(1,2)
f.set_size_inches((10, 3))

mask2d = (mask[..., np.newaxis] & np.ones_like(K.f_obs, dtype = np.bool))

ax = axArr[0]
ax.fill_between(K.l_obs,
                np.ma.masked_array(K.f_obs / 1.e-16, mask = ~mask2d).max(axis = 1), 
                np.ma.masked_array(K.f_obs / 1.e-16, mask = ~mask2d).min(axis = 1), 
                edgecolor='gray', facecolor='lightgray')
ax.set_ylabel(r'$F_{obs}\ [10^{-16} erg/s/cm^2/\AA]$')
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
#ax.set_xlim([3800, 6850])
ax.grid()

#for xmin in ax.xaxis.get_minorticklocs():
#    ax.axvline(x = xmin, ls = ':', c = 'grey')


ax = axArr[1]
ax.fill_between(K.l_obs,
                np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d).max(axis = 1), 
                np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d).min(axis = 1), 
                edgecolor='gray', facecolor='lightgray')
ax.set_ylabel(r'$f_{obs}$')
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
#ax.set_xlim([3800, 6850])
ax.grid()

#for xmin in ax.xaxis.get_minorticklocs():
#    ax.axvline(x = xmin, ls = ':', c = 'grey')

f.subplots_adjust(left=0.07, bottom=0.1, top=0.95, wspace=0.2, hspace=0)
#plt.setp(axArr[0].get_xticklabels(), rotation = 45)
#plt.setp(axArr[1].get_xticklabels(), rotation = 45)
f.savefig('%s/%s-exampleSpectraFill.%s' % (output_dir, K.califaID, output_fmt))
