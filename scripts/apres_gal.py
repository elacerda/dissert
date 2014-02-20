#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import PCAlifa as PCA

galimg_dir = '/home/lacerda/CALIFA/images'
galimg_fmt = 'jpg'

output_fmt = 'png'

fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = [3800, 6850])
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

K = P.K

prop = {
    'arr'   : [ K.at_flux__yx, np.log10(K.aZ_flux__yx / 0.019), K.A_V__yx, K.v_0__yx, K.v_d__yx ],
    'label' : [ r'$\log\ t\ [yr]$', r'$\log\ Z\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$' ],
    'name'  : [ 'at_flux', 'aZ_flux', 'AV', 'v0', 'vd' ]
}

f, axArr = plt.subplots(3, 3)
f.set_size_inches((15, 13))

for ax in f.axes:
    ax.set_axis_off()

galimg = plt.imread('%s/%s.%s' % (galimg_dir, K.califaID, galimg_fmt))

ax = axArr[0,1]
ax.set_axis_on()
plt.setp(ax.get_xticklabels(), visible = False)
plt.setp(ax.get_yticklabels(), visible = False)

ax.imshow(galimg)

ax = axArr[1,0]
ax.set_axis_on()
fobs_norm__yx = K.zoneToYX(K.fobs_norm / 1.e-16, extensive = False)
ax.set_title(r'$F_{\lambda 5635}\ [10^{-16} erg/s/cm^2/\AA]$')
im = ax.imshow(fobs_norm__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
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
f.savefig('%s/%s-apresent.%s' % (output_dir, K.califaID, output_fmt))
