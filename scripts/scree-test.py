#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import PCAlifa as PCA

output_fmt = 'png'

fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

lc = [3800, 6850]
#lc = [3840, 6840]

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = lc)
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

var__z = ((P.f_obs_norm__zl - P.f_obs_norm__zl.mean(axis = 0)) ** 2.0).sum(axis = 1) / (K.N_zone - 1)
euclid_dist__z = np.sqrt(np.abs(var__z)) # desvio padrão

print '%s ed__z.mean():%.5f\teigVal__k.sum(): %.5f\tsqrt(eigVal_k).mean(): %.5f' % (K.califaID, euclid_dist__z.mean(), P.eigVal_obs_norm__k.sum(), (np.sqrt(np.abs(P.eigVal_obs_norm__k))).mean())

plt.plot(PCs, eigval_obs[:maxPCs], 'k+--', label = '$F_{obs}$')
plt.plot(PCs, eigval_obs_norm[:maxPCs], 'k^-', label = '$f_{obs}$')
plt.plot(PCs, eigval_syn_norm[:maxPCs], 'k*-', label = '$f_{syn}$')
plt.legend()
plt.ylim([0, eigval_obs_norm[1] * 1.1])
plt.xlim([1, maxPCs])
plt.xticks(range(1,maxPCs + 1))
plt.title(r'%s - %s' % (K.galaxyName, K.califaID))
#plt.title(r'%s - %s ($\Lambda\ =$ %.4f)' % (K.galaxyName, K.califaID, Lambda))
plt.xlabel(r'PC')
plt.ylabel(r'Var. [$\%%$]')
plt.grid()
f.savefig('%s/%s-screetest.%s' % (output_dir, K.califaID, output_fmt))
