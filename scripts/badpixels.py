#!/usr/bin/python
import PCAlifa as PCA
import numpy as np
from matplotlib import pyplot as plt

P = PCA.PCAlifa(fitsFile = '/home/lacerda/CALIFA/gal_fits/K0277/K0277_synthesis_eBR_v20_q036.d13c512.ps03.k2.mC.CCM.Bgsd61.fits', quantilQFlag = 0.9, lc = [3800,6850])

K = P.K

f_flag_masked = K.f_flag[P.maskLambdaConstrains, :]

for i in range(K.N_zone):
    badpix = np.where(f_flag_masked[:,i] > 0)[0]
    # Force badpix flag = 1
    f_flag_masked[badpix, i] = 1.

histo = f_flag_masked.sum(axis = 1) / K.N_zone

plt.figure(figsize = (20, 11.25), dpi = 100)
plt.title('Bad Pixels')
plt.plot(K.l_obs[P.maskLambdaConstrains], histo, label = 'Correct masked intervals')
plt.plot(K.l_obs, P.histo, label = 'Today masked intervals')
plt.axhline(y = 0.9, color = 'r', lw = 1.0, ls = '--', label = '0.9 quantil')
plt.legend()
plt.savefig('badpix.png')
