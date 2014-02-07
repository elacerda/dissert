#!/usr/bin/python
import sys
import numpy as np
from pycasso import fitsQ3DataCube

fitsfile = sys.argv[1]

K = fitsQ3DataCube(fitsfile)

total = K.zoneArea_pix.size
onepix = np.where(K.zoneArea_pix == 1)[0].size
morethanten = np.where(K.zoneArea_pix > 10)[0].size

onepixperc = onepix / np.double(total)
morethantenperc = morethanten / np.double(total)

print '%s & %s & $%d$ & $%d$ & $%.2f$ & $%d$ & $%.2f$ \\\\' % (K.galaxyName, K.califaID, total, onepix, onepixperc, morethanten, morethantenperc)
