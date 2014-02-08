#!/usr/bin/python
import sys
import os
import pylab
from pylab import sqrt
from matplotlib import cm
from atpy import Table
import numpy as np

outformat = 'pdf'

def set_eps_output_1():
    # From http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    fig_width_pt = 448.07378
    inches_per_pt = 1.0 / 72.27
    golden_mean = (sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean * 0.85
    fig_size = (fig_width, fig_height)
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 8,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.subplot.hspace': .5,
              'figure.subplot.bottom': 0.12,
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

def set_eps_output_2x1():
    # From http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    fig_width_pt = 448.07378
    inches_per_pt = 1.0 / 72.27
    golden_mean = (sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean
    fig_size = (fig_width, fig_height)
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 8,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.subplot.hspace': .5,
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

# FIXME: take only ten percent of points
def oneInTen(arr):
    return (np.random.rand(len(arr)) < 0.1) & arr

set_eps_output_1()
# Carregar arquivo FITS com os dados.
from pycasso import fitsQ3DataCube
K = fitsQ3DataCube('/home/lacerda/CALIFA/gal_fits/K0277/K0277_synthesis_eBR_v20_q036.d13c512.ps03.k2.mC.CCM.Bgsd61.fits')

# Converter zonas para imagem.
at_image = K.zoneToYX(K.at_flux__z, extensive = False)

# Desenhar o mapa.
import matplotlib.pyplot as plt
pylab.imshow(at_image, origin = 'lower', interpolation = 'nearest')
cb = pylab.colorbar()
cb.set_label(r'$\langle \log\ t\rangle_L [anos]$')
pylab.xlabel(r'Pixels')
pylab.title(r'NGC\ 2916')
pylab.savefig('../figuras/at_flux_zone.' + outformat, format = outformat)

bins = np.arange(0, 26, 1)
bin_center = 0.5 * (bins[1:] + bins[:-1])
at_rad = K.radialProfile(at_image, bins, rad_scale = 1.0)

pylab.clf()
pylab.plot(bin_center, at_rad)
pylab.ylabel(r'$\langle \log\ t\rangle_L [anos]$')
pylab.xlabel(r'radius [arcsec]')
pylab.title(r'NGC\ 2916')
pylab.savefig('../figuras/at_flux_radprof.' + outformat, format = outformat)
