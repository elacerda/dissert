# coding=UTF-8
import sys
import pylab
from pylab import sqrt
from atpy import Table
import os

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin:/opt/local/bin'

filters = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z']

table_names = {}
table_names['FUV'] = 'galex1500.txt'
table_names['NUV'] = 'galex2500.txt'
table_names['u'] = 'u_SDSS.txt'
table_names['g'] = 'g_SDSS.txt'
table_names['r'] = 'r_SDSS.txt'
table_names['i'] = 'i_SDSS.txt'
table_names['z'] = 'z_SDSS.txt'

filter_line = {
               'FUV': 'b-',
               'NUV': 'r-',
               'u': 'b--',
               'g': 'g--',
               'r': 'y--',
               'i': 'r--',
               'z': 'k--'
}

# From http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
fig_width_pt = 448.07378
inches_per_pt = 1.0 / 72.27
golden_mean = (sqrt(5) - 1.0) / 2.0
fig_width = fig_width_pt * inches_per_pt
fig_height = fig_width * golden_mean * 0.85
fig_size = (fig_width, fig_height)
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family': 'serif',
          'figure.subplot.bottom': 0.12,
          'figure.figsize': fig_size}
pylab.rcParams.update(params)

pylab.figure(1)
pylab.axis([1100, 10000, 0, 0.7])
pylab.xlabel('Comprimento de onda [\AA]')
pylab.ylabel('Efici\^encia')
for f in filters:
    t = Table('filters/' + table_names[f], type='ascii')
    pylab.plot(t.col1, t.col2, filter_line[f], label=f)

pylab.legend()
#pylab.show()
#exit()
pylab.savefig('../doc/figuras/galex-filters.eps', format='eps')

