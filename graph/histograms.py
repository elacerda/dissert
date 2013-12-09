import sys
import os
import pylab
from pylab import sqrt
from matplotlib import cm
from matplotlib import gridspec
from atpy import Table
import numpy as np 

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin:/opt/local/bin'

debug = False
outformat = 'eps'


def set_eps_output_1():
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
              'title.fontsize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

def set_eps_output_5x2():
    # From http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    fig_width_pt = 448.07378
    inches_per_pt = 1.0 / 72.27
    golden_mean = (sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean * 2.0
    fig_size = (fig_width, fig_height)
    params = {'backend': 'ps',
              'axes.labelsize': 8,
              'axes.titlesize': 8,
              'text.fontsize': 8,
              'title.horizontalalignment': 'left',
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.subplot.hspace': .5,
              'figure.subplot.wspace': .15,
              'figure.subplot.left': 0.05,
              'figure.subplot.right': 0.95,
              'figure.subplot.top': 0.95,
              'figure.subplot.bottom': 0.05,
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

def randomFraction(arr, f):
    return (np.random.rand(len(arr)) < f)


tablename = sys.argv[1]
t = Table(tablename)
t.add_column('nii_ha', t.nii_6584_flux / t.halpha_flux)

sample_filter = (t.redshift < 0.17) & (t.redshift > 0.04)
sample_filter &= (t.m_r < 17.77)


# Completeness in redshift
bounds = [0.04, 0.17, -23.5, -18.5]

redshift = t.redshift[sample_filter]
Mr = t.r[sample_filter]

set_eps_output_1()
pylab.figure()
pylab.axis(bounds)
pylab.xlabel('Redshift')
pylab.ylabel('$M_r$')
pylab.hexbin(redshift, Mr, extent=bounds, gridsize=50, cmap=cm.Oranges)
pylab.colorbar()

if debug:
    pylab.show()
else:    
    pylab.savefig('../doc/figuras/completeness-volume.' + outformat, format=outformat)



# Histograms

cols = ['FUV', 'NUV', 'nii_ha', 'halpha_ew',
        'at_flux', 'at_mass', 'am_flux', 'am_mass', 'mcor_gal', 'AV']

filter = {'FUV': t.FUV <> -999,
          'NUV': t.NUV <> -999,
          'nii_ha': (t.nii_6584_flux > 0) & (t.halpha_flux > 0),
          'halpha_ew': t.halpha_ew > 0,
          'at_flux': (t.at_flux <> -999) ,
          'at_mass': t.at_mass <> -999,
          'am_flux': t.am_flux <> -999,
          'am_mass': t.am_mass <> -999,
          'mcor_gal': t.mcor_gal <> -999,
          'AV': t.AV <> -999
          }


op =     {'FUV':None,
          'NUV': None,
          'nii_ha': np.log10,
          'halpha_ew': np.log10,
          'at_flux': None,
          'at_mass': None,
          'am_flux': None,
          'am_mass': None,
          'mcor_gal': None,
          'AV': None
          }

range =  {'FUV': (-21.0, -12.0),
          'NUV':  (-21.0, -12.0),
          'nii_ha': (-1.0, 0.5),
          'halpha_ew': (-1.5, 2.0),
          'at_flux': (7.0, 10.25),
          'at_mass': (9.0, 10.25),
          'am_flux': (0, 2.5),
          'am_mass': (0, 2.5),
          'mcor_gal': (9.0, 12.0),
          'AV': (-0.5, 1.5)
          }

title =  {'FUV': '\\textbf{(a)}',
          'NUV': '\\textbf{(b)}',
          'at_flux': '\\textbf{(c)}',
          'at_mass': '\\textbf{(d)}',
          'am_flux': '\\textbf{(e)}',
          'am_mass': '\\textbf{(f)}',
          'mcor_gal': '\\textbf{(g)}',
          'AV': '\\textbf{(h)}',
          'nii_ha': '\\textbf{(i)}',
          'halpha_ew': '\\textbf{(j)}',
          }

label =  {'FUV': '$\\mathrm{FUV}$',
          'NUV': '$\\mathrm{NUV}$',
          'at_flux': '$\\langle\\log(t_{\\star}/[\\mathrm{a}])\\rangle_F$',
          'at_mass': '$\\langle\\log(t_{\\star}/[\\mathrm{a}])\\rangle_M$',
          'am_flux': '$\\langle Z_{\\star}/[Z_{\\odot}] \\rangle_F$',
          'am_mass': '$\\langle Z_{\\star}/[Z_{\\odot}] \\rangle_M$',
          'mcor_gal': '$\\log(M_{\\star}/[M_{\\odot}])$',
          'AV': '$A_V\\,[\\mathrm{mag}]$',
          'nii_ha': '$\\log([\mathrm{N\\,\\textsc{ii}}]/\\mathrm{H}\\alpha)$',
          'halpha_ew': '$\\log(W_{\\mathrm{H}\\alpha}/[\\mathrm{\\AA}])$',
          }

position = {'FUV': 521,
          'NUV': 522,
          'at_flux': 523,
          'at_mass': 524,
          'am_flux': 525,
          'am_mass': 526,
          'mcor_gal': 527,
          'AV': 528,
          'nii_ha': 529,
          'halpha_ew': 520,
          }


set_eps_output_5x2()
pylab.figure()
for col in cols:
    pylab.subplot(position[col])
    if op[col] is not None:
        x = op[col](t.data[col])
    else:
        x = t.data[col]
    x1 = x[filter[col] & sample_filter]
    h1, e1 = pylab.histogram(x1, 50, range[col])
    pylab.title(title[col], ha='left', x=0.05, y=0.75)
    pylab.xlabel(label[col])
    pylab.ylim(0.0, 1.1)
    pylab.setp(pylab.gca(), 'yticklabels', [])
    pylab.plot(e1[:-1], h1*1.0/np.max(h1), 'k-')

if debug:
    pylab.show()
else:    
    pylab.savefig('../doc/figuras/histogram-sample.' + outformat, format=outformat)

