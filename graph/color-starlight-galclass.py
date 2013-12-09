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
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.subplot.bottom': 0.12,
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

def set_eps_output_3x2():
    # From http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    fig_width_pt = 448.07378
    inches_per_pt = 1.0 / 72.27
    golden_mean = (sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean * 1.5
    fig_size = (fig_width, fig_height)
    params = {'backend': 'ps',
              'axes.labelsize': 8,
              'text.fontsize': 8,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.subplot.hspace': .2,
              'figure.subplot.wspace': .2,
              'figure.subplot.left': 0.05,
              'figure.subplot.right': 0.95,
              'figure.subplot.top': 0.90,
              'figure.subplot.bottom': 0.05,
              'figure.figsize': fig_size}
    pylab.rcParams.update(params)

def randomFraction(arr, f):
    return (np.random.rand(len(arr)) < f)


tablename = sys.argv[1]
t = Table(tablename)

sample = (t.redshift > 0.04) 
sample &= (t.redshift < 0.17)
sample &= (t.m_r < 17.77)

NUV_r = t.NUV[sample] - t.r[sample]
g_r = t.g[sample] - t.r[sample]


logN2Ha = np.log10(t.nii_6584_flux[sample] / t.halpha_flux[sample])
WHa = t.halpha_ew[sample]
WN2 = t.nii_6584_ew[sample]

types = ['starforming',
         'sAGN',
         'wAGN',
         'retired',
         'passive'
]

typemask = {
            'starforming': (logN2Ha < -0.4) & (WHa > 3.0),
            'sAGN': (logN2Ha > -0.4) & (WHa > 6.0),
            'wAGN': (logN2Ha > -0.4) & (WHa < 6.0) & (WHa > 3.0),
            'retired': (WHa < 3.0) & (WHa > 0.5) & (WN2 > 0.5),
            'passive': (WHa < 0.5) & (WN2 < 0.5)
            }

position = {
            'starforming': 322,
            'sAGN': 323,
            'wAGN': 324,
            'retired': 325,
            'passive': 326
}

color = {
         'starforming': 'blue',
         'sAGN': 'lightgreen',
         'wAGN': 'darkgreen',
         'retired': 'black',
         'passive': 'red'
         }

label = {
         'starforming': '{\\bf (b)} SF',
         'sAGN': '{\\bf (c)} sAGN',
         'wAGN': '{\\bf (d)} wAGN',
         'retired': '{\\bf (e)} RG',
         'passive': '{\\bf (f)} PG'
         }

label2 = {
         'starforming': 'SF',
         'sAGN': 'sAGN',
         'wAGN': 'wAGN',
         'retired': 'RG',
         'passive': 'PG'
         }


bounds = [.5,7,.2,.9]

# Color - UV color for all the galaxy classes.
set_eps_output_1()
pylab.figure()
for type in types:
    mask = typemask[type]
    if type == 'all': continue
    # Plot less points
    fraction = randomFraction(NUV_r[mask], 0.1)
    pylab.axis(bounds)
    pylab.xlabel('$\\mathrm{NUV} - r$')
    pylab.ylabel('$g - r$')
    pylab.scatter(NUV_r[mask][fraction], g_r[mask][fraction],
                  c=color[type], marker='o', edgecolor='None', s=1, label=label2[type])
    
pylab.legend(loc='lower right', markerscale=5, numpoints=1)

if debug:
    pylab.show()
else:    
    pylab.savefig('../doc/figuras/uvcolor-color-class.' + outformat, format=outformat)


# mcor_gal : massa em estrelas
# at_flux: idade ponderada em fluxo
# at_mass: idade em massa
# am_flux: metalicidade em fluxo
# am_mass: metalicidade em massa
# AV: extincao
vmin = {}
vmax = {}
vmin_h = {}
vmax_h = {}
param = {}

vmin['at_flux'] = 7.5
vmax['at_flux'] = 10.0
vmin_h['at_flux'] = 7.0
vmax_h['at_flux'] = 10.5
param['at_flux'] = 'M\\\'edia do logaritmo da idade estelar [a] ponderada em fluxo'

vmin['at_mass'] = 9.4
vmax['at_mass'] = 10.2
vmin_h['at_mass'] = 9.0
vmax_h['at_mass'] = 10.4
param['at_mass'] = 'M\\\'edia do logaritmo da idade estelar [a] ponderada em massa'

vmin['am_flux'] = 0.0
vmax['am_flux'] = 2.0
vmin_h['am_flux'] = 0.0
vmax_h['am_flux'] = 2.0
param['am_flux'] = 'Metalicidade estelar m\\\'edia [$Z_{\odot}$] ponderada em fluxo'

vmin['am_mass'] = 0.0
vmax['am_mass'] = 2.0
vmin_h['am_mass'] = 0.0
vmax_h['am_mass'] = 2.0
param['am_mass'] = 'Metalicidade estelar m\\\'edia [$Z_{\odot}$] ponderada em massa'

vmin['mcor_gal'] = 9.0
vmax['mcor_gal'] = 12.0
vmin_h['mcor_gal'] = 9.0
vmax_h['mcor_gal'] = 12.5
param['mcor_gal'] = 'Logaritmo da massa estelar [$M_{\odot}$]'

vmin['AV'] = -0.25
vmax['AV'] = 1.0
vmin_h['AV'] = -0.4
vmax_h['AV'] = 1.2
param['AV'] = 'Extin\\c{c}\\~ao por poeira [magnitude]'

vmin['halpha_ew'] = -1.0
vmax['halpha_ew'] = 2.0
vmin_h['halpha_ew'] = -1.5
vmax_h['halpha_ew'] = 2.5
param['halpha_ew'] = 'Logaritmo da largura equivalente de $\\mathrm{H}_{\\alpha}$ [\AA]'

# Hack: plot log(EWHa)
t.halpha_ew = np.log10(t.halpha_ew)

#if debug:
#    for col in label.keys():
#        pylab.figure()
#        #pylab.hist(t.data[col][sample], range=(vmin[col], vmax[col]))
#        pylab.hist(t.data[col][sample])
#        pylab.title(col)
#        pylab.show()
#    #exit()

set_eps_output_3x2()

for p in param.keys():
    pylab.figure()
    pylab.suptitle(param[p])
    ax = pylab.subplot(321)
    for type in types:
        x = t.data[p][sample][typemask[type]]
        h, e = pylab.histogram(x, 40, (vmin_h[p], vmax_h[p]))
        pylab.title('\\textbf{(a)} Histograma', fontsize=8, ha='left', x=0.06, y=0.83)
        pylab.ylim(0.0, 1.2)
        pylab.setp(pylab.gca(), 'yticklabels', [])
        pylab.plot(e[:-1], h*1.0/np.max(h), color=color[type])
        
    for type in types:
        mask = typemask[type]
        h, ex, ey = np.histogram2d(NUV_r[mask], g_r[mask], bins=30,
                                   range=[[bounds[0], bounds[1]],[bounds[2], bounds[3]]])

        # Plot less points
        fraction = randomFraction(NUV_r[mask], 0.05)

        z = t.data[p][sample]
        ax = pylab.subplot(position[type])
        pylab.axis(bounds)
        pylab.text(4.75, .25, label[type], fontsize=8)
        pylab.scatter(NUV_r[mask][fraction], g_r[mask][fraction], c=z[mask][fraction],
            marker='o', edgecolor='None', s=1, cmap=cm.spectral,
            vmin=vmin[p], vmax=vmax[p], zorder=1)
        pylab.colorbar()

        pylab.contour(h.T, extent=bounds, colors='black', linewidths=0.5)

    if debug:
        pylab.show()
    else:    
        pylab.savefig('../doc/figuras/uvcolor-color-' + p + '-byclass.' + outformat, format=outformat)
