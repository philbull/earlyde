#!/usr/bin/env python
"""
Plot w(z) curves for a few different parameter values.
"""
import numpy as np
import pylab as P
import wztanh as model
from matplotlib import ticker
import matplotlib.gridspec as gridspec

MODE = 'tracker'
#MODE = 'mocker'

col_hrx = '#D92E1C'; col_hrx_light = '#FF9B99'
#col_hizrx = '#E3B4B0'; col_hizrx_light = '#E8C2BF'
#col_hizrx = '#81A871'; col_hizrx_light = '#81A871'
col_hizrx = col_hizrax_light = '#E9AD2B'
col_cv = '#00529B'; col_cv_light = '#85C5FF'

if MODE == 'tracker':
    colours = ['#D92E1Cff', '#D92E1C99', '#D92E1C44', '#D92E1C44']
else:
    colours = ['#00529Bff', '#00529Bbb', '#00529B88', '#00529B44']


z = np.linspace(0., 10., 500)
a = 1./(1. + z)

# Set base parameters
pbase = {
    'omegaB':   0.0493,
    'omegaM':   0.3153,
    'omegaK':   0.0,
    'h':        0.6736,
}

lcdm = {
    'w0':       -0.999999,
    'winf':     -0.999999,
    'zc':       2000.0,
    'deltaz':   0.5,
}

tracker1 = {
    'w0':       -0.9,
    'winf':     -0.7,
    'zc':       2.0,
    'deltaz':   1.5,
}
tracker2 = {
    'w0':       -0.9,
    'winf':     -0.7,
    'zc':       2.0,
    'deltaz':   0.5,
}
tracker3 = {
    'w0':       -0.9,
    'winf':     -0.7,
    'zc':       4.0,
    'deltaz':   0.5,
}
tracker4 = {
    'w0':       -0.9,
    'winf':     -1.1,
    'zc':       4.0,
    'deltaz':   0.5,
}

mocker1 = {
    'w0':       -0.9,
    'Cpow':     0.5,
    'mocker':   True
}
mocker2 = {
    'w0':       -0.9,
    'Cpow':     1.,
    'mocker':   True
}
mocker3 = {
    'w0':       -0.9,
    'Cpow':     2.,
    'mocker':   True
}
mocker4 = {
    'w0':       -0.9,
    'Cpow':     3.,
    'mocker':   True
}


# Collect parameter sets together
if MODE == 'tracker':
    paramsets = {
        '($z_c=2$, $\Delta z = 1.5$)':  tracker1,
        '($z_c=2$, $\Delta z = 0.5$)':  tracker2,
        '($z_c=4$, $\Delta z = 0.5$, $\Delta w = +0.2$)':  tracker3,
        '($z_c=4$, $\Delta z = 0.5$, $\Delta w = -0.2$)':  tracker4,
    }
else:
    paramsets = {
        '($C = 0.5$)':   mocker1,
        '($C = 1$)':     mocker2,
        '($C = 2$)':     mocker3,
        '($C = 3$)':     mocker4,
    }

# Get LCDM values
lcdm.update(pbase)
hz_lcdm = model.Hz(a, lcdm)
ode_lcdm = model.omegaDE(a, lcdm)
w_lcdm = model.wz(a, lcdm)


# Loop over parameter sets and plot observables
#ax1, ax2 = P.subplot(211), P.subplot(212)
gs = gridspec.GridSpec(2, 1, hspace=0.01, left=0.23, top=0.96, right=0.97)
ax1 = P.subplot(gs[0, 0])
ax2 = P.subplot(gs[1, 0])

for i, key in enumerate(paramsets):
    
    # Get parameter set
    p = paramsets[key]
    p.update(pbase)
    
    # Calculate observables
    hz = model.Hz(a, p)
    #ode = model.omegaDE(a, p)
    w = model.wz(a, p)
    
    # Plot observables
    ls = 'solid'
    if i == 3: ls = 'dashed'
    ax1.plot(z, w, lw=1.8, label=key, color=colours[i], ls=ls)
    ax2.plot(z, hz/hz_lcdm - 1., lw=1.8, label=key, color=colours[i], ls=ls)


#ax1.axhline(-1., lw=1.5, color='k', ls='dashed')
ax2.axhline(0., lw=1.5, color='k', ls='dashed')

# Common x axes
for ax in (ax1, ax2):
    ax.set_xlim((0., 9.1))
ax2.set_xlabel("$z$", fontsize=15)

ax1.set_ylabel(r"$w(z)$", fontsize=15, labelpad=10)
ax2.set_ylabel(r"$\Delta H(z) / H_{\Lambda{\rm CDM}}(z)$", 
               fontsize=15, labelpad=10)

if MODE == 'tracker':
    ax1.set_ylim((-1.2, -0.5))
    ax2.set_ylim((-0.01, 0.048))
else:
    ax1.set_ylim((-1.0, 0.15))
    ax2.set_ylim((-0.01, 0.117))

ax1.legend(loc='upper left', frameon=False, prop={'size': 11.5}, ncol=2)

for ax in (ax1, ax2):
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    
    ax.tick_params(axis='both', which='major', labelsize=16, size=8., 
                        width=1.5, pad=5.)
    ax.tick_params(axis='both', which='minor', labelsize=16, size=5., 
                        width=1.5, pad=5.)

ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax1.tick_params(axis='x', which='both', labelbottom=False, direction='in')

P.gcf().set_size_inches((7., 8.))
#P.tight_layout()

if MODE == 'tracker':
    P.savefig("pub_bkgd_tracker.pdf")
else:
    P.savefig("pub_bkgd_mocker.pdf")
P.show()
