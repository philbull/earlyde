#!/usr/bin/env python
"""
Plot percentile bounds on rho_DE(z).
"""
import numpy as np
import pylab as P
import wztanh as model
from load_data_files import load_chain
from matplotlib import ticker
import sys, os 

MODE = 'zmax'
#MODE = 'zmin'
#figname = "pub_bounds_cvallz.pdf"


if MODE == 'zmax':
    # zmax increasing
    tmpl = "chains/final_wztanh_seed200_cmb_lss_cvallz-%s"
    expts = ['%02d' % (i+1) for i in range(15)]
    zmin = 0.1
    zmax = np.arange(0.3, 5.91, 0.4)
    print(zmax)
    _z = zmax
else:
    # zmin increasing
    tmpl = "chains/final_wztanh_seed200_cmb_lss_cvallz-%s"
    expts = ['15', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    expts = [s.lower() for s in expts]
    zmin = np.arange(0.1, 5.91, 0.4)
    zmax = 5.90
    print(zmin)
    _z = zmin         
    

fnames = [tmpl % e for e in expts]
labels = expts

#    'CVALLZ-01':    "Fisher-full-gCVALLZ-0.10-0.30",
#    'CVALLZ-02':    "Fisher-full-gCVALLZ-0.10-0.70",
#    'CVALLZ-03':    "Fisher-full-gCVALLZ-0.10-1.10",
#    'CVALLZ-04':    "Fisher-full-gCVALLZ-0.10-1.50",
#    'CVALLZ-05':    "Fisher-full-gCVALLZ-0.10-1.90",
#    'CVALLZ-06':    "Fisher-full-gCVALLZ-0.10-2.30",
#    'CVALLZ-07':    "Fisher-full-gCVALLZ-0.10-2.70",
#    'CVALLZ-08':    "Fisher-full-gCVALLZ-0.10-3.10",
#    'CVALLZ-09':    "Fisher-full-gCVALLZ-0.10-3.50",
#    'CVALLZ-10':    "Fisher-full-gCVALLZ-0.10-3.90",
#    'CVALLZ-11':    "Fisher-full-gCVALLZ-0.10-4.30",
#    'CVALLZ-12':    "Fisher-full-gCVALLZ-0.10-4.70",
#    'CVALLZ-13':    "Fisher-full-gCVALLZ-0.10-5.10",
#    'CVALLZ-14':    "Fisher-full-gCVALLZ-0.10-5.50",
#    'CVALLZ-15':    "Fisher-full-gCVALLZ-0.10-5.90",
#    'CVALLZ-A':     "Fisher-full-gCVALLZ-0.50-5.90",
#    'CVALLZ-B':     "Fisher-full-gCVALLZ-0.90-5.90",
#    'CVALLZ-C':     "Fisher-full-gCVALLZ-1.30-5.90",
#    'CVALLZ-D':     "Fisher-full-gCVALLZ-1.70-5.90",
#    'CVALLZ-E':     "Fisher-full-gCVALLZ-2.10-5.90",
#    'CVALLZ-F':     "Fisher-full-gCVALLZ-2.50-5.90",
#    'CVALLZ-G':     "Fisher-full-gCVALLZ-2.90-5.90",
#    'CVALLZ-H':     "Fisher-full-gCVALLZ-3.30-5.90",
#    'CVALLZ-I':     "Fisher-full-gCVALLZ-3.70-5.90",
#    'CVALLZ-J':     "Fisher-full-gCVALLZ-4.10-5.90",
#    'CVALLZ-K':     "Fisher-full-gCVALLZ-4.50-5.90",
#    'CVALLZ-L':     "Fisher-full-gCVALLZ-4.90-5.90",
#    'CVALLZ-M':     "Fisher-full-gCVALLZ-5.30-5.90",
#    'CVALLZ-N':     "Fisher-full-gCVALLZ-5.70-5.90",


# Load data and select random samples
vals_z02 = []
vals_z1 = []
vals_z2 = []
vals_z3 = []
for j, fn in enumerate(fnames):
    
    # Custom redshift sampling
    z = np.concatenate(( [0.,], 
                         np.linspace(0.001, 1., 40), 
                         np.linspace(1.05, 3., 34), 
                         np.linspace(3.15, 10., 25) ))
    a = 1./(1.+z)
    pcts = [2.5, 16., 50., 84., 97.5]
    
    try:
        pct = np.load("%s.pctcache.npy" % fn)
        print("%s: Loaded percentiles from cache" % fn)
    except:
        raise ValueError("No cache of percentiles was found for %s" % fn)
    
    # Get 95% percentiles at different redshifts
    lbl = labels[j]
    print(lbl)
    
    # Redshifts
    vals_z02.append((pct[-1] - pct[0])[np.argmin((np.abs(z - 0.2)))])
    vals_z1.append((pct[-1] - pct[0])[np.argmin((np.abs(z - 1.)))])
    vals_z2.append((pct[-1] - pct[0])[np.argmin((np.abs(z - 2.)))])
    vals_z3.append((pct[-1] - pct[0])[np.argmin((np.abs(z - 3.)))])

# Plot values
col = '#865390'
P.plot(_z, vals_z02, color=col, marker='.', alpha=1.0, label="$z = 0.2$")
P.plot(_z, vals_z1, color=col, marker='.', alpha=0.7, label="$z = 1.0$")
P.plot(_z, vals_z2, color=col, marker='.', alpha=0.5, label="$z = 2.0$")
P.plot(_z, vals_z3, color=col, marker='.', alpha=0.2, label="$z = 3.0$")

P.xlim(0., 6.)
P.ylim(0., 0.5)

if MODE == 'zmax':
    P.text(3.8, 0.45, r"$z_{\rm min} = %1.1f$" % zmin, fontsize=14)
    P.xlabel(r"$z_{\rm max}$", fontsize=15)
else:
    P.text(3.8, 0.45, r"$z_{\rm max} = %1.1f$" % zmax, fontsize=14)
    P.xlabel(r"$z_{\rm min}$", fontsize=15)

#P.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#P.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
P.gca().xaxis.set_major_locator(ticker.MultipleLocator(1.0))
P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
P.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.02))

if MODE == 'zmin':
    P.legend(loc='upper left', frameon=False, fontsize=15)

P.ylabel(r"$\sigma(\rho_{\rm DE}(z) / \rho_{\rm crit}(z=0))$", fontsize=15, labelpad=10)

P.gca().tick_params(axis='both', which='major', labelsize=14, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=14, size=5., 
                    width=1.5, pad=8.)

P.tight_layout()
if MODE == 'zmax':
    P.savefig("pub_bounds_cvallz_zmax.pdf")
else:
    P.savefig("pub_bounds_cvallz_zmin.pdf")
P.show()
