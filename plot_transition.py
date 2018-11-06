#!/usr/bin/env python
"""
Plot percentile bounds of 
"""
import numpy as np
import pylab as P
import wztanh as model
from matplotlib import ticker

# FIXME: Plot this out to z ~ 1090, to see how big a fraction DE can be at 
# the surface of last scattering

fnames = [ "new_wztanh_cmb_lss.dat", 
           "new_wztanh_cmb_lss_desi.dat", 
           #"new_wztanh_cmb_lss_desi_cosvis.dat",
           "new_wztanh_cmb_lss_desi_cosvis_3pbwedge.dat",
           "new_wztanh_cmb_lss_desi_cosvis_nowedge.dat" ]
#fnames = [ "new_mocker_cmb_lss.dat", 
#           "new_mocker_cmb_lss_desi.dat", 
#           #"new_mocker_cmb_lss_desi_cosvis.dat",
#           "new_mocker_cmb_lss_desi_cosvis_3pbwedge.dat",
#           "new_mocker_cmb_lss_desi_cosvis_nowedge.dat" ]

labels = [ 'CMB + LSS', 
           '  + DESI', 
           '  + DESI + Stage 2 (primary beam wedge)',
           '  + DESI + Stage 2 (no wedge)' ]

colours = [ '#FFE800', '#FF9F00', '#F43434', '#34CBF4' ]
colours = [ '#FFE28F', '#FFB45C', '#E76E6E', '#82D0F1' ]

colours = [ '#FFCB8F', '#E76E6E', '#C3E2FE', '#5291D1' ]

def load_chain(fname):
    """
    Load MCMC chain output from file
    """
    # Load header
    f = open(fname, 'r')
    hdr = f.readline()
    f.close()
    names = hdr[2:-1].split(' ')
    
    # Load data
    dat = np.genfromtxt(fname).T
    data = {}
    for i, n in enumerate(names):
        data[n] = dat[i]
    return data


# Load data and select random samples
for j, fn in enumerate(fnames):
    
    """
    try:
        pct = np.load("%s.pctcache.npy" % fn)
        print "%s: Loaded percentiles from cache" % fn
    except:
        # Load samples
        print "%s: Loading samples" % fn
        dat = load_chain(fn)
        
        ode = []
        for i in range(10000, dat['h'].size):
            pp = {key: dat[key][i] for key in dat.keys()}
            if 'mocker' in fn: pp['mocker'] = True # FIXME
            OmegaDE0 = 1. - pp['omegaM'] - model.omegaR(pp) - 0.
            ode.append( OmegaDE0 * model.omegaDE(a, pp) )
    """
    
    dat = load_chain(fn)
    P.hist(dat['zc'][10000:], bins=40, label=labels[j], alpha=0.4)
    
    # Plot 95% percentiles
    #lbl = labels[j]
    #P.fill_between(z, pct[0], pct[-1], color=colours[j], alpha=0.85, 
    #               label=lbl, lw=2.5)
    #P.plot(z, pct[2], ls='solid', color=colours[j], lw=1.8)


#P.plot(z, plcdm['omegaM'] * a**-3., 'k-', lw=1.8)
#P.plot(z, 0.10 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)
#P.plot(z, 0.01 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)

#P.xlim((0., 7.))
#P.ylim((0., 3.))

P.xlabel("$z_c$", fontsize=15)
#P.ylabel(r"$\rho_{\rm DE}(z) / \rho_{\rm crit}(z=0)$", fontsize=15, labelpad=10)

P.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

P.gca().tick_params(axis='both', which='major', labelsize=18, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=18, size=5., 
                    width=1.5, pad=8.)

leg = P.legend(loc='upper right', frameon=False)

P.tight_layout()
#P.savefig("pub_bounds_tanh.pdf")
#P.savefig("pub_bounds_mocker.pdf")
P.show()
