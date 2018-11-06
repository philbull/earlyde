#!/usr/bin/env python
"""
Plot percentile bounds on rho_DE(z) for HIRAX and DESI.
"""
import numpy as np
import pylab as P
import wztanh as model
from matplotlib import ticker

# FIXME: Plot this out to z ~ 1090, to see how big a fraction DE can be at 
# the surface of last scattering

fnames = [ "new_wztanh_cmb_lss.dat", 
           "new_wztanh_cmb_lss_desi_hetdex.dat", 
           #"new_wztanh_cmb_lss_newexpt_horizwedge.dat",
           #"new_wztanh_cmb_lss_newexpt_3pbwedge.dat",
           #"new_wztanh_cmb_lss_newexpt_nowedge.dat",
           "new_wztanh_cmb_lss_desi.dat", ]
#fnames = [ "new_mocker_cmb_lss.dat", 
#           "new_mocker_cmb_lss_desi.dat", 
#           #"new_mocker_cmb_lss_desi_cosvis.dat",
#           "new_mocker_cmb_lss_desi_cosvis_3pbwedge.dat",
#           "new_mocker_cmb_lss_desi_cosvis_nowedge.dat" ]

labels = [ 'CMB + LSS', 
           '  + HETDEX + DESI',
           #'  + HIRAX (horizon wedge)',
          #r'  + HIRAX ($3\times$ PB wedge)', 
           #'  + HIRAX (no wedge)',
           '  + DESI' ]

#colours = [ '#FFE800', '#FF9F00', '#F43434', '#34CBF4' ]
#colours = [ '#FFE28F', '#FFB45C', '#E76E6E', '#82D0F1' ]

#colours = [ '#FFCB8F', 'g', '#E64D3D', '#5291D1', '#555555' ]
colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00', '#555555' ]

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
        data[n] = dat[i][20000:] # trim burn-in
    return data


# Load data and select random samples
for j, fn in enumerate(fnames):
    
    z = np.linspace(0., 10., 100)
    a = 1./(1.+z)
    pcts = [2.5, 16., 50., 84., 97.5]
    
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
        
        # Get percentiles of w(z) or OmegaDE(z) in each z bin
        #wz = [model.wz(a, {key: dat[key][i] for key in dat.keys()}) 
        #      for i in range(10000, dat['h'].size)] # dat['h'].size
        print dat['h'].size
        pct = np.percentile(ode, pcts, axis=0)
        
        # Save results
        np.save("%s.pctcache" % fn, pct)
    
    # Plot 95% percentiles
    lbl = labels[j]
    if j == 0:
        P.fill_between(z, pct[0], pct[-1], color=colours[j], alpha=0.6, 
                       label=lbl, lw=2.5)
    else:
        dashes = [1.5,1.5] if 'DESI' in lbl else []
        P.plot(z, pct[0], color=colours[j], lw=2.5, label=lbl, dashes=dashes)
        P.plot(z, pct[-1], color=colours[j], lw=2.5, dashes=dashes)

# LCDM curves
plcdm = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
                    w0=-1., winf=-1., zc=1e5, deltaz=0.5)
OmegaDE0 = 1. - plcdm['omegaM'] - model.omegaR(plcdm) - 0.
P.plot(z, OmegaDE0*model.omegaDE(a, plcdm), 'k-', lw=1.8, alpha=1.)

#P.plot(z, plcdm['omegaM'] * a**-3., 'k-', lw=1.8)
#P.plot(z, 0.10 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)
#P.plot(z, 0.01 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)

#P.text(0.1, 0.17, r"$1\% \times \rho_M(z)/\rho_{{\rm crit}, 0}$", fontsize=13, alpha=0.85)

P.xlim((0., 7.))
P.ylim((0., 1.8))

P.xlabel("$z$", fontsize=15)
P.ylabel(r"$\rho_{\rm DE}(z) / \rho_{\rm crit}(z=0)$", fontsize=15, labelpad=10)

P.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

P.gca().tick_params(axis='both', which='major', labelsize=18, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=18, size=5., 
                    width=1.5, pad=8.)

# Re-order legend
handles, labels = P.gca().get_legend_handles_labels()
handles = [handles[-1], handles[0], handles[1], handles[2], handles[3]]
labels = [labels[-1], labels[0], labels[1], labels[2], labels[3]]

leg = P.legend(handles, labels, loc='upper left', frameon=False)

P.tight_layout()
#P.savefig("pub_bounds_tanh_hirax.pdf")
#P.savefig("pub_bounds_mocker.pdf")
P.show()
