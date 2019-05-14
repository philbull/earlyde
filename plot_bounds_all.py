#!/usr/bin/env python
"""
Plot percentile bounds on rho_DE(z) for HIRAX and DESI, now with HETDEX included.
"""
import numpy as np
import pylab as P
import wztanh as model
from load_data_files import load_chain
from matplotlib import ticker
import sys, os 

# Get arguments
if len(sys.argv) != 3:
    print("Requires two arguments: idx [full|lowz]")
    sys.exit(1)
idx = int(sys.argv[1])
MODE = str(sys.argv[2]).lower()
if MODE not in ['lowz', 'full']:
    print("Second argument must be 'full' or 'lowz'")
    sys.exit(1)

legend = True

if idx == 0:
    tmpl = "chains/final_wztanh_cmb_lss%s.dat"
    expts = ["", "_desi", "_desi_hetdex", "_hirax_pbw", "_hetdex_hirax_pbw", "_cvlowz",]
    colours = [ '#E1E1E1', '#555555', '#863AB7', '#E6773D', '#F5BC00', 'g']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
               '  + DESI',
               '  + HETDEX + DESI',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + HETDEX + HIRAX ($3\times$ PB wedge)',
               '  + CV-limited low-z' ]

if idx == 1:
    tmpl = "chains/final_wztanh_cmb_lss%s.dat"
    expts = ["", "_hirax_hw", "_hirax_pbw", "_hirax_nw", "_cvlowz"]
    colours = [ '#E1E1E1', '#555555', '#863AB7', '#E6773D', '#F5BC00', 'g' ]
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + HIRAX (horizon wedge)',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + HIRAX (no wedge)',
               '  + CV-limited low-z' ]

if idx == 2:
    tmpl = "chains/final_wztanh_cmb_lss%s.dat"
    fnames = [ #"chains/final_wztanh_cmb_lss.dat", 
               #"chains/final_wztanh_cmb_lss_cvlowz.dat", 
               #"chains/final_wztanh_seed16_cmb_lss_cvlowz.dat",
               #"chains/final_wztanh_seed17_cmb_lss_cvlowz.dat",
               #"chains/final_wztanh_seed18_cmb_lss_cvlowz.dat",
               #"chains/final_wztanh_seed19_cmb_lss_cvlowz.dat",
               "chains/final_wztanh_seed20_cmb_lss_cvlowz.dat",
               "chains/final_wztanh_seed21_cmb_lss_cvlowz.dat",
             ]
    colours = [ '#E1E1E1', '#555555', '#863AB7', '#E6773D', '#F5BC00', 'g', 'm']
    labels = [ #'CMB + LSS', 
               #'  + CV-limited low-z (15)',
               #'  + CV-limited low-z (16)',
               #'  + CV-limited low-z (17)',
               #'  + CV-limited low-z (18)',
               #'  + CV-limited low-z (19)',
               '  + CV-limited low-z (20)',
               '  + CV-limited low-z (21)',
             ]

if idx == 3:
    tmpl = "chains/final_wztanh_seed21_cmb_lss%s"
    expts = ["", "_desi", "_desi_hetdex",
             #"_hirax_pbw", 
             "_hetdex_hirax_pbw", 
             "_cvlowz",]
    colours = [ '#E1E1E1', '#555555', '#863AB7', '#E6773D', '#F5BC00', 'g']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
               '  + DESI',
               '  + HETDEX + DESI',
              #r'  + HIRAX ($3\times$PB wedge)',
              r'  + HETDEX + HIRAX ($3\times$ PB wedge)',
               '  + CV-limited low-z' ]

if idx == 4:
    tmpl = "chains/final_wztanh_seed21_cmb_lss%s"
    expts = [ #"", 
             #"_desi",
             "_hirax_nw", "_hirax_pbw", "_hirax_hw",]
    colours = [ '#E1E1E1', '#555555', '#863AB7', '#E6773D', '#F5BC00', 'g']
    fnames = [tmpl % e for e in expts]
    labels = [ #'CMB + LSS', 
               #'  + DESI',
              r'  + HIRAX (no wedge)',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + HIRAX (horiz. wedge)',
              ]

if idx == 5:
    figname = "pub_bounds_tanh_hirax_%s.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/final_wztanh_seed21_cmb_lss%s"
    expts = [ "", 
             "_hirax_hw", "_hirax_pbw", "_hirax_nw",]
    colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + HIRAX (horiz. wedge)',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + HIRAX (no wedge)',
              ]

if idx == 6:
    #figname = "pub_bounds_tanh_hirax_%s.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/final_wztanh_seed21_cmb_lss%s"
    expts = [ "", 
             '_cosvis_nw', '_cosvis_pbw', '_cosvis_hw', 
             '_hizrax_nw', '_hizrax_pbw', '_hizrax_hw',]
    colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00', 'm', 'c', 'y']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + CosVis (no wedge)',
              r'  + CosVis ($3\times$PB wedge)',
              r'  + CosVis (horiz. wedge)',
              r'  + HIRAX high-z (no wedge)',
              r'  + HIRAX high-z ($3\times$PB wedge)',
              r'  + HIRAX high-z (horiz. wedge)',
              ]


# Load data and select random samples
for j, fn in enumerate(fnames):
    
    if "_hizrax_hw" in fn:
        fn = fn.replace("21", "23")
        
    
    #z = np.linspace(0., 10., 800)
    #if '20' in fn or '21' in fn: 
    #    z = np.linspace(0., 10., 150)
    #z = np.concatenate( ([0.,], np.logspace(-3., 1., 99)) )
    
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
        # Load samples
        print("%s: Loading samples" % fn)
        dat = load_chain(fn)
        
        ode = []
        for i in range(200000, dat['h'].size, 2):
            pp = {key: dat[key][i] for key in dat.keys()}
            if 'mocker' in fn: pp['mocker'] = True # FIXME
            OmegaDE0 = 1. - pp['omegaM'] - model.omegaR(pp) - 0.
            ode.append( OmegaDE0 * model.omegaDE(a, pp) )
        
        # Get percentiles of w(z) or OmegaDE(z) in each z bin
        #wz = [model.wz(a, {key: dat[key][i] for key in dat.keys()}) 
        #      for i in range(10000, dat['h'].size)] # dat['h'].size
        print(dat['h'].size)
        pct = np.percentile(ode, pcts, axis=0)
        
        # Save results
        np.save("%s.pctcache" % fn, pct)
    
    # Plot 95% percentiles
    lbl = labels[j]
    print(lbl)
    if j == 0:
        P.fill_between(z, pct[0], pct[-1], color=colours[j], alpha=0.6, 
                       label=lbl, lw=2.5)
    else:
        dashes = [1.5,1.5] if j == 1 else []
        P.plot(z, pct[0], color=colours[j], lw=2.5, label=lbl, dashes=dashes)
        P.plot(z, pct[-1], color=colours[j], lw=2.5, dashes=dashes)

# LCDM curves
plcdm = model.pdict(h=0.6727, omegaM=0.3166, omegaK=0.0, omegaB=0.04941, 
                    w0=-1., winf=-1., zc=1e5, deltaz=0.5)
OmegaDE0 = 1. - plcdm['omegaM'] - model.omegaR(plcdm) - 0.
P.plot(z, OmegaDE0*model.omegaDE(a, plcdm), 'k-', lw=1.8, alpha=1.)

#P.plot(z, plcdm['omegaM'] * a**-3., 'k-', lw=1.8)
#P.plot(z, 0.10 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)
#P.plot(z, 0.01 * plcdm['omegaM'] * a**-3., 'k--', lw=1.8, alpha=0.4)

#P.text(0.1, 0.17, r"$1\% \times \rho_M(z)/\rho_{{\rm crit}, 0}$", fontsize=13, alpha=0.85)


if MODE == 'lowz':
    P.xlim((0., 2.2))
    P.ylim((0.5, 0.9))
    
    P.gca().xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
else:
    P.xlim((0., 7.))
    P.ylim((0., 1.8))
    
    P.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

P.xlabel("$z$", fontsize=15)
P.ylabel(r"$\rho_{\rm DE}(z) / \rho_{\rm crit}(z=0)$", fontsize=15, labelpad=10)

P.gca().tick_params(axis='both', which='major', labelsize=18, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=18, size=5., 
                    width=1.5, pad=8.)

# Re-order legend
#handles, labels = P.gca().get_legend_handles_labels()
#handles = [handles[-1], handles[0], handles[1], handles[2], handles[3]]
#labels = [labels[-1], labels[0], labels[1], labels[2], labels[3]]
#leg = P.legend(handles, labels, loc='upper left', frameon=False)
if legend:
    # Re-order legend so CMB+LSS is always first
    _handles, _labels = P.gca().get_legend_handles_labels()
    handles = [_handles[-1],]
    labels = [_labels[-1],]
    for i in range(0, len(_handles)-1):
        handles += [_handles[i],]
        labels += [_labels[i],]
    leg = P.legend(handles, labels, loc='upper left', frameon=False)
    

P.tight_layout()
#P.savefig("pub_bounds_tanh_hetdex.pdf")
#P.savefig("pub_bounds_mocker.pdf")
try:
    P.savefig(figname)
except:
    pass
P.show()
