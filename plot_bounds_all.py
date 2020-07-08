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

# Get arguments
if len(sys.argv) != 3:
    print("Requires two arguments: idx [full|lowz]")
    sys.exit(1)
idx = int(sys.argv[1])
MODE = str(sys.argv[2]).lower()
if MODE not in ['lowz', 'full']:
    print("Second argument must be 'full' or 'lowz'")
    sys.exit(1)

SAMPLE_SEED = 22
BURNIN = 80 * 500 # Leaves 80 * 1500 samples
#NSAMP = 80 * 1500
NSAMP = 80 * 1000

legend = True

# Default dash pattern
dash = [False for j in range(10)]
dash[1] = True
LW = 2.5

if idx == 0:
    tmpl = "chains/final_wztanh_seed30_cmb_lss%s"
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
    tmpl = "chains/final_wztanh_cmb_lss%s"
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
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ "", 
             #"_desi",
             "_hirax_nw", "_hirax_pbw", "_hirax_hw",]
    colours = [ '#E1E1E1', '#863AB7', '#E6773D', '#F5BC00',] # '#555555'
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
               #'  + DESI',
              r'  + HIRAX (no wedge)',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + HIRAX (horiz. wedge)',
              ]

if idx == 5:
    figname = "pub_bounds_tanh_hirax_%s.200.pdf" % MODE
    legend = True
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
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
    #figname = "pub_bounds_tanh_hirax_%s.200.pdf" % MODE
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

if idx == 7:
    figname = "pub_bounds_tanh_compare_desi_%s.200.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ '', 
              '_desi_hetdex',
              '_hirax_pbw',
              #'_hizrax_pbw',
              '_cosvis_pbw',
              #'_cvlowz',
              #'_cvlowz_hetdex',
              #'_cvlowz_hizrax_pbw',
              ]
             
    colours = [ '#E1E1E1', '#555753', '#E6773D', '#3465a4']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + DESI + HETDEX',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + Stage 2 ($3\times$PB wedge)',
              #r'  + CV lim. (low-z)',
              #r'  + CV lim. (low-z) + HETDEX',
              #r'  + CV lim. (low-z) + HIRAX (high-z)',
              ]

if idx == 8:
    figname = "pub_bounds_tanh_compare_galsurv_%s.200.pdf" % MODE
    #if MODE == 'lowz': legend = False
    
    dash = [False for j in range(10)]
    dash[1] = True
    LW = 2.1
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ '', 
              #'_desi',
              '_cvlowz',
              '_cvlowz_hetdex',
              #'_cvlowz_hizrax_pbw',
              '_cvlowz_desi',
              '_cvlowz_cosvis_pbw',
              ]
             
    colours = [ '#E1E1E1', '#64a164', '#555753', '#fcaf3e', '#ad7fa8', '#3465a4',]
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS',
              #r'  + DESI', 
              r'  + CV-lim. (low-z)',
              r'  + CV-lim. (low-z) + HETDEX',
              #r'  + CV-lim. (low-z) + HIRAX high-z ($3\times$PB)',
              r'  + CV-lim. (low-z) + DESI',
              r'  + CV-lim. (low-z) + Stage 2 ($3\times$PB)',
              ]


if idx == 9:
    figname = "pub_bounds_tanh_compare_lowzhighz_%s.200.pdf" % MODE
    #if MODE == 'lowz': legend = False
    YLIM = (0.4, 1.1)
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ '', 
              '_cvlowz',
              '_spectel',
              '_cosvis_pbw',
              '_cvlowz_cosvis_pbw',
              ]
             
    colours = [ '#E1E1E1', 'r', '#64a164', '#3465a4', '#ad7fa8',]
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + CV-lim. (low-z)',
              r'  + Next-gen. spec-z',
              r'  + Stage 2 ($3\times$PB wedge)',
              r'  + CV-lim. (low-z) + Stage 2 ($3\times$PB wedge)',
              ]


if idx == 10:
    figname = "pub_bounds_mocker_compare_desi_%s.200.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    tmpl = "chains/final_mocker_seed200_cmb_lss%s"
    expts = [ '', 
              '_desi_hetdex',
              '_hirax_pbw',
              '_cosvis_pbw',
              ]
    
    colours = [ '#E1E1E1', '#555753', '#E6773D', '#3465a4']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + DESI + HETDEX',
              r'  + HIRAX ($3\times$PB wedge)',
              r'  + Stage 2 ($3\times$PB wedge)',
             ]


if idx == 11:
    figname = "pub_bounds_tanh_compare_spectel_%s.200.pdf" % MODE
    #if MODE == 'lowz': legend = False
    YLIM = (0.4, 1.1)
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ '', 
              '_cvlowz',
              '_cosvis_pbw',
              '_cvlowz_cosvis_pbw',
              '_spectel',
              '_spectel_desi',
              ]
             
    colours = [ '#E1E1E1', '#64a164', '#3465a4', '#ad7fa8', 'r', 'g']
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + CV-lim. (low-z)',
              r'  + Stage 2 ($3\times$PB wedge)',
              r'  + CV-lim. (low-z) + Stage 2 ($3\times$PB wedge)',
              r'  + SpecTel',
              r'  + SpecTel + DESI'
              ]

if idx == 12:
    figname = "pub_bounds_tanh_compare_cvallz.pdf"
    YLIM = (0.4, 1.1)
    
    dash = [False for j in range(40)]
    dash[1] = True
    
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss_cvallz-%s"
    #expts = ['%02d' % (i+1) for i in range(15)]
    expts = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    expts = [s.lower() for s in expts]
             
    colours = [ '#E1E1E1', '#64a164', '#3465a4', '#ad7fa8', 'r', 'g']
    colours += colours + colours + colours
    colours += colours
    fnames = [tmpl % e for e in expts]
    labels = expts


if idx == 13:
    figname = "pub_bounds_tanh_compare_highz_%s.200.pdf" % MODE
    #if MODE == 'lowz': legend = False
    YLIM = (0.2, 1.4)
    
    dash = [False for j in range(10)]
    dash[4] = True
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    expts = [ '', 
              '_cosvis_hw',
              '_cosvis_pbw',
              '_cosvis_nw',
              '_spectel',
              ]
             
    colours = [ '#E1E1E1', '#092f61', '#3465a4', '#84a9d9', '#11bdbd'] 
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + Stage 2 (horiz. wedge)',
              r'  + Stage 2 ($3\times$PB wedge)',
              r'  + Stage 2 (no wedge)',
              r'  + Next-gen. spec-z',
              ]


if idx == -1:
    P.title("Seed 200 Sel 22") # 21, 30
    #figname = "pub_bounds_tanh_hirax_%s.200.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
    #tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
    expts = [
        '_cosvis_hw', 
        '_cosvis_nw', 
        '_cosvis_pbw', 
        '_cvlowz_cosvis_pbw', 
        '_cvlowz', 
        '_cvlowz_hetdex', 
        '_cvlowz_hizrax_pbw', 
        '', 
        '_desi', 
        '_desi_hetdex', 
        '_desi_hirax_pbw', 
        '_desi_hizrax_hw', 
        '_desi_hizrax_nw', 
        '_desi_hizrax_pbw', 
        '_hetdex_hirax_pbw', 
        '_hirax_hw', 
        '_hirax_nw', 
        '_hirax_pbw', 
        '_hizrax_hw', 
        '_hizrax_nw', 
        '_hizrax_pbw',
    ]
             
    colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00', 'm', 'c', 'y',]
    colours += colours + colours + colours
    fnames = [tmpl % e for e in expts]
    labels = expts

# Load data and select random samples
for j, fn in enumerate(fnames):
    
    #if '_cosvis_pbw' in fn:
    #    fn = fn.replace("_seed21", "")
        
    
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
        print(dat['h'].shape)
        #print("Walkers:", dat['walker'][1200000:1200000+10])
        
        # Choose random subset of samples
        #sample_idxs = range(1200000, dat['h'].size, 1)
        np.random.seed(SAMPLE_SEED)
        sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
        print("Selecting %d samples out of %d available (burnin %d discarded)" \
              % (sample_idxs.size, dat['h'].size - BURNIN, BURNIN))
        
        # Calculate model for selected samples
        for i in sample_idxs:
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
        dashes = [1.5,1.5] if dash[j] else []
        P.plot(z, pct[0], color=colours[j], lw=LW, label=lbl, dashes=dashes)
        P.plot(z, pct[-1], color=colours[j], lw=LW, dashes=dashes)

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
    try:
        P.ylim(YLIM)
    except:
        P.ylim((0.2, 1.5))
    
    P.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    P.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
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
#P.savefig("pub_bounds_tanh_hetdex.200.pdf")
#P.savefig("pub_bounds_mocker.200.pdf")
try:
    P.savefig(figname)
    print("Figure saved as", figname)
except:
    pass
P.show()
