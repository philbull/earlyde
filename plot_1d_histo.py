#!/usr/bin/env python
"""
Plot 1D histogram for a parameter.
"""
import numpy as np
import pylab as plt
import wztanh as model
from load_data_files import load_chain
from matplotlib import ticker
import sys, os 
import corner

# Get arguments
if len(sys.argv) != 2:
    print("Requires one argument: idx")
    sys.exit(1)
idx = int(sys.argv[1])

SAMPLE_SEED = 88 #44
BURNIN = 200000 #1200000
NSAMP = 150000 #50000

PARAM = 'h'

legend = True

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
    tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
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
    figname = "pub_histo_tanh_hirax_%s.pdf" % MODE
    legend = True
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
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
    #figname = "pub_histo_tanh_hirax_%s.pdf" % MODE
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
    figname = "pub_histo_tanh_compare_desi_%s.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
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
    figname = "pub_histo_tanh_compare_galsurv_%s.pdf" % MODE
    #if MODE == 'lowz': legend = False
    
    dash = [False for j in range(10)]
    dash[1] = True
    LW = 2.1
    
    tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
    expts = [ '', 
              #'_desi',
              '_cvlowz',
              '_cvlowz_hetdex',
              '_cvlowz_hizrax_pbw',
              '_cvlowz_cosvis_pbw',
              ]
             
    colours = [ '#E1E1E1', '#000000', '#555753', '#fcaf3e', '#ad7fa8', '#3465a4',]
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS',
              #r'  + DESI', 
              r'  + CV-lim. (low-z)',
              r'  + CV-lim. (low-z) + HETDEX',
              r'  + CV-lim. (low-z) + HIRAX high-z ($3\times$PB)',
              r'  + CV-lim. (low-z) + Stage 2 ($3\times$PB)',
              ]


if idx == 9:
    figname = "pub_histo_tanh_compare_lowzhighz_%s.pdf" % MODE
    #if MODE == 'lowz': legend = False
    YLIM = (0.4, 1.1)
    
    dash = [False for j in range(10)]
    dash[1] = True
    
    
    tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
    expts = [ '', 
              '_cvlowz',
              '_cosvis_pbw',
              '_cvlowz_cosvis_pbw',
              ]
             
    colours = [ '#E1E1E1', '#000000', '#3465a4', '#ad7fa8',]
    fnames = [tmpl % e for e in expts]
    labels = [ 'CMB + LSS', 
              r'  + CV-lim. (low-z)',
              r'  + Stage 2 ($3\times$PB wedge)',
              r'  + CV-lim. (low-z) + Stage 2 ($3\times$PB wedge)',
              ]


if idx == -1:
    plt.title("Seed 100 Sel 22") # 21, 30
    #figname = "pub_histo_tanh_hirax_%s.pdf" % MODE
    if MODE == 'lowz': legend = False
    
    tmpl = "chains/seed100_selseed22/final_wztanh_seed100_cmb_lss%s"
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

# Create figure
ax = plt.subplot(111)

# Load data and select random samples
for j, fn in enumerate(fnames[:3]):
    
    # Load samples
    print("%s: Loading samples" % fn)
    dat = load_chain(fn)
    
    # Choose random subset of samples
    np.random.seed(SAMPLE_SEED)
    sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
    print("Selecting %d samples out of %d available (burnin %d discarded)" \
          % (sample_idxs.size, dat['h'].size - BURNIN, BURNIN))
    
    ax.hist(dat[PARAM][sample_idxs], bins=160, range=(0.62, 0.73), 
            density=True, histtype='step', color=colours[j],
            label=labels[j])

"""
if MODE == 'lowz':
    plt.xlim((0., 2.2))
    plt.ylim((-1.5, -0.2))
    
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
else:
    plt.xlim((0., 1170.))
    try:
        plt.ylim(YLIM)
    except:
        plt.ylim((-1.5, -0.2))
    
    plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    #plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    #plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
"""

plt.xlabel("$%s$" % PARAM, fontsize=15)
#plt.ylabel(r"$w(z)$", fontsize=15, labelpad=10)

plt.gca().tick_params(axis='both', which='major', labelsize=18, size=8., 
                      width=1.5, pad=8.)
plt.gca().tick_params(axis='both', which='minor', labelsize=18, size=5., 
                      width=1.5, pad=8.)

plt.gca().tick_params(axis='y', which='major', labelsize=18, size=5., 
                      width=1.5, pad=8., label1On=False)

# Re-order legend
#handles, labels = plt.gca().get_legend_handles_labels()
#handles = [handles[-1], handles[0], handles[1], handles[2], handles[3]]
#labels = [labels[-1], labels[0], labels[1], labels[2], labels[3]]
#leg = plt.legend(handles, labels, loc='upper left', frameon=False)

plt.legend(loc='upper left', frameon=False) 
"""
if legend:
    # Re-order legend so CMB+LSS is always first
    _handles, _labels = plt.gca().get_legend_handles_labels()
    handles = [_handles[-1],]
    labels = [_labels[-1],]
    for i in range(0, len(_handles)-1):
        handles += [_handles[i],]
        labels += [_labels[i],]
    leg = plt.legend(handles, labels, loc='upper left', frameon=False) 
"""
plt.tight_layout()
try:
    plt.savefig(figname)
    print("Figure saved as", figname)
except:
    pass
plt.show()
