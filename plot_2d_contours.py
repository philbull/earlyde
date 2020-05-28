#!/usr/bin/env python
"""
Plot 2D contours for a pair of parameters.
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
    print("Requires one arguments: idx")
    sys.exit(1)
idx = int(sys.argv[1])

SAMPLE_SEED = 88 #44
BURNIN = 200000 #1200000
NSAMP = 150000 #50000

legend = False

plt.title("Seed 100") # 21, 30
#figname = "pub_bounds_tanh_hirax_%s.pdf" % MODE
#if MODE == 'lowz': legend = False

tmpl = "chains/seed100_selseed22/final_wztanh_seed100_cmb_lss%s"
#tmpl = "chains/final_wztanh_seed100_cmb_lss%s"
expts = [
    '_cosvis_hw',           # 0
    '_cosvis_nw',           # 1
    '_cosvis_pbw',          # 2
    '_cvlowz_cosvis_pbw',   # 3
    '_cvlowz',              # 4
    '_cvlowz_hetdex',       # 5
    '_cvlowz_hizrax_pbw',   # 6
    '',                     # 7
    '_desi',                # 8
    '_desi_hetdex',         # 9
    '_desi_hirax_pbw',      # 10
    '_desi_hizrax_hw',      # 11
    '_desi_hizrax_nw',      # 12
    '_desi_hizrax_pbw',     # 13
    '_hetdex_hirax_pbw',    # 14
    '_hirax_hw',            # 15
    '_hirax_nw',            # 16
    '_hirax_pbw',           # 17
    '_hizrax_hw',           # 18
    '_hizrax_nw',           # 19
    '_hizrax_pbw',          # 20
]
         
colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00', 'm', 'c', 'y',]
colours += colours + colours + colours
fnames = [tmpl % e for e in expts]
labels = expts

# Create figure
fig = None

# Load data and select random samples
for j, fn in enumerate(fnames[:1]):
    
    # Load samples
    print("%s: Loading samples" % fn)
    dat = load_chain(fn)
    
    # Choose random subset of samples
    np.random.seed(SAMPLE_SEED)
    sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
    print("Selecting %d samples out of %d available (burnin %d discarded)" \
          % (sample_idxs.size, dat['h'].size - BURNIN, BURNIN))
    
    #p1, p2 = PARAMS
    #d = np.vstack( [dat[p1][sample_idxs], dat[p2][sample_idxs]] ).T
    params = [_p for _p in dat.keys()]
    d = np.vstack( [dat[_p][sample_idxs] for _p in dat.keys()] ).T
    corner.corner(d, labels=params, fig=fig, plot_datapoints=False, plot_contours=False)
    if fig is None: fig = plt.gcf()

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

P.xlabel("$z$", fontsize=15)
P.ylabel(r"$w(z)$", fontsize=15, labelpad=10)

P.gca().tick_params(axis='both', which='major', labelsize=18, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=18, size=5., 
                    width=1.5, pad=8.)

"""
# Re-order legend
#handles, labels = plt.gca().get_legend_handles_labels()
#handles = [handles[-1], handles[0], handles[1], handles[2], handles[3]]
#labels = [labels[-1], labels[0], labels[1], labels[2], labels[3]]
#leg = plt.legend(handles, labels, loc='upper left', frameon=False)

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
    #plt.savefig(figname)
    plt.savefig("corner.pdf")
    #print("Figure **NOT** saved as", figname)
except:
    pass
plt.show()
