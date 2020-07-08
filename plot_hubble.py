#!/usr/bin/env python
"""
Plot 1D histogram for a parameter.
"""
import numpy as np
import pylab as plt
import wztanh as model
from load_data_files import load_chain
from matplotlib import ticker
import os 
import corner

SAMPLE_SEED = 22 #44
BURNIN = 80 * 500 # Leaves 80 * 1500 samples #200000 #1200000
NSAMP = 80 * 1500 #150000 #50000

PARAM = 'h'

legend = True

figname = "pub_h0_tanh_all.200.pdf"

tmpl = "chains/final_wztanh_seed200_cmb_lss%s"
expts = [
    '',
    
    '_desi', 
    '_desi_hetdex', 
    '_desi_hirax_pbw', 
    
    '_hizrax_hw', 
    '_hizrax_pbw',
    '_hizrax_nw', 
    
    '_hirax_hw', 
    '_hirax_pbw', 
    '_hirax_nw', 
    
    '_cosvis_hw', 
    '_cosvis_pbw', 
    '_cosvis_nw', 
    
    '_cvlowz', 
    '_cvlowz_hetdex', 
    '_cvlowz_hizrax_pbw', 
    '_cvlowz_cosvis_pbw', 
    
    #'_desi_hizrax_hw', 
    #'_desi_hizrax_nw', 
    #'_desi_hizrax_pbw', 
    #'_hetdex_hirax_pbw', 
]

labels = [
    'CMB + LSS', 
    '  + DESI',
    '  + DESI + HETDEX', 
   r'  + DESI + HIRAX ($3\times$ PB)', 
    
    '  + HIRAX high-z (horiz. wedge)', 
   r'  + HIRAX high-z ($3\times$ PB)',
    '  + HIRAX high-z (no wedge)',
     
    '  + HIRAX (horiz. wedge)', 
   r'  + HIRAX ($3\times$ PB)',
    '  + HIRAX (no wedge)', 
    
    '  + Stage 2 (horiz. wedge)', 
   r'  + Stage 2 ($3\times$ PB)',
    '  + Stage 2 (no wedge)',
    
    '  + CV-lim. low-z', 
    '  + CV-lim. low-z + HETDEX', 
   r'  + CV-lim. low-z + HIRAX high-z ($3\times$ PB)', 
   r'  + CV-lim. low-z + Stage 2 ($3\times$ PB)', 

]

colours = ['k',  
           'gray', 'gray', 'gray',   
           '#D92E1C', '#D92E1C', '#D92E1C', 
           '#E6773D', '#E6773D', '#E6773D', 
           '#3465a4', '#3465a4', '#3465a4',
           '#64a164', '#64a164', '#64a164', '#64a164',]

#colours = [ '#E1E1E1', '#D92E1C', '#E6773D', '#F5BC00', 'm', 'c', 'y',]
#colours += colours + colours + colours
fnames = [tmpl % e for e in expts]

# Create figure
ax = plt.subplot(111)

# Load data and select random samples
std_vals = []
for j, fn in enumerate(fnames):
    
    # Load samples
    print("%s: Loading samples" % fn)
    try:
        dat = load_chain(fn)
    except:
        continue
    
    # Choose random subset of samples
    np.random.seed(SAMPLE_SEED)
    sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
    print("Selecting %d samples out of %d available (burnin %d discarded)" \
          % (sample_idxs.size, dat['h'].size - BURNIN, BURNIN))
    
    std = 100. * np.std(dat[PARAM][sample_idxs])
    print(labels[j])
    plt.plot(j, std, marker='o', color=colours[j], label=labels[j])

plt.xlim((-0.8, 16.6))
plt.ylim((0., 2.9))
#plt.xlabel("$%s$", fontsize=15)
plt.ylabel(r"$\sigma(H_0)$ [km/s/Mpc]", fontsize=15, labelpad=10)

plt.gca().tick_params(axis='both', which='major', labelsize=16, size=7., 
                      width=1.5, pad=8., right=True, direction='in')
plt.gca().tick_params(axis='both', which='minor', labelsize=18, size=4., 
                      width=1.5, pad=8., right=True, direction='in')

plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(1.0))
plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
#plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

plt.xticks(ticks=[])
#plt.xticks(ticks=np.arange(len(expts), dtype=int), labels=labels, rotation=45, size=10)

plt.fill_betweenx([0., 3.], [-1.,-1.], [0.5, 0.5], color='gray', alpha=0.1, zorder=50, lw=0)
plt.fill_betweenx([0., 3.], [3.5,3.5], [6.5, 6.5], color='gray', alpha=0.1, zorder=50, lw=0)
plt.fill_betweenx([0., 3.], [9.5,9.5], [12.5, 12.5], color='gray', alpha=0.1, zorder=50, lw=0)

# Labels
plt.text(0.0, 1.9, r"  CMB + LSS",   color='k', fontsize=11, rotation=90, ha='center')
plt.text(2.0, 1.9, r"  + DESI",     color='gray', fontsize=11, rotation=90, ha='center')
plt.text(5.0, 1.9, r"  + HIRAX high-z", color='#D92E1C', fontsize=11, rotation=90, ha='center')
plt.text(8.0, 1.9, r"  + HIRAX",    color='#E6773D', fontsize=11, rotation=90, ha='center')
plt.text(11., 1.9, r"  + Stage 2", color='#3465a4', fontsize=11, rotation=90, ha='center')
plt.text(14.5, 1.9, r"  + CV-lim. low-z", color='#64a164', fontsize=11, rotation=90, ha='center')

#plt.legend(loc='upper left')
plt.tight_layout()
try:
    plt.savefig(figname)
    print("Figure saved as", figname)
except:
    pass
plt.show()
