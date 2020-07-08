#!/usr/bin/env python
"""
Plot superposed 1D marginal histograms for two experiments.
"""
import numpy as np
import pylab as plt
from load_data_files import load_chain

fname1 = "chains/final_wztanh_seed200_cmb_lss"
fname2 = "chains/final_wztanh_seed200_cmb_lss_cosvis_pbw"

# Load data
d1 = load_chain(fname1)
d2 = load_chain(fname2)


# Reshape into per-walker chains
nw = np.max(d1['walker']).astype(int) + 1
for k in d1.keys():
    d1[k] = d1[k].reshape((-1, nw))
for k in d2.keys():
    d2[k] = d2[k].reshape((-1, nw))
#d2 = d2.reshape((d2.shape[0], -1, nw))
print(d1['walker'].shape, d2['walker'].shape)

# Create Figure
fig = plt.figure(constrained_layout=True)

gs = fig.add_gridspec(2, 3)
ax11 = fig.add_subplot(gs[0,0])
ax12 = fig.add_subplot(gs[0,1])
ax13 = fig.add_subplot(gs[0,2])
ax21 = fig.add_subplot(gs[1,0])
ax22 = fig.add_subplot(gs[1,1])
ax23 = fig.add_subplot(gs[1,2])
axes = [ax11, ax12, ax13, ax21, ax22, ax23]

# Get parameter list
params = ['omegaM', 'w0', 'h', 'deltaw', 'zc', 'deltaz'] # 'omegaB'
pnames = [r'$\Omega_{\rm M}$', r'$w_0$', r'$h$', 
          r'$\Delta w$', r'$z_c$', r'$\Delta z$']
fiducial = [0.3166, -1., 0.6727, 0., None, None]
# r'$\Omega_{\rm b}$' 0.04941

l1 = r'CMB + LSS'
l2 = r' + Stage 2 ($3\times$PB)'

# Loop over parameters
for j, p in enumerate(params):
    
    # FIXME: Due to a bug, 'deltaw' is actually '-deltaw' in runs prior to 
    # 2020-06-15.
    sgn = 1.
    if p == 'deltaw': sgn = -1.
    
    # Get percentiles for w0
    if p == 'w0':
        pct = np.percentile(sgn*d2[p][-1000:].flatten(), [2.5, 50., 97.5])
        print("w0 percentiles (2.5%, 97.5%): ", pct)
        print("w0 (95% CL):", np.median(sgn*d2[p][-1000:].flatten()), pct[-1] - pct[0])
    
    # Get min/max for both experiments
    vmin = np.min([d1[p][-1000:].min(), d2[p][-1000:].min()])
    vmax = np.max([d1[p][-1000:].max(), d2[p][-1000:].max()])
    
    # Plot histograms
    axes[j].hist(sgn*d1[p][-1000:].flatten(), bins=70, density=True, alpha=0.4, 
                 range=(vmin, vmax), color='gray', label=l1)
    axes[j].hist(sgn*d2[p][-1000:].flatten(), bins=70, density=True, alpha=0.4, 
                 range=(vmin, vmax), color='#3465a4', label=l2)
    
    # Plot fiducial values
    if fiducial[j] is not None:
        axes[j].axvline(fiducial[j], color='k', ls='dashed', alpha=0.7, lw=1.5)
    
    axes[j].set_xlabel(pnames[j], fontsize=16)
    axes[j].set_xlim((vmin, vmax))
    #axes[j].set_yscale('log')
    
    # Add legend
    if j == 1:
        axes[j].legend(loc='upper left', frameon=False, fontsize=13)
    
    # Set tick labels and settings
    axes[j].tick_params(axis='both', which='major', labelsize=16, size=7., 
                        width=1.5, pad=8., right=True, direction='in')
    axes[j].tick_params(axis='both', which='minor', labelsize=18, size=4., 
                        width=1.5, pad=8., right=True, direction='in')
    axes[j].tick_params(axis='y', labelleft=False)
    
plt.gcf().set_size_inches((16., 7.))
plt.savefig("pub_marginals_stage2.pdf")
plt.show()
