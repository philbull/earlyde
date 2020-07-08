#!/usr/bin/env python
"""
Plot 2D scatter plots of Delta w and w_0.
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


# Get parameter list
params = ['omegaM', 'w0', 'h', 'deltaw', 'zc', 'deltaz'] # 'omegaB'
pnames = [r'$\Omega_{\rm M}$', r'$w_0$', r'$h$', 
          r'$\Delta w$', r'$z_c$', r'$\Delta z$']
fiducial = [0.3166, -1., 0.6727, None, None, None]
# r'$\Omega_{\rm b}$' 0.04941

l1 = r'CMB + LSS'
l2 = r' + Stage 2 ($3\times$PB)'

deltaz1 = d1['deltaz'][-1000:].flatten()
deltaz2 = d2['deltaz'][-1000:].flatten()

# Plot scatter plot
plt.subplot(111)

plt.plot(d1['w0'][-1000:].flatten()[np.logical_and(deltaz1 > 1., deltaz1 < 2.)], 
         d1['deltaw'][-1000:].flatten()[np.logical_and(deltaz1 > 1., deltaz1 < 2.)], 
         'r.', alpha=0.3)

plt.plot(d2['w0'][-1000:].flatten()[np.logical_and(deltaz2 > 1., deltaz2 < 2.)], 
         d2['deltaw'][-1000:].flatten()[np.logical_and(deltaz2 > 1., deltaz2 < 2.)], 
         'b.', alpha=0.3)


plt.xlabel("$w_0$", fontsize=16)
plt.ylabel(r"$\Delta w \approx w_a$", fontsize=16)

plt.legend(loc='upper left', frameon=False, fontsize=13)
    
# Set tick labels and settings
plt.tick_params(axis='both', which='major', labelsize=16, size=7., 
                    width=1.5, pad=8., right=True, direction='in')
plt.tick_params(axis='both', which='minor', labelsize=18, size=4., 
                    width=1.5, pad=8., right=True, direction='in')

plt.gcf().set_size_inches((7., 5.))
#plt.savefig("pub_scatter_wparams.pdf")
plt.show()
