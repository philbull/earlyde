#!/usr/bin/env python
"""
Plot 2D scatter plot of multiple MCMC chains to see if they overlap.
"""
import numpy as np
import pylab as P
#import wztanh as model
from load_data_files import load_chain
import corner

nmin = 100000 #-100000
nmax = -1 #200000

fnames = [ #"chains/final_wztanh_cmb_lss", 
           #    "chains/final_wztanh_cmb_lss_cvlowz", 
           #    "chains/final_wztanh_seed16_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed17_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed18_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed19_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed20_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed21_cmb_lss_cvlowz",
               "chains/final_wztanh_cmb_lss_hirax_nw",
               "chains/final_wztanh_seed21_cmb_lss_hirax_nw",
             ]
"""
colours = [ #'#E1E1E1', 
            #'#555555', 
            '#863AB7', 
            '#E6773D', 
            '#F5BC00', 
            'g', 
            'm'
          ]
"""
colours=['r', 'b', 'm', 'c']
labels = [ #'CMB + LSS', 
           #'  + CV-limited low-z (15)',
           #'  + CV-limited low-z (16)',
           #'  + CV-limited low-z (17)',
           #'  + CV-limited low-z (18)',
           #'  + CV-limited low-z (19)',
           #'  + CV-limited low-z (20)',
           #'  + CV-limited low-z (21)',
           'HIRAX no wedge (10)',
           'HIRAX no wedge (21)',
         ]
ranges = {
    'omegaM':   (0.30, 0.33),
    'w0':       (-1.8, -0.8),
    'h':        (0.66, 0.685),
    'omegaB':   (0.048, 0.052),
    'winf':     (-2., 0.),
    'zc':       (0., 10.),
}


#fig = None
axes = [P.subplot(321), P.subplot(322), P.subplot(323), 
        P.subplot(324), P.subplot(325), P.subplot(326), ]
for i, fname in enumerate(fnames):
    
    # Load data
    dat = load_chain(fname, burnin=100000)
    print(dat.keys())
    print(dat['w0'].shape)
    
    lbls = [key for key in dat.keys()]
    lbls.remove('logl')
    lbls.remove('walker')
    data = np.array([dat[key][nmin:nmax] for key in lbls]).T
    
    for j, pname in enumerate(['omegaM', 'w0', 'h', 'omegaB', 'winf', 'zc',]):
        
        if i == 0:
            axes[j].set_xlabel(pname)
        
        if j == 0:
            axes[j].hist(data[:,j], bins=80, color=colours[i], histtype='step', 
                         label=labels[i], density=True, range=ranges[pname])
        else:
            axes[j].hist(data[:,j], bins=80, color=colours[i], histtype='step', 
                         density=True, range=ranges[pname])
        
        
    """
    if fig is None:
        fig = corner.corner(data, labels=lbls,
                            quantiles=[0.16, 0.5, 0.84], bins=80,
                            show_titles=True, color=colours[i],
                            title_kwargs={"fontsize": 12},
                            plot_density=False, plot_datapoints=False)
    else:
        corner.corner(data, labels=lbls,
                      quantiles=[0.16, 0.5, 0.84], bins=80,
                      show_titles=True, color=colours[i],
                      title_kwargs={"fontsize": 12}, fig=fig, 
                      plot_density=False, plot_datapoints=False)
    """
    
    # w0, winf
    #x = dat['w0'][nmin:nmax]
    #y = dat['h'][nmin:nmax]
    
    #P.plot(x, y, color=colours[i], ls='none', marker='.', alpha=0.15, label=labels[i])

axes[0].legend(loc='upper left')


#P.xlabel("$w_0$", fontsize=15)
#P.ylabel("$h$", fontsize=15)

#P.axvline(-1., color='k', ls='dashed')

#P.legend(loc='upper left', frameon=False)
P.tight_layout()

P.show()
