#!/usr/bin/env python
"""
Plot 2D scatter plot of multiple MCMC chains to see if they overlap.
"""
import numpy as np
import pylab as P
#import wztanh as model
from load_data_files import load_chain
import corner

#nmin = 0 #100000 #-100000
#nmax = -1 #200000

fnames = [ #"chains/final_wztanh_cmb_lss", 
           #    "chains/final_wztanh_cmb_lss_cvlowz", 
           #    "chains/final_wztanh_seed16_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed17_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed18_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed19_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed20_cmb_lss_cvlowz",
               #"chains/final_wztanh_seed21_cmb_lss_cvlowz",
               #"chains/final_wztanh_cmb_lss_cosvis_pbw",
               "chains/final_wztanh_seed21_cmb_lss_cosvis_pbw",
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
           #'CosVis PBW (10)',
           'CosVis PBW (21)',
         ]
ranges = {
    'omegaM':   (0.30, 0.33),
    'w0':       (-1.8, -0.8),
    'h':        (0.66, 0.685),
    'omegaB':   (0.048, 0.052),
    'winf':     (-2., 0.),
    'zc':       (0., 10.),
    'logl':     (-23., -4.),
}


"""
P.subplot(111)
for i, fname in enumerate(fnames):
    
    # Load data
    dat = load_chain(fname, burnin=100000)

    P.hist(dat['logl'], bins=80, color=colours[i], histtype='step', 
           label=labels[i], density=True, range=ranges['logl'])

P.tight_layout()
P.show()
exit()
"""


fig = None
#axes = [P.subplot(321), P.subplot(322), P.subplot(323), 
#        P.subplot(324), P.subplot(325), P.subplot(326), ]
for i, fname in enumerate(fnames):
    
    # Load data
    dat = load_chain(fname, burnin=1500000)
    print(dat.keys())
    print(dat['w0'].shape)
    
    #idxs = np.where(dat['logl'] > np.max(dat['logl']) - 5.)
    
    lbls = [key for key in dat.keys()]
    print(lbls)
    lbls.remove('logl')
    lbls.remove('walker')
    
    
    data = np.array([dat[key] for key in lbls]).T
    
    """
    for j, pname in enumerate(['omegaM', 'w0', 'h', 'omegaB', 'winf', 'zc',]):
        
        if i == 0:
            axes[j].set_xlabel(pname)
        
        if j == 0:
            axes[j].hist(dat[pname], bins=80, color=colours[i], 
                         histtype='step', label=labels[i], density=True, 
                         range=ranges[pname])
        else:
            axes[j].hist(dat[pname], bins=80, color=colours[i], 
                         histtype='step', density=True, 
                         range=ranges[pname])
    """ 
    
    if fig is None:
        fig = corner.corner(data, labels=lbls,
                            quantiles=[0.16, 0.5, 0.84], bins=80,
                            show_titles=True, color=colours[i],
                            title_kwargs={"fontsize": 12},
                            plot_density=False, plot_datapoints=False, 
                            levels=[0.68, 0.95])
    else:
        corner.corner(data, labels=lbls,
                      quantiles=[0.16, 0.5, 0.84], bins=80,
                      show_titles=True, color=colours[i],
                      title_kwargs={"fontsize": 12}, fig=fig, 
                      plot_density=False, plot_datapoints=False,
                      levels=[0.68, 0.95])
    

#axes[0].legend(loc='upper left')

P.tight_layout()

P.show()
