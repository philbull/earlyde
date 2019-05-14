#!/usr/bin/env python
"""
Plot errors on observables that are used in the simulated likelihoods.
"""
import numpy as np
import pylab as P
import wztanh as model
from load_data_files import *
from matplotlib import ticker

col_hrx = '#D92E1C'; col_hrx_light = '#FF9B99'
#col_hizrx = '#E3B4B0'; col_hizrx_light = '#E8C2BF'
#col_hizrx = '#81A871'; col_hizrx_light = '#81A871'
col_hizrx = col_hizrax_light = '#E9AD2B'
col_cv = '#00529B'; col_cv_light = '#85C5FF'
col_cvlowz = '#3E893E'

ii = 0 # H(z)
#ii = 1 # D_A(z)

# CMB: Planck 2015, Gaussianised
pl15_mean, pl15_icov = load_planck_data("planck_derived_fisher_distances", 
                                      params=['omegabh2', 'omegamh2', 'DAstar'])

# DESI, Stage 2, and HIRAX Fisher matrices, (H, D_A)
fisher_desi = ["Fisher-full-gDESI_CV_desicv",]
fisher_cosvis = [ "Fisher-full-iCosVis256x256_2yr_nowedge",
                  "Fisher-full-iCosVis256x256_2yr_3pbwedge", 
                  "Fisher-full-iCosVis256x256_2yr_wedge" ]
fisher_hirax = [ "Fisher-full-iHIRAX_2yr",
                 "Fisher-full-iHIRAX_2yr_3pbwedge",
                 "Fisher-full-iHIRAX_2yr_horizwedge" ]
fisher_hizrax = [ "Fisher-full-iHIRAX_highz_2yr",
                  "Fisher-full-iHIRAX_highz_2yr_3pbwedge",
                  "Fisher-full-iHIRAX_highz_2yr_horizwedge" ]
fisher_hetdex = ["Fisher-full-HETDEXdz03",]
fisher_cvlowz = ["Fisher-full-gCVLOWZ"]

fishers = fisher_desi + fisher_cosvis + fisher_hirax + fisher_hizrax \
        + fisher_hetdex + fisher_cvlowz
names = ['DESI', 'CosVis_nw', 'CosVis_3pb', 'CosVis_hw', 
         'HIRAX_nw', 'HIRAX_3pb', 'HIRAX_hw', 
         'HIZRAX_nw', 'HIZRAX_3pb', 'HIZRAX_hw',
         'HETDEXdz03', 'CVLOWZ']

# Get std deviations for Fisher matrices
zc = {}; frac_std = {}
for i, fname in enumerate(fishers):
    expt_zc, expt_mean, expt_icov = load_fisher_data(fname)
    expt_cov = np.linalg.inv(expt_icov)
    expt_std = np.sqrt( np.diag(expt_cov) )
    
    zc[names[i]] = expt_zc
    print(expt_zc)
    frac_std[names[i]] = (expt_std[:expt_zc.size] / expt_mean[:expt_zc.size], 
                          expt_std[expt_zc.size:] / expt_mean[expt_zc.size:])
                          # sigma_H/H, sigma_DA/DA


# Plot errors
P.subplot(111)

# CVLOWZ
_zc = zc['CVLOWZ']
mu = frac_std['CVLOWZ'][ii]
P.plot(_zc, mu, marker='o', color=col_cvlowz, lw=1.8, alpha=0.8)

# HIRAX high-z
_zc = zc['HIZRAX_nw']
mu = frac_std['HIZRAX_3pb'][ii]
dp = frac_std['HIZRAX_hw'][ii] - mu
dm = mu - frac_std['HIZRAX_nw'][ii]

P.fill_between(_zc, mu - dm, mu + dp, color=col_hizrx, alpha=0.2)
#P.errorbar(_zc, mu, yerr=(dm, dp), 
#           marker='o', color=col_hizrx, ls='none', ms=5., 
#           elinewidth=1.5, ecolor=col_hizrx_light)
P.plot(_zc, mu, marker='o', color=col_hizrx, ms=5.)
P.plot(_zc, mu - dm, color=col_hizrx, marker='o', ms=3.)
P.plot(_zc, mu + dp, color=col_hizrx, marker='o', ms=3.)


# HIRAX
_zc = zc['HIRAX_nw']
mu = frac_std['HIRAX_3pb'][ii]
dp = frac_std['HIRAX_hw'][ii] - mu
dm = mu - frac_std['HIRAX_nw'][ii]

P.fill_between(_zc, mu - dm, mu + dp, color=col_hrx, alpha=0.2)
P.errorbar(_zc, mu, yerr=(dm, dp), 
           marker='o', color=col_hrx, ls='none', ms=5., 
           elinewidth=1.5, ecolor='none')
P.plot(_zc, mu - dm, color=col_hrx, marker='o', ms=3.)
P.plot(_zc, mu + dp, color=col_hrx, marker='o', ms=3.)

# CosVis
_zc = zc['CosVis_nw']
mu = frac_std['CosVis_3pb'][ii]
dp = frac_std['CosVis_hw'][ii] - mu
dm = mu - frac_std['CosVis_nw'][ii]

P.fill_between(_zc, mu - dm, mu + dp, color=col_cv, alpha=0.2)
P.errorbar(_zc, mu, yerr=(dm, dp), 
           marker='o', color=col_cv, ls='none', ms=5., 
           elinewidth=0.0, ecolor='none') #col_cv_light)
P.plot(_zc, mu - dm, color=col_cv, marker='o', ms=3.)
P.plot(_zc, mu + dp, color=col_cv, marker='o', ms=3.)


# DESI
P.plot(zc['DESI'], frac_std['DESI'][ii], 'ko-', label="DESI")

# HETDEX
P.plot(zc['HETDEXdz03'], frac_std['HETDEXdz03'][ii], marker='o', 
       color="#863AB7",  label="HETDEX")

P.plot(-_zc, mu, marker='o', color=col_cvlowz, lw=1.8, alpha=0.8,
       label="CV-lim low-$z$") # for legend

P.plot(-_zc, mu + dp, color=col_hrx, marker='o', ms=5., label="HIRAX") # for legend
P.plot(-_zc, mu + dp, color=col_hizrx, marker='o', ms=5., 
       label="HIRAX high-z") # for legend
P.plot(-_zc, mu + dp, color=col_cv, marker='o', ms=5., label="Stage 2") # for legend

P.xlim((0., 6.5))
P.yscale('log')
P.xlabel("$z$", fontsize=15)
if ii == 0:
    #P.ylim((5e-4, 0.03))
    #P.ylim((5e-4, 0.2))
    P.ylim((5e-4, 7.))
    P.ylabel(r"$\sigma_H / H$", fontsize=15, labelpad=10)
    P.legend(loc='upper left', frameon=False, prop={'size': 'x-large'}, ncol=2)
else:
    #P.ylim((5e-4, 0.8))
    P.ylim((5e-4, 7.))
    P.ylabel(r"$\sigma_{D_A} / D_A$", fontsize=15, labelpad=10)
    #P.legend(loc='upper left', frameon=False, prop={'size': 'large'}, ncol=2)

P.gca().xaxis.set_major_locator(ticker.MultipleLocator(2.0))
P.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

P.gca().tick_params(axis='both', which='major', labelsize=16, size=8., 
                    width=1.5, pad=8.)
P.gca().tick_params(axis='both', which='minor', labelsize=16, size=5., 
                    width=1.5, pad=8.)

P.gcf().set_size_inches((7., 5.))

P.tight_layout()
if ii == 0:
    P.savefig("pub_error_H.pdf")
else:
    P.savefig("pub_error_DA.pdf")
P.show()
