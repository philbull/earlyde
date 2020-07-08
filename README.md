# earlyde
Fit early-ish dark energy models to background distance/expansion rate data.
This code was used to produce the plots in 
[arXiv:2007.02865](https://arxiv.org/abs/2007.02865).

Requirements
============
 * numpy, scipy, matplotlib
 * pyccl (for some test functions)
 * emcee

Usage
=====
This code is meant to be used with Fisher matrices output by 
[RadioFisher](https://gitlab.com/radio-fisher/bao21cm). Some example 
matrices are included in this repository. These matrices model the covariances 
for the angular diameter distance and expansion rate for a set of redshift 
bins. Nuisance parameters, such as bias and growth rate, may also be included.

You can, however, modify the code to use whatever likelihood function you 
want. Simplified Planck data (i.e. using shift parameters; see 
[arXiv:2007.02865](https://arxiv.org/abs/2007.02865)) are also included here 
for example.

To run this code, use the `fit_wztanh.py` script. This runs an `emcee` MCMC 
chain using particular experimental settings. It's best to inspect the code 
to see what settings can be changed. The code does run on multiple cores.

The dark energy models themselves are defined in `wztanh.py`.

A variety of plotting scripts are also provided. The main one is 
`plot_bounds_all.py`. The `plot_obs_errors.py` script is also useful for 
inspecting the constraints implied by the various Fisher matrices.

Notes
=====
 - For all of the results generated before 2020-06-15, there was a minus sign 
   bug in the definition of the Tracker model 'deltaw' parameter (if used). It 
   works out so that 'deltaw' in these runs is actually '-deltaw'. This bug has 
   now been fixed (but double-check the code).
