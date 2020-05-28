#!/usr/bin/env python
"""
Run MCMC on tanh model.
"""
import numpy as np
import pylab as P
import copy, time, sys
import emcee
import wztanh as model
from load_data_files import *

np.random.seed(200)

# MCMC sampler settings
NTHREADS = 28
NSAMPLES = 150000
NWALKERS = 80
MOCKER = False

CHAIN_FILE_ROOT = "chains/final_wztanh_seed200"

# Which likelihoods to use
use_lss = True
use_planck2015 = False
use_planck = True

# Choose which Fisher-based likelihoods to use
expt_sets = [
    [],                         # 0
    ['DESI',],                  # 1
    ['DESI', 'HETDEX'],         # 2
    ['HIRAX_nw'],               # 3
    ['HIRAX_pbw'],              # 4
    ['HIRAX_hw'],               # 5
    ['DESI', 'HIRAX_pbw'],      # 6
    ['HETDEX', 'HIRAX_pbw'],    # 7
    ['CosVis_nw',],             # 8
    ['CosVis_pbw',],            # 9
    ['CosVis_hw',],             # 10
    ['HIZRAX_nw',],             # 11
    ['HIZRAX_pbw',],            # 12
    ['HIZRAX_hw',],             # 13
    ['DESI', 'HIZRAX_nw'],      # 14
    ['DESI', 'HIZRAX_pbw'],     # 15
    ['DESI', 'HIZRAX_hw'],      # 16
    ['CVLOWZ',],                # 17
    ['CVLOWZ', 'HETDEX'],       # 18
    ['CVLOWZ', 'HIZRAX_pbw',],  # 19
    ['CVLOWZ', 'CosVis_pbw']    # 20
]
use_expts = [] #['DESI', 'HETDEX']


# Available Fisher matrix likelihoods
fisher_mats = {
    'DESI':         "Fisher-full-gDESI_CV_desicv",
    'HETDEX':       "Fisher-full-HETDEXdz03",
    'CosVis_nw':    "Fisher-full-iCosVis256x256_2yr_nowedge",
    'CosVis_pbw':   "Fisher-full-iCosVis256x256_2yr_3pbwedge", 
    'CosVis_hw':    "Fisher-full-iCosVis256x256_2yr_wedge", 
    'HIRAX_nw':     "Fisher-full-iHIRAX_2yr",
    'HIRAX_pbw':    "Fisher-full-iHIRAX_2yr_3pbwedge",
    'HIRAX_hw':     "Fisher-full-iHIRAX_2yr_horizwedge",
    'HIZRAX_nw':    "Fisher-full-iHIRAX_highz_2yr",
    'HIZRAX_pbw':   "Fisher-full-iHIRAX_highz_2yr_3pbwedge",
    'HIZRAX_hw':    "Fisher-full-iHIRAX_highz_2yr_horizwedge",
    'CVLOWZ':       "Fisher-full-gCVLOWZ"
}

# Sound horizon at drag epoch
# Planck 2015 LCDM Mnu = 0.06 eV, Neff = 3.046, see p6 of 1411.1074
#r_d = 147.33 #147.49 # \pm 0.59 Mpc
# Planck 2018, arXiv:1807.06209 (Table 2; TT,TE,EE+lowE+lensing)
#r_d = 147.09 # +/- 0.26

# Planck 2018, arXiv:1807.06209 (Table 2; TT,TE,EE+lowE)
r_d = 147.05 # +/- 0.30

# Planck distance scale, from Table 2 of arXiv:1807.06209 (TT,TE,EE+lowE+lensing)
#rstar, sigma_rstar = 144.43, 0.26 # Mpc
#thstar, sigma_thstar = 1.04110e-2, 0.00031e-2 # radians
#zstar, sigma_zstar = 1089.92, 0.25

# Planck distance scale, from Table 2 of arXiv:1807.06209 (TT,TE,EE+lowE)
rstar, sigma_rstar = 144.39, 0.30 # Mpc
thstar, sigma_thstar = 1.04109e-2, 0.00030e-2 # radians
zstar, sigma_zstar = 1089.95, 0.27
#planck_dm = rstar / thstar / (1. + zstar) # Mpc
#planck_sigma_dm = planck_dm * np.sqrt(  (sigma_rstar/rstar)**2. \
#                                      + (sigma_thstar/thstar)**2. \
#                                      + (sigma_zstar/(1.+zstar)) )

# Low-z distance data
# - Combined BAO data from arXiv:1411.1074 (Table II)
# - Combined Lya-Lya and Lya-QSO constraints from arXiv:1904.03430 (Eqs. 50-51)
lss_data = [
    ('6dFGS', 'DV', 0.106, 3.047, 0.137),
    ('MGS', 'DV', 0.15, 4.480, 0.168),
    ('BOSS LOWZ', 'DV', 0.32, 8.467, 0.167),
    ('LyaF auto+QSO', ('DM', 'DH'), 2.34, (37.0, 9.00), (1.3, 0.22, -0.40)),
    ('BOSS CMASS comb', ('DM', 'DH'), 0.57, (14.945, 20.75), (0.210, 0.73, -0.52)),
    #('BOSS CMASS', 'DM', 0.57, 14.945, 0.210),
    #('BOSS CMASS', 'DH', 0.57, 20.75, 0.73),
    #('LyaF auto', 'DM', 2.34, 37.675, 2.171),
    #('LyaF auto', 'DH', 2.34, 9.18, 0.28),
    #('LyaF-QSO', 'DM', 2.36, 36.288, 1.344),
    #('LyaF-QSO', 'DH', 2.36, 9.00, 0.30),
    ('CMB approx', 'DM', zstar, 94.51, np.sqrt(0.004264)), # values ignored
    # FIXME: Keep CMB line for z value only
]


def fisher_like(expt_name, zc, dh, dm, verbose=False):
    """
    Calculate Fisher matrix-based likelihood, assuming H(z) and D_A(z) are the 
    parameters.
    """
    # Get pre-loaded Fisher data for this experiment
    expt_zc, expt_mean, expt_icov = expt_data[expt_name]
    
    # Convert distances in this model into H and D_A
    # Units: expects H ~ 100 km/s/Mpc, D_A ~ Gpc
    Hz = [model.C / dh[np.where(zc==_z)][0] / 1e2 for _z in expt_zc]
    DAz = [dm[np.where(zc==_z)][0] / (1.+_z) / 1e3 for _z in expt_zc]
    mean = np.concatenate((Hz, DAz))
    x = mean - expt_mean
    
    logL = -0.5 * np.dot(x, np.dot(expt_icov, x))
    if verbose: print("\t%16s: %3.3f" % ("%s Fisher" % expt_name, logL))
    return logL


def loglike(pvals, pnames, params0, priors, verbose=False):
    """
    Evaluate total log-likelihood for the set of input parameter values.
    """
    # Build parameter dictionary from input list
    p = copy.copy(params0)
    for i in range(len(pnames)): p[pnames[i]] = pvals[i]
    
    # Apply prior ranges if specified
    for pn in priors.keys():
        if pn not in p.keys(): continue # Skip unidentified parameters
        if pn not in pnames: continue # Only worry about sampled params
        pmin, pmax = priors[pn]
        if p[pn] < pmin: return -np.inf
        if p[pn] > pmax: return -np.inf
    
    # Enforce non-flat w(z) FIXME
    #if np.abs(p['w0'] - p['winf']) < 0.1: return -np.inf
    
    # Collect all available datapoints
    dname, dtype, zc, dval, derr = zip(*lss_data)
    
    # Make sure external experiments have distances calculated
    for e in expt_data.keys():
        zc = np.concatenate((zc, expt_data[e][0]))
    
    # Calculate model values for these datapoints
    dv, dm, dh = model.lss_distances(np.array(zc), p)
    model_calc = {'DV': dv, 'DM': dm, 'DH': dh}
    
    # Calculate simple independent Gaussian log-likelihoods for LSS
    logL = 0
    for i in range(len(dname)):
        if 'CMB' in dname[i]: continue # Skip CMB for now
        
        # Pair of values with covariance
        if isinstance(dtype[i], tuple):
            assert len(dtype[i]) == 2, "Only supports 2x2 covariances for now"
            dval_vec = np.array(dval[i])
            dist_vec = np.array([model_calc[dt][i] for dt in dtype[i]]) / r_d
            
            # Build covariance matrix
            std1, std2, rho_corr = derr[i]
            covmat = np.zeros((2,2))
            covmat[0,0] = std1**2.
            covmat[1,1] = std2**2.
            covmat[0,1] = covmat[1,0] = std1 * std2 * rho_corr
            
            # Calculate likelihood
            x = dval_vec - dist_vec
            _logL = -0.5 * np.dot(x.T, np.dot(np.linalg.inv(covmat), x))
            
        else:
            # Single value likelihood
            dist = model_calc[dtype[i]] # Get distance measure for this data point
            _logL = -0.5*(dval[i] - dist[i]/r_d)**2. / derr[i]**2.
            
        if verbose: print("\t%16s: %3.3f" % (dname[i], _logL))
        if use_lss: logL += _logL
    
    # Calculate log likelihoods for Fisher-based experiments
    for ename in use_expts:
        logL += fisher_like(ename, zc, dh, dm, verbose=verbose)
    
    # Find where CMB item is in the list of calculated distances
    cmb_idxs = [i for i in range(len(dname)) if 'CMB' in dname[i]]
    if len(cmb_idxs) != 1: raise KeyError("Should be only 1 CMB item in lss_data.")
    idx = cmb_idxs[0]
    assert np.abs(zc[idx] - 1090.) < 10., "CMB datapoint found at wrong redshift!"
    
    # Planck 2015 Gaussianised likelihood
    # Assumed order is 'omegabh2', 'omegamh2', 'DAstar'
    if use_planck2015:
        h = p['h']
        model_vec = np.array([p['omegaB']*h**2., p['omegaM']*h**2., dm[idx]])
        x = model_vec - pl15_mean
        _logL = -0.5 * np.dot(x, np.dot(pl15_icov, x).T)
        if verbose: print("\t%16s: %3.3f" % ("Planck 2015 CMB", _logL))
        if use_planck2015: logL += _logL
    
    # Planck 2018 Gaussianised likelihood
    # Assumed order is 'omegabh2', 'omegamh2', 'DAstar'
    h = p['h']
    model_vec = np.array([p['omegaB']*h**2., p['omegaM']*h**2., dm[idx]])
    x = model_vec - pl18_mean
    _logL = -0.5 * np.dot(x, np.dot(pl18_icov, x).T)
    if verbose: print("\t%16s: %3.3f" % ("Planck 2018 CMB", _logL))
    if use_planck: logL += _logL
    
    return logL


def run_mcmc(pnames, params0, priors):
    """
    Run MCMC sampler for a given model.
    
    Parameters
    ----------
    pnames : list of str
        Names of parameters to sample. These must exist as keys in the params0 
        dictionary.
    
    params0 : dict
        Dictionary containing full set of parameters, including default values. 
        If a parameter is not being sampled, it will be fixed to its value in 
        this dict.
    
    priors : dict
        Dictionary containing prior ranges for a sub-set of parameters.
    """
    # Get initial parameter values and starting value of log-likelihood
    p0 = np.array([params0[pp] for pp in pnames])
    ndim = p0.size
    logl0 = loglike(p0, pnames, params0, priors)
    
    # Get random initial positions for walkers (best-fit values x some O(1) factor)
    p0 = np.outer(np.ones(NWALKERS), p0)
    p0 *= np.random.normal(loc=1., scale=0.02, size=p0.shape)
    p0 += np.random.normal(loc=0., scale=1e-5, size=p0.shape)
    # Don't allow initial parameter values to get stuck at zero
    
    # Initialise emcee sampler and write header of chain file
    sampler = emcee.EnsembleSampler(NWALKERS, ndim, loglike, 
                                    args=(pnames, params0, priors), 
                                    threads=NTHREADS)
    f = open(CHAIN_FILE, "w")
    f.write("# %s %s %s\n" % ("walker", "logl", " ".join(pnames)))
    f.close()

    # Iterate over samples
    nsteps = NSAMPLES
    tstart = time.time()
    print("Starting %d samples with %d walkers and %d threads." \
           % (nsteps, NWALKERS, NTHREADS))
    
    for i, state in enumerate(sampler.sample(p0, iterations=nsteps)):
        # Save current sample to disk
        position = state.coords
        prob = state.log_prob #result[1]
        f = open(CHAIN_FILE, "a")
        for k in range(NWALKERS):
            pvals = " ".join(["%s" % x for x in position[k]])
            f.write("%d %f %s\n" % (k, prob[k], pvals))
        f.close()
        
        # Print status
        if (i+1) % 50 == 0:
            print("Step %d / %d done in %3.1f sec" \
                    % (i+1, nsteps, time.time() - tstart))
            print("  ", ", ".join([pn for pn in pnames]))
            print("  ", ", ".join(["%3.3f" % pv for pv in position[0]]))
            print("  ", "%3.3e" % prob[k])
        tstart = time.time()
    
    print("Done.")


if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        use_expts = expt_sets[int(sys.argv[1])]
        print("Adding experiments:", use_expts)
    else:
        print("MUST PASS AN ARGUMENT!")
        sys.exit(1)

    # Construct chain filename
    CHAIN_FILE = "%s" % CHAIN_FILE_ROOT
    if use_planck: CHAIN_FILE += "_cmb"
    if use_lss: CHAIN_FILE += "_lss"
    for expt in use_expts: CHAIN_FILE += "_%s" % expt.lower()
    CHAIN_FILE += ".dat"
    
    # CMB: Planck 2015, Gaussianised
    pl15_mean, pl15_icov = load_planck_data("planck_derived_fisher_distances", 
                                      params=['omegabh2', 'omegamh2', 'DAstar'])
    pl18_mean, pl18_icov = load_planck_data(
                                 "planck_derived_fisher_distances_2018_omegak", 
                                      params=['omegabh2', 'omegamh2', 'DAstar'] )

    # Get data for all specified Fisher-based experiments
    expt_data = {e: load_fisher_data(fisher_mats[e]) for e in use_expts}
    
    # Set (uniform) prior ranges
    priors = {
        #'w0':       (-1., -0.1), # FIXME: Enable for Mocker model
        'w0':       (-2., -0.1), # FIXME: Enable for tanh model
        #'winf':     (-2., -0.1),
        'deltaw':   (-2., 2.),
        'zc':       (-0.2, 10.), #(-0.2, 1000.),
        'deltaz':   (0.01, 10.), #(0.01, 3.), #(0.01, 5.),
        'omegaB':   (0.01, 0.1),
        'omegaM':   (0.2, 0.4),
        'omegaK':   (-0.2, 0.2),
        'h':        (0.5, 0.8),
        'Cpow':     (1.3, 1.7),
    }
    
    # Set initial/default parameter values
    params0 = {
        'w0':       -0.99,
        #'winf':     -0.8,
        'deltaw':   0.1,
        'zc':       2.0, #1e5
        'deltaz':   0.5,
        'omegaB':   0.0493,
        'omegaM':   0.3153,
        'omegaK':   0.0,
        'h':        0.6736,
        'Cpow':     1.5,
    }
    
    # Set which parameters should be sampled
    #pnames = ['h', 'omegaB', 'omegaM', 'w0',]
    ##pnames = ['omegaM', 'w0', 'h', 'omegaB', 'winf', 'zc', 'deltaz']
    pnames = ['omegaM', 'w0', 'h', 'omegaB', 'deltaw', 'zc', 'deltaz']
    #pnames = ['omegaM', 'w0', 'h', 'omegaB', 'omegaK', 'winf', 'zc', 'deltaz']
    #pnames = ['omegaM', 'w0', 'omegaB', 'winf', 'zc', 'deltaz']
    if MOCKER: pnames = ['h', 'omegaB', 'omegaM', 'w0', 'Cpow']
    
    pvals = [params0[pn] for pn in pnames]
    if MOCKER: params0['mocker'] = True
    
    #omegaM, omegaK, w0, h, deltaz, omegaB, winf, z_eq, zc
    logL = loglike(pvals, pnames, params0, priors, verbose=True)
    print("\t%16s: %3.3f" % ("TOTAL", logL))
    print("\tOutput file:", CHAIN_FILE)
    
    # Run the MCMC
    run_mcmc(pnames, params0, priors)
    
    
