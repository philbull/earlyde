#!/usr/bin/python
"""
Run MCMC on tanh model.
"""
import numpy as np
import pylab as P
import copy, time
import emcee
import wztanh as model
import classy

# MCMC sampler settings
NTHREADS = 4
NSAMPLES = 1500
NWALKERS = 40
#CHAIN_FILE = "wcdm_cmbonly.dat"
#CHAIN_FILE = "wcdm_cmb_lss.dat"
#CHAIN_FILE = "wcdm_cmb_lss_lya.dat"
#CHAIN_FILE = "wztanh_zc100_cmb_lss.dat"
CHAIN_FILE = "wztanh_cmb_lss_MADEUP.dat"
CHAIN_FILE = "wzmocker_cmb_lss.dat"

np.random.seed(10)


def load_planck_data(froot, params=['omegabh2', 'omegamh2', 'DAstar']):
    """
    Load Gaussianised Planck data from file.
    The parameters will be ordered in the output mean and C^-1 arrays according 
    to their order in the 'params' list.
    """
    f_mean = "%s.mean.dat" % froot
    f_cov = "%s.cov.dat" % froot
    
    # Load header
    f = open(f_mean, 'r')
    hdr = f.readline()
    f.close()
    names = hdr[2:-1].split(' ')
    
    # Load mean, cov, inverse cov
    mean = np.genfromtxt(f_mean)
    cov = np.genfromtxt(f_cov)
    
    # Select only certain parameters
    # We want to *marginalise* over the non-included parameters, which means 
    # we only keep the elements of the *covariance* (not C^-1) for the 
    # parameters that we do want to include.
    new_mean = np.array([mean[names.index(pn)] for pn in params])
    new_cov = np.zeros((len(params), len(params)))
    for i, pni in enumerate(params):
        for j, pnj in enumerate(params):
            new_cov[i,j] = cov[names.index(pni), names.index(pnj)]
    
    # The correct covariance
    icov = np.linalg.inv(new_cov)
    return new_mean, icov


# Data
r_d = 147.33 #147.49 # \pm 0.59 Mpc [Planck LCDM Mnu and N_eff, see p6 of 1411.1074]
lss_data = [
    ('6dFGS', 'DV', 0.106, 3.047, 0.137),
    ('MGS', 'DV', 0.15, 4.480, 0.168),
    ('BOSS LOWZ', 'DV', 0.32, 8.467, 0.167),
    ('BOSS CMASS', 'DM', 0.57, 14.945, 0.210),
    ('BOSS CMASS', 'DH', 0.57, 20.75, 0.73),
    #('LyaF auto', 'DM', 2.34, 37.675, 2.171),
    #('LyaF auto', 'DH', 2.34, 9.18, 0.28),
    #('LyaF-QSO', 'DM', 2.36, 36.288, 1.344),
    #('LyaF-QSO', 'DH', 2.36, 9.00, 0.30),
    
    #('MADEUP', 'DH', 4.0, 4.789, 0.004789), # FIXME: Just took LCDM at z=4
    #('MADEUP', 'DM', 2.36, 36.288, 0.0344), # FIXME
    
    ('CMB approx', 'DM', 1090., 94.51, np.sqrt(0.004264)) # FIXME
]

"""
# Audberg et al.
cmb_data = np.array([0.02245, 0.1386, 94.33]) # omega_b, omega_cb, D_M(1090)/r_d
cmb_cov = np.array([ [ 1.286e-7, -6.033e-7,  1.443e-5],
                     [-6.033e-7,  7.542e-6, -3.605e-5],
                     [ 1.443e-5, -3.605e-5,  0.004264] ])
cmb_icov = np.linalg.inv(cmb_cov)

# Planck 2015 (1D marginals)
plnck_rstar = (144.61, 0.49)
plnck_theta_s = (1.04105, 0.00046) # theta_* = r_s(z*) / D_A(z*)
plnck_omegab = (0.02222, 0.00023)
plnck_omegam = (0.1426, 0.0020)
"""

# Planck 2015 Gaussianised
pl15_mean, pl15_icov = load_planck_data("planck_derived_fisher_distances", 
                                      params=['omegabh2', 'omegamh2', 'DAstar'])


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
    
    # Calculate model values for these datapoints
    # FIXME: Ignoring the covariance between H ad D_A for now
    dv, dm, dh = model.lss_distances(np.array(zc), p)
    model_calc = {'DV': dv, 'DM': dm, 'DH': dh}
    
    # Calculate simple independent Gaussian log-likelihoods for LSS
    logL = 0
    for i in range(len(dname)):
        if 'CMB' in dname[i]: continue # Skip CMB for now
        dist = model_calc[dtype[i]] # Get distance measure for this data point
        _logL = -0.5*(dval[i] - dist[i]/r_d)**2. / derr[i]**2.
        if verbose: print "\t%10s: %3.3f" % (dname[i], _logL)
        logL += _logL
    
    # Find where CMB item is in the list of calculated distances
    cmb_idxs = [i for i in range(len(dname)) if 'CMB' in dname[i]]
    if len(cmb_idxs) != 1: raise KeyError("Should be only 1 CMB item in lss_data.")
    idx = cmb_idxs[0]
    assert np.abs(zc[idx] - 1090.) < 10., "CMB datapoint found at wrong redshift!"
    
    """
    # Audberg Planck values
    # Construct CMB parameter model vector: omega_b, omega_cb, DM(1090)/r_d
    h = p['h']
    cmb_model = np.array([p['omegaB']*h**2., p['omegaM']*h**2., dm[idx]/r_d])
    _logL = -0.5 * np.dot(cmb_model, np.dot(cmb_icov, cmb_model))
    if verbose: print "\t%10s: %3.3f" % ("CMB", _logL)
    logL += _logL / 1e4 # FIXME: Seriously down-weighting CMB
    # FIXME: Cut out the CMB!
    """
    """
    # Planck 2015 1D marginal CMB
    h = p['h']
    theta_s_model = plnck_rstar[0] / dm[idx]
    _logL = -0.5 * ((theta_s_model*100. - plnck_theta_s[0]) / plnck_theta_s[1])**2.
    _logL += -0.5* ((p['omegaM']*h**2. - plnck_omegam[0]) / plnck_omegam[1])**2.
    _logL += -0.5* ((p['omegaB']*h**2. - plnck_omegab[0]) / plnck_omegab[1])**2.
    """
    #print (theta_s_model*100. - plnck_theta_s[0]) / plnck_theta_s[1]
    
    # Planck 2015 Gaussianised likelihood
    # Assumed order is 'omegabh2', 'omegamh2', 'DAstar'
    h = p['h']
    model_vec = np.array([p['omegaB']*h**2., p['omegaM']*h**2., dm[idx]])
    x = model_vec - pl15_mean
    _logL = -0.5 * np.dot(x, np.dot(pl15_icov, x).T)
    if verbose: print "\t%10s: %3.3f" % ("CMB", _logL)
    logL += _logL
    
    """
    # FIXME: Test CLASS
    class_params = {
        "output"    : "",
        "T_cmb"     : 2.725,
        "h"         : p['h'],
        "Omega_cdm" : p['omegaM']-p['omegaB'],
        "Omega_b"   : p['omegaB'],
        "A_s"       : 2e-9,
        "n_s"       : 1.0,
        #"w0_fld"    : -1., #p['w0'],
        #"wa_fld"    : 0.0,
        "Omega_k"  : 0.0,
        "N_ur"     : 3.0, # normally 3
        "N_ncdm"   : 1,   # 1 species
    }
    cosm = classy.Class()
    cosm.set(class_params)
    cosm.compute()
    model_vec = np.array([p['omegaB']*h**2., p['omegaM']*h**2., 
                          (1.+1090.)*cosm.angular_distance(1090.)])
    x = model_vec - pl15_mean
    print cosm.angular_distance(1090.)*(1.+1090.), dm[idx]
    cosm.struct_cleanup()
    
    _logL = -0.5 * np.dot(x, np.dot(pl15_icov, x).T)
    if verbose: print "\t%10s: %3.3f" % ("CMB CLASS", _logL)
    logL += _logL
    """
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
    p0 *= np.random.normal(loc=1., scale=0.00005, size=p0.shape)
    
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
    print "Starting %d samples with %d walkers and %d threads." \
           % (nsteps, NWALKERS, NTHREADS)
    
    for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
        # Save current sample to disk
        position = result[0]
        prob = result[1]
        f = open(CHAIN_FILE, "a")
        for k in range(NWALKERS):
            pvals = " ".join(["%s" % x for x in position[k]])
            f.write("%d %f %s\n" % (k, prob[k], pvals))
        f.close()
        
        # Print status
        if (i+1) % 50 == 0:
            print "Step %d / %d done in %3.1f sec" \
                    % (i+1, nsteps, time.time() - tstart)
            print "  ", ", ".join([pn for pn in pnames])
            print "  ", ", ".join(["%3.3f" % pv for pv in position[0]])
            print "  ", "%3.3e" % prob[k]
        tstart = time.time()
    
    print "Done."


if __name__ == '__main__':
    
    # Set (uniform) prior ranges
    priors = {
        'w0':       (-2., -0.1),
        'winf':     (-2., -0.1),
        'zc':       (-0.2, 10.), #(-2., 1000.),
        'deltaz':   (0.01, 3.), #(0.01, 5.),
        'omegaB':   (0.01, 0.1),
        'omegaM':   (0.25, 0.36),
        'omegaK':   (-0.2, 0.2),
        'h':        (0.5, 0.8),
        'Cpow':     (1.3, 1.7),
    }
    
    # Set initial/default parameter values
    params0 = {
        'w0':       -1.,
        'winf':     -0.8,
        'zc':       2.0, #1e5
        'deltaz':   0.5,
        'omegaB':   0.045,
        'omegaM':   0.3183,
        'omegaK':   0.0,
        'h':        0.6704,
        'Cpow':     1.5,
    }
    
    # Set which parameters should be sampled
    #pnames = ['h', 'omegaB', 'omegaM', 'w0',]
    #pnames = ['omegaM', 'w0', 'h', 'omegaB', 'winf', 'zc', 'deltaz']
    #pnames = ['omegaM', 'w0', 'h', 'omegaB', 'omegaK', 'winf', 'zc', 'deltaz']
    #pnames = ['omegaM', 'w0', 'omegaB', 'winf', 'zc', 'deltaz']
    pnames = ['h', 'omegaB', 'omegaM', 'w0', 'Cpow']
    
    pvals = [params0[pn] for pn in pnames]
    
    params0['mocker'] = True
    
    #omegaM, omegaK, w0, h, deltaz, omegaB, winf, z_eq, zc
    #logL = loglike(pvals, pnames, params0, priors, verbose=True)
    #print logL
    
    # Run the MCMC
    run_mcmc(pnames, params0, priors)
    
    
