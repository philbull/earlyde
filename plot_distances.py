#!/usr/bin/env python
"""
Plot distance measure samples and data points
"""
import numpy as np
import pylab as P
import wztanh as model
import load_data_files
import copy

EXPT = 4

fnames = [ "final2_wztanh_cmb_lss.dat", 
           "final_wztanh_cmb_lss_desi.dat", 
           "final2_wztanh_cmb_lss_hirax_hw.dat",
           "final2_wztanh_cmb_lss_hirax_pbw.dat",
           "final2_wztanh_cmb_lss_hirax_nw.dat" ]
labels = [ 'tanh CMB + LSS', 
           'tanh CMB + LSS + DESI', 
           'tanh CMB + LSS + HIRAX (horiz. wedge)',
           'tanh1 CMB + LSS + HIRAX (3PB wedge)',
           'tanh CMB + LSS + HIRAX (nPB wedge)' ]

def load_chain(fname):
    """
    Load MCMC chain output from file
    """
    print("Loading", fname)
    # Load header
    f = open(fname, 'r')
    hdr = f.readline()
    f.close()
    names = hdr[2:-1].split(' ')
    
    # Load data
    dat = np.genfromtxt(fname).T
    data = {}
    for i, n in enumerate(names):
        data[n] = dat[i]
    return data


# Load data and select random samples
dat = load_chain(fnames[EXPT])
np.random.seed(10)
print(dat.keys())

# Get random set of samples
idxs = np.random.randint(low=100000, high=dat['h'].size, size=500)

# Scale factor array
#z = np.logspace(-4., np.log10(1090.), 500)
z = np.linspace(0., 10., 500)
a = 1./(1.+z)

# Calculate distance measures for a random subset of parameters
Hz = []; DAz = []; wz = []
for i in idxs:
    p = {key: dat[key][i] for key in dat.keys()}
    p['omegaK'] = 0.
    #p['mocker'] = True
    #P.plot(z, model.Hz(a, p) / (1.+z), 'b-', lw=1.8, alpha=0.1)
    Hz.append( model.Hz(a, p) )
    DAz.append( model.DAz(a, p) )
    wz.append( model.wz(a, p) )
    omegaK = 0.
    #OmegaDE0 = 1. - p['omegaM'] - model.omegaR(p) - omegaK
    #P.plot(z, OmegaDE0*model.omegaDE(a, p), 'b-', lw=1.8, alpha=0.1)
Hz = np.array(Hz)
DAz = np.array(DAz)
wz = np.array(wz)

pcts = [2.5, 16., 50., 84., 97.5]
print(dat['h'].size)

pct = np.percentile(Hz/(1.+z), pcts, axis=0)
pct_da = np.percentile(DAz, pcts, axis=0)
pct_wz = np.percentile(wz, pcts, axis=0)
#for i in range(pct.shape[0]):
#    P.plot(z, pct[i], 'b-', lw=1.8, alpha=1. - np.abs(pcts[i] - 50.)/70.)

plcdm = model.pdict(h=0.6727, omegaM=0.3166, omegaK=0.0, omegaB=0.04941, 
                    w0=-1., winf=-1., zc=1e5, deltaz=0.5)
#OmegaDE0 = 1. - plcdm['omegaM'] - model.omegaR(plcdm) - omegaK
#P.plot(z, OmegaDE0*model.omegaDE(a, plcdm), 'k-', lw=1.8, alpha=1.)

# Some other model (within bounds)
ptanh = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
                    w0=-1., winf=-0.9, zc=3, deltaz=2.5)
ptanh2 = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
                     w0=-1., winf=-0.9, zc=1, deltaz=2.5)


# Load Fisher matrix
print("Loading Fisher data")
if EXPT == 2: 
    expt_zc, expt_mean, expt_icov \
        = load_data_files.load_fisher_data("Fisher-full-iHIRAX_2yr_horizwedge")

if EXPT == 3: 
    expt_zc, expt_mean, expt_icov \
        = load_data_files.load_fisher_data("Fisher-full-iHIRAX_2yr_3pbwedge")

if EXPT == 4:
    expt_zc, expt_mean, expt_icov \
        = load_data_files.load_fisher_data("Fisher-full-iHIRAX_2yr")

expt_cov = np.linalg.inv(expt_icov)

# Load Planck covariance
pl18_mean, pl18_icov = load_data_files.load_planck_data(
                                 "planck_derived_fisher_distances_2018_omegak", 
                                      params=['omegabh2', 'omegamh2', 'DAstar'] )
pl18_cov = np.linalg.inv(pl18_icov)

print("Planck std:", np.sqrt(np.diag(pl18_cov)))


print(expt_zc)
print(expt_mean)
print(expt_icov.shape)

# Plot H(z)
P.subplot(121)
P.title(labels[EXPT])

# MCMC constraints
yy = model.Hz(a, plcdm) / (1.+z)
P.fill_between(z, pct[0] - yy, pct[-1] - yy, color='r', alpha=0.2, lw=1.)
P.fill_between(z, pct[1] - yy, pct[-2] - yy, color='r', alpha=0.4, lw=1.)
P.plot(z, pct[2] - yy, color='r', lw=1.5)

# LCDM theory model
#P.plot(z, yy, 'k-', lw=1.8, alpha=1.)

P.plot(z, model.Hz(a, ptanh) / (1.+z) - yy, 'c-', lw=1.8, alpha=1.)
P.plot(z, model.Hz(a, ptanh2) / (1.+z) - yy, 'm-', lw=1.8, alpha=1.)

# Input data from Fisher matrix
for i in range(expt_zc.size):
    yyc = model.Hz(1./(1.+expt_zc[i]), plcdm) / (1.+expt_zc[i])
    Hz_fisher = expt_mean[i] * 1e2 / (1. + expt_zc[i])
    Hz_err = np.sqrt(np.diag(expt_cov))[i] * 1e2 / (1. + expt_zc[i])
    P.errorbar(expt_zc[i], Hz_fisher - yyc, yerr=Hz_err, color='k', marker='o')

P.xlabel("z", fontsize=15)
P.ylabel("H(z)", fontsize=15)
#P.ylim((58., 77.))
P.ylim((-0.4, 0.6))
#P.xscale('log')


# Plot D_A(z)
P.subplot(122)

## MCMC constraints
#P.fill_between(z, pct_da[0], pct_da[-1], color='r', alpha=0.2, lw=1.)
#P.fill_between(z, pct_da[1], pct_da[-2], color='r', alpha=0.4, lw=1.)
#P.plot(z, pct_da[2], color='r', lw=1.5)

## LCDM theory model
#P.plot(z, model.DAz(a, plcdm), 'k-', lw=1.8, alpha=1.)

# MCMC constraints
P.fill_between(z, pct_wz[0], pct_wz[-1], color='r', alpha=0.2, lw=1.)
P.fill_between(z, pct_wz[1], pct_wz[-2], color='r', alpha=0.4, lw=1.)
P.plot(z, pct_wz[2], color='r', lw=1.5)

# LCDM theory model
P.plot(z, model.wz(a, plcdm), 'k-', lw=1.8, alpha=1.)

# Some other model
P.plot(z, model.wz(a, ptanh), 'b-', lw=1.8, alpha=1.)
P.plot(z, model.wz(a, ptanh2), 'm-', lw=1.8, alpha=1.)

P.ylim((-1.2, -0.7))

P.tight_layout()
P.show()

