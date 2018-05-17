#!/usr/bin/python
"""
Test that the angular diameter distance to the CMB is calculated accurately.
"""
import numpy as np
import pylab as P
import pyccl as ccl
import wztanh as model

# Default parameters
params = {
    'w0':       -1.,
    'winf':     -1.,
    'zc':       2.0,
    'deltaz':   0.5,
    'omegaM':   0.315,
    'omegaK':   0.0,
    'omegaB':   0.045,
    'h':        0.6731,
    'z_eq':     3393.,
}

# Initialise matching CCL cosmology
p = params
cosmo = ccl.Cosmology(Omega_c=p['omegaM']-p['omegaB'], Omega_b=p['omegaB'], 
                      h=p['h'], Omega_k=0., w0=-1., A_s=2.1e-9, n_s=0.96,
                      m_nu=0.06)

# Plot low-redshift D_A(z) and H(z)
z = np.logspace(-3., 3.1, 1000)
a = 1. / (1. + z)

# CCL quantities
ccl_DA = a * ccl.comoving_angular_distance(cosmo, a)
ccl_Hz = (100. * p['h']) * ccl.h_over_h0(cosmo, a)

# My quantities
this_Hz = model.Hz(a, params)
this_DA = model.DAz(a, params)

#-------------------------------------------------------------------------------

# Load CCL benchmark file
dat = np.genfromtxt("CCL_TEST_CHI.dat").T
bench_z = dat[0]
#bench_chi = dat[1] # flat, no massive neutrinos, 3 massless
bench_chi = dat[3] # flat, 1 massive neutrino, 0.1 eV

#P.title("massless")
#P.title("0.1 eV")

params_bench = {
    'w0':       -1.,
    'winf':     -1.,
    'zc':       2.0,
    'deltaz':   0.5,
    'omegaM':   0.30,
    'omegaK':   0.0,
    'omegaB':   0.05,
    'h':        0.7,
    'z_eq':     3393.,
}
bench_a = 1./(1. + bench_z)
this_DA_bench = model.DAz(bench_a, params_bench)

p = params_bench
cosmo_bench = ccl.Cosmology(Omega_c=p['omegaM']-p['omegaB'], Omega_b=p['omegaB'], 
                            h=p['h'], Omega_k=0., w0=-1., A_s=2.1e-9, n_s=0.96,
                            m_nu=0.1)
ccl_DA_bench = bench_a * ccl.comoving_angular_distance(cosmo_bench, bench_a)


#-------------------------------------------------------------------------------
# CLASS calculations
import classy

p = params
class_params = {
    "output"    : "",
    "T_cmb"     : 2.725,
    "h"         : p['h'],
    "Omega_cdm" : p['omegaM']-p['omegaB'],
    "Omega_b"   : p['omegaB'],
    "A_s"       : 2e-9,
    "n_s"       : 1.0,
    # "w0_fld"    : -1.0,
    # "wa_fld"    : 0.0,
    "back_integration_stepsize" : 1.0e-3, # Req. for agreement at <1e-5 with CCL.
    "Omega_k"  : 0.0,
    "N_ur"     : 2.0, # normally 3
    "N_ncdm"   : 1,   # 1 species
    "m_ncdm"   : 0.06
}

cosm = classy.Class()
cosm.set(class_params)
cosm.compute()
print cosm.angular_distance(1090.)
print cosm.Omega_g()*cosm.h()**2. #, cosmo.Omega0_m()

#help(cosm)
exit()

z_class = cosm.get_background()["z"][:]
a_class = 1./(1. + z_class)
DA_class = cosm.get_background()["ang.diam.dist."][:]
#d_co = cosm.get_background()["comov. dist."][:]
cosm.struct_cleanup()

# Calculate using our code
this_DA_class = model.DAz(a_class, params)


#-------------------------------------------------------------------------------

# Plot fractional difference
P.subplot(111)

P.axhline(0., color='k', alpha=0.2, lw=1.8)

P.plot(z, this_Hz/ccl_Hz - 1., 'b-', lw=1.8, label="H(z)")
P.plot(z, this_DA/ccl_DA - 1., 'r-', lw=1.8, label="D_A(z)")

P.plot(bench_z, this_DA_bench/(bench_a*bench_chi) - 1., 
       'c--', lw=1.8, label="D_A This vs CLASS")

P.plot(bench_z, ccl_DA_bench/(bench_a*bench_chi) - 1., 
       'm--', lw=1.8, label="D_A CCL vs CLASS")

P.plot(z_class, this_DA_class/DA_class - 1., 'y-', lw=1.8, label="This vs. CLASS")

#P.ylim((-0.02, 0.02))
P.ylim((-0.002, 0.002))

P.legend(loc='upper right')
P.xscale('log')
P.tight_layout()
P.show()
