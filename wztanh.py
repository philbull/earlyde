#!/usr/bin/python
"""
Integrate a simple tanh model for the equation of state of DE.
"""
import numpy as np
import pylab as P
import copy
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

C = 2.99792458e5 # km/s
ZMAX = 1600. # Max. redshift to integrate distances to
NSAMP_OMEGADE = 500
NSAMP_CHI = 1000

default_params = {
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
params = default_params

# Data
r_d = 147.33 #147.49 # \pm 0.59 Mpc [Planck LCDM Mnu and N_eff, see p6 of 1411.1074]
data = [
    ('6dFGS', 'DV', 0.106, 3.047, 0.137),
    ('MGS', 'DV', 0.15, 4.480, 0.168),
    ('BOSS LOWZ', 'DV', 0.32, 8.467, 0.167),
    ('BOSS CMASS', 'DM', 0.57, 14.945, 0.210),
    ('BOSS CMASS', 'DH', 0.57, 20.75, 0.73),
    ('LyaF auto', 'DM', 2.34, 37.675, 2.171),
    ('LyaF auto', 'DH', 2.34, 9.18, 0.28),
    ('LyaF-QSO', 'DM', 2.36, 36.288, 1.344),
    ('LyaF-QSO', 'DH', 2.36, 9.00, 0.30),
    ('CMB approx', 'DM', 1090., 94.51, np.sqrt(0.004264)) # FIXME
]

def pdict(**kwargs):
    """
    Build a parameter dictionary from named kwargs.
    """
    p = copy.deepcopy(default_params)
    for key in kwargs.keys(): p[key] = kwargs[key]
    return p

def wz(a, params):
    """
    Equation of state, parametrised as a tanh model that interpolates between 
    w0 at z~0 and winf as z -> infty.
    """
    z = 1./a - 1.
    
    # Mocker model
    if 'mocker' in params.keys() and params['mocker']:
        w0 = params['w0']; Cpow = params['Cpow']
        w = w0 * a**Cpow / (w0*(1. - a**Cpow) + 1.)
        return w
    
    # tanh model
    w0 = params['w0']
    zc = params['zc']
    deltaz = params['deltaz']
    if 'winf' in params.keys():
        winf = params['winf']
    elif 'deltaw' in params.keys():
        winf = w0 - params['deltaw']
    else:
        raise ValueError("Parameters 'winf' or 'deltaw' must be specified.")
    return w0 + 0.5 * (winf - w0) * (np.tanh((z - zc) / deltaz) + 1.)
    
def omegaR(params):
    """
    Fractional density of radiation (incl. neutrinos).
    """
    Neff = 3.046
    #OmegaR = OmegaM / (1. + z_eq) #* (1. + Neff*(7./8.)*(4./11.)**(4./3.))
    Omega_g = 2.471165e-05 / params['h']**2. # FIXME: Copied this number from CCL
    OmegaR = Omega_g * (1. + Neff*(7./8.)*(4./11.)**(4./3.))
    return OmegaR

def omegaDE(a=None, params=None):
    """
    Fractional density of dark energy, Omega_DE(a) = rho_DE(a) / rho_DE(a=1).
    """
    # Scale factor sample points, in reverse order
    aa = np.logspace(np.log10(1./(1.+ZMAX)), 0., NSAMP_OMEGADE)
    
    # Integrate to get argument of exponential
    w = wz(aa, params)
    integ = -3.*(1. + w)/aa
    log_rho = cumtrapz(integ[::-1], aa[::-1], initial=0.)[::-1]
    
    # Interpolate and return Omega_DE(a) at requested scale factor values
    interp_logrho = interp1d(aa, log_rho, kind='linear', bounds_error=False)
    if a is not None:
        # Interpolate directly
        return np.exp(interp_logrho(a))
    else:
        fn = lambda a: np.exp(interp_logrho(a))
        return fn

def Hz(a, params):
    """
    Hubble rate, H(z), in km/s/Mpc units.
    """
    h = params['h']
    OmegaM = params['omegaM']
    OmegaK = params['omegaK']
    #z_eq = params['z_eq']
    
    # Massive neutrinos
    Omega_numass = (0.06/93.) / h**2.
    #OmegaM += Omega_numass
    
    # Radiation fractional density at z=0
    OmegaR = omegaR(params)
    
    # Set OmegaDE,0 using a constraint
    OmegaDE0 = 1. - OmegaM - OmegaR - OmegaK
    
    # Dimensionless Hubble rate squared
    E2 = OmegaR * a**-4. + OmegaM * a**-3. + OmegaK * a**-2. \
       + OmegaDE0 * omegaDE(a, params)
    
    # Return H(z) in km/s/Mpc units
    return (h * 100.) * np.sqrt(E2)

def chi(a, params):
    """
    Comoving radial distance, in Mpc.
    """
    # Ensure denser sampling near a=1, where errors can get large
    aa = np.concatenate(
            (np.logspace(np.log10(1./(1.+ZMAX)), np.log10(0.7), NSAMP_CHI/2),
             np.logspace(np.log10(0.71), 0., NSAMP_CHI/2)) )
    
    integ = 1./(aa**2.) / (Hz(aa, params)/100.) # Convert H -> H/(100 km/s)
    intchi = -cumtrapz(integ[::-1], aa[::-1], initial=0.)[::-1]
    
    # Interpolate unscaled integral for chi
    interp_intchi = interp1d(aa, intchi, kind='linear', bounds_error=False)
    
    return (C / 100.) * interp_intchi(a)
    
def DAz(a, params):
    """
    Angular diameter distance, in Mpc.
    """
    ok = params['omegaK']
    d_hub = C / (100. * params['h']) # C / H_0 in Mpc
    
    # Get comoving radial dist. and transform according to sign of curvature
    _chi = chi(a, params)
    if np.abs(ok) < 1e-5: 
        return a * _chi # Flat
    if ok > 0.:
        return d_hub * a * np.sinh(np.sqrt(ok)*_chi/d_hub) / np.sqrt(ok) # Open
    if ok < 0.:
        return d_hub * a * np.sin(np.sqrt(-ok)*_chi/d_hub) / np.sqrt(-ok) # Closed
    raise ValueError("OmegaK calculation invalid.")


def lss_distances(z, params):
    """
    Calculates large-scale structure distances D_V(z), D_M(z), and D_H(z).
    """
    a = 1. / (1. + z)
    hz = Hz(a, params)
    da = DAz(a, params)
    
    # Calculate distance combinations
    DM = da / a # Comoving angular diameter distance
    DH = C / hz # Hubble distance
    DV = (z * DH * DM**2.)**(1./3.)
    return DV, DM, DH
    

def compare_with_ccl():
    """
    Compare some results with CCL to check for accuracy.
    """
    # Set-up CCL
    import pyccl as ccl
    p = default_params
    OmegaB = 0.045
    cosmo1 = ccl.Cosmology(Omega_c=p['omegaM']-OmegaB, Omega_b=OmegaB, h=p['h'], 
                           Omega_k=0., w0=-1., wa=0., A_s=2.1e-9, n_s=0.96)
    cosmo2 = ccl.Cosmology(Omega_c=p['omegaM']-OmegaB, Omega_b=OmegaB, h=p['h'], 
                           Omega_k=0.02, w0=-1., wa=0., A_s=2.1e-9, n_s=0.96)
    cosmo3 = ccl.Cosmology(Omega_c=p['omegaM']-OmegaB, Omega_b=OmegaB, h=p['h'], 
                           Omega_k=0.02, w0=-0.8, wa=0., A_s=2.1e-9, n_s=0.96,
                           m_nu=0.)
    
    # Set z_eq to include massless neutrinos
    p['z_eq'] = p['omegaM'] \
              / (cosmo1.params['Omega_g'] + cosmo1.params['Omega_n_rel']) - 1.
    print(p['z_eq'])
    
    # Expansion rate
    a = np.logspace(-3., 0., 400)
    ccl_H1 = (100. * p['h']) * ccl.h_over_h0(cosmo1, a)
    ccl_H2 = (100. * p['h']) * ccl.h_over_h0(cosmo2, a)
    ccl_H3 = (100. * p['h']) * ccl.h_over_h0(cosmo3, a)
    
    this_H1 = Hz(a, pdict(w0=-1., winf=-1., omegaK=0., zc=2., deltaz=0.5))
    this_H2 = Hz(a, pdict(w0=-1., winf=-1., omegaK=0.02, zc=2., deltaz=0.5))
    this_H3 = Hz(a, pdict(w0=-0.8, winf=-0.8, omegaK=0.02, zc=2., deltaz=0.5))
    
    # Angular diameter distance
    ccl_DA1 = ccl.luminosity_distance(cosmo1, a) * a**2.
    ccl_DA2 = ccl.luminosity_distance(cosmo2, a) * a**2.
    ccl_DA3 = ccl.luminosity_distance(cosmo3, a) * a**2.
    
    this_DA1 = DAz(a, pdict(w0=-1., winf=-1., omegaK=0., zc=2., deltaz=0.5))
    this_DA2 = DAz(a, pdict(w0=-1., winf=-1., omegaK=0.02, zc=2., deltaz=0.5))
    this_DA3 = DAz(a, pdict(w0=-0.8, winf=-0.8, omegaK=0.02, zc=2., deltaz=0.5))
    
    # Plot comparisons
    P.subplot(121)
    P.plot(a, this_H1 / ccl_H1 - 1., 'k-', lw=1.8)
    P.plot(a, this_H2 / ccl_H2 - 1., 'r-', lw=1.8)
    P.plot(a, this_H3 / ccl_H3 - 1., 'b-', lw=1.8)
    P.xscale('log')
    
    P.subplot(122)
    P.plot(a, this_DA1 / ccl_DA1 - 1., 'k-', lw=1.8)
    P.plot(a, this_DA2 / ccl_DA2 - 1., 'r-', lw=1.8)
    P.plot(a, this_DA3 / ccl_DA3 - 1., 'b-', lw=1.8)
    P.xscale('log')
    
    P.tight_layout()
    P.show()


def plot_lss_distances():
    """
    Plot D_V, D_M, and D_H, as in Fig. 1 of 1411.1074.
    """
    P.subplot(111)
    colours = {'DV': 'b', 'DM': 'r', 'DH': 'g'}

    for dp in data:
        dname, dtype, zc, dval, derr = dp # Data have already been divided by r_d
        fac = np.sqrt(zc) # Scaling factor from Fig. 1 of 1411.1074.
        if dtype == 'DH': fac /= zc
        
        P.errorbar(zc, dval/fac, yerr=derr/fac, color=colours[dtype], label=dname, 
                   marker='.', mew=1.8, capsize=5.)

    # Plot theory curves
    #z = np.logspace(-3., np.log10(3.), 200)
    z = np.logspace(-3., np.log10(1500.), 500)
    fac = r_d * np.sqrt(z)

    # LambdaCDM
    DV, DM, DH = lss_distances(z, pdict(w0=-1., winf=-1., zc=2., deltaz=0.5))
    P.plot(z, DV/fac, lw=1.8, color=colours['DV'])
    P.plot(z, DM/fac, lw=1.8, color=colours['DM'])
    P.plot(z, z*DH/fac, lw=1.8, color=colours['DH'])

    DV, DM, DH = lss_distances(z, pdict(w0=-0.5, winf=-0.5, 
                                        zc=2., deltaz=0.5, omegaK=0.02))
    P.plot(z, DV/fac, lw=1.8, color=colours['DV'], ls='dashed', alpha=0.5)
    P.plot(z, DM/fac, lw=1.8, color=colours['DM'], ls='dashed', alpha=0.5)
    P.plot(z, z*DH/fac, lw=1.8, color=colours['DH'], ls='dashed', alpha=0.5)

    P.xlabel("$z$", fontsize=16.)
    P.ylabel(r"${\rm Distance} / r_d \sqrt{z}$", fontsize=16.)

    #P.xlim((0.08, 3.))
    P.xlim((0.08, 1500.))
    P.ylim((-0.5, 30.))
    P.xscale('log')

    P.legend(loc='upper right', frameon=False, ncol=2)

    P.tight_layout()
    P.show()


def plot_lss_deviations():
    """
    Plot D_V, D_M, and D_H, relative to a LambdaCDM model
    """
    P.subplot(111)
    colours = {'DV': 'b', 'DM': 'r', 'DH': 'g'}

    for dp in data:
        dname, dtype, zc, dval, derr = dp # Data have already been divided by r_d
        
        DV, DM, DH = lss_distances(zc, pdict(w0=-1., winf=-1., zc=2., deltaz=0.5))
        if dtype == 'DV': fac = DV / r_d
        if dtype == 'DM': fac = DM / r_d
        if dtype == 'DH': fac = DH / r_d
        
        P.errorbar(zc, dval/fac - 1., yerr=derr/fac, color=colours[dtype], 
                   label=dname, marker='.', mew=1.8, capsize=5.)

    # Plot theory curves
    #z = np.logspace(-3., np.log10(3.), 200)
    z = np.logspace(-3., np.log10(1500.), 500)
    
    P.axhline(0., color='k', lw=1.8, alpha=0.4)
    
    # LambdaCDM
    DV0, DM0, DH0 = lss_distances(z, pdict(w0=-1., winf=-1., zc=0., deltaz=0.5))
    
    # Alternative model (1)
    pmodel = pdict(w0=-1., winf=-0.5, zc=0.3, deltaz=0.2, omegaK=0.03)
    DV, DM, DH = lss_distances(z, pmodel)
    P.plot(z, DV/DV0 - 1., lw=1.8, color=colours['DV'], ls='dashed', alpha=0.5)
    P.plot(z, DM/DM0 - 1., lw=1.8, color=colours['DM'], ls='dashed', alpha=0.5)
    P.plot(z, DH/DH0 - 1., lw=1.8, color=colours['DH'], ls='dashed', alpha=0.5)
    
    # Alternative model (1)
    #pmodel = pdict(w0=-1., winf=-0.5, zc=0.3, deltaz=0.2, omegaK=0.0)
    
    pmodel = pdict(omegaM = 0.323, omegaK = -0.199, w0 = -1.298, h = 0.743,
                   deltaz = 4.748, omegaB = 0.100, winf = -0.297, 
                   z_eq = 13434.643, zc = 147.897)
    
    DV, DM, DH = lss_distances(z, pmodel)
    P.plot(z, DV/DV0 - 1., lw=1.8, color=colours['DV'], ls='dotted', alpha=0.5)
    P.plot(z, DM/DM0 - 1., lw=1.8, color=colours['DM'], ls='dotted', alpha=0.5)
    P.plot(z, DH/DH0 - 1., lw=1.8, color=colours['DH'], ls='dotted', alpha=0.5)
    
    # FIXME: Planck 2015 CMB
    plcdm = pdict(w0=-1., winf=-1., zc=0., deltaz=0.5, zeq=3394.)
    DVs, DMs, DHs = lss_distances(1090., plcdm)
    ob = (0.02222, 0.00023)
    om = (0.1426, 0.0020)
    rstar = (144.61, 0.49)
    theta_s = (1.04105, 0.00046) # theta_* = r_s(z*) / D_A(z*)
    #theta_s = (1.04085, x) # theta_* = r_s(z*) / D_A(z*)
    
    zq = 3393.
    #DM = da / a
    
    theta_s_model = 144.61 / DMs
    print("Delta(100 theta_s) =", (theta_s_model*100. - theta_s[0]))
    print("Frac(100 theta_s) =", (theta_s_model*100. - theta_s[0])/theta_s[1])
    #print("D_A,model =", DAz(1./(1.+1090.), plcdm))
    
    import pyccl as ccl
    cosmo = ccl.Cosmology(Omega_c=plcdm['omegaM']-plcdm['omegaB'], 
                          Omega_b=plcdm['omegaB'], h=plcdm['h'], 
                          Omega_k=0., w0=-1., wa=0., A_s=2.1e-9, n_s=0.96,
                          m_nu=0.06)
    zstar = 1090.
    da_co = ccl.comoving_angular_distance(cosmo, 1./(1.+zstar))
    theta_s_ccl = 144.61 / da_co
    print("Delta(100 theta_s) =", (theta_s_ccl*100. - theta_s[0]))
    print("Frac(100 theta_s) =", (theta_s_ccl*100. - theta_s[0])/theta_s[1])
    
    exit()
    
    
    """
    # Try CCL plot
    import pyccl as ccl
    p = default_params
    OmegaB = 0.045
    cosmo = ccl.Cosmology(Omega_c=p['omegaM']-OmegaB, Omega_b=OmegaB, h=p['h'], 
                           Omega_k=0.02, w0=-1., wa=0., A_s=2.1e-9, n_s=0.96, 
                           m_nu=0.06)
    P.plot(z, ccl.comoving_angular_distance(cosmo, 1./(1.+z)) / DM0 - 1., 'm-')
    """
    
    # Physical matter fraction
    cmb_om = 0.1386
    cmb_om_err = 0.019 * 0.1386 # 1.9%
    this_om = pmodel['omegaM'] * pmodel['h']**2.
    P.text(30., -0.1, "$\Delta\omega_m/\sigma_{\omega_m} = %+3.3f$"
                        % ((this_om - cmb_om)/cmb_om_err), fontsize=13)
    
    P.xlabel("$z$", fontsize=16.)
    P.ylabel(r"$\Delta D / D$", fontsize=16.)

    #P.xlim((0.08, 3.))
    P.xlim((0.08, 1500.))
    #P.ylim((-0.5, 30.))
    P.xscale('log')

    P.legend(loc='upper right', frameon=False, ncol=2)

    P.tight_layout()
    P.show()


if __name__ == '__main__':
    plot_lss_deviations()
    
