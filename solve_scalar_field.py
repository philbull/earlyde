#!/usr/bin/env python
"""
Solve Klein-Gordon equation for quintessence scalar field.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import pylab as plt

# Critical density (PDG astrophysical constants)
RHO_CRIT = 2.77536627e11 # h^2 Msun/Mpc^3

params = {
    'h':        0.67,
    'omega_m':  0.33,
    'omega_k':  0.,
    'U0':       1.,     # Potential amplitude, V_0 in units of rho_cr,0
    'lambda':   0.1,    # Scalar potential exponent, dimensionless
    'y0':       0.,     # dphi/dt in units of sqrt(rho_cr,0)
    
}


def U(phi, params):
    """
    Scalar field potential, U(phi) = V(phi) / rho_cr,0.
    """
    return params['U0'] * (np.exp(phi * params['lambda']) + 0.3)

def dUdphi(phi, params):
    """
    Gradient of scalar field potential with respect to scalar field value, 
        
        dU(phi)/dphi ~ dim'less
    """
    #return U(phi, params) / params['lambda']
    return params['U0'] * np.exp(phi * params['lambda']) * params['lambda']


def phi0(params):
    """
    Solve V(phi_0) given dphi/dt|_0 and Omega_phi, i.e. use:
    
        rho_phi(z=0) = (phidot_0^2)/2 + V(phi_0)
        (where all terms are in units of the critical density at z=0)
    
    Returns
    -------
    phi_0 : float
        Value of phi_0 for given initial conditions (phi is dimensionless).
    """
    # Find omega_phi0
    omega_phi0 = 1. - params['omega_m'] - params['omega_k']
    y0 = params['y0']
    
    # Minimise to find value of phi_0 that matches expected scalar field 
    # energy density at z=0 
    min_fn = lambda p: (0.5*y0**2. + U(p, params) - omega_phi0)**2.
    res = minimize_scalar(min_fn, method='brent')
    return res.x


def E(fields, params):
    """
    Dimensionless expansion rate, E(a), defined as:
    
        E^2(a) = Omega_m a^-3 + Omega_K a^-2 + (rho_phi(a) / rho_cr,0)
    """
    a, phi, y = fields
    rho_phi = 0.5*y**2. + U(phi, params) # = rho_phi(a) / rho_cr,0
    
    Esq = params['omega_m']*a**-3. + params['omega_k']*a**-2. + rho_phi
    return np.sqrt(Esq)


def ode_system(tau, fields, params):
    """
    Dimensionless ODE system (from friedmann and Klein-Gordon eqns):
        
        da/dtau = a E(a)
        dphi/dtau = y
        dy/dtau = -3 y E - dU/dphi
    
    We have defined an evolution variable tau = sqrt(rho_cr,0) * t, and the 
    dimensionless quantities:
    
        V(phi) = rho_cr,0 * U(phi)
        dphi/dt = sqrt(rho_cr,0) * y
        rho_phi = rho_cr,0 (y^2 / 2 + U(phi))
    """
    a, phi, y = fields
    Ea = E(fields, params)
    
    # Define ODE right-hand-sides
    da_dtau = a * Ea
    dphi_dtau = y
    dy_dtau = -3.*y*Ea - dUdphi(phi, params)
    
    # Return RHS for ODE
    return [da_dtau, dphi_dtau, dy_dtau]


if __name__ == '__main__':
    
    # Find initial value of phi
    print("phi0 = {:3.3f}".format(phi0(params)))
    
    # Set ODE time range and initial conditions
    tau = np.linspace(0., 1., 1000)[::-1]
    tau_range = (np.max(tau), np.min(tau))
    fields0 = [1., phi0(params), params['y0']] # a_0, phi_0, y_0
    
    # Solve ODE by integrating backwards in evolution parameter, tau
    ode_sys = lambda t, f: ode_system(t, f, params)
    res = solve_ivp(ode_sys, tau_range, fields0, t_eval=tau)
    a, phi, y = res['y'] # Extract field solutions
        
    print("Shape of results:", res['y'].shape)
    
    # Eqn. of state for scalar field
    w = (0.5*y**2. - U(phi, params)) / (0.5*y**2. + U(phi, params))
    
    # Plot solutions as fn. of scale factor
    plt.subplot(221)
    plt.plot(a, -phi, 'b-', lw=1.8, label="$-\phi(a)$")
    plt.plot(a, y, 'r-', lw=1.8, label="$y(a)$")
    plt.plot(a, U(phi, params), 'g--', lw=1.8, label="$U(a)$")
    
    plt.yscale('log')
    plt.xlabel('a')
    plt.legend(loc='upper right')
    
    # Plot potential as fn. of phi
    plt.subplot(222)
    #p = np.linspace(np.min(phi), np.max(phi), 200)
    p = np.linspace(np.min(phi), 0.1, 200)
    plt.plot(p, U(p, params), 'r-', lw=1.8)
    plt.axvline(phi0(params), color='r', ls='dashed', alpha=0.4)
    
    plt.yscale('log')
    plt.ylabel("$U(\phi)$")
    plt.xlabel("$\phi$")
    
    
    # Plot density as fn. of scale factor
    plt.subplot(223)
    
    plt.plot(a, params['omega_m']*a**-3., 'r-', lw=1.8, 
             label=r"$\rho_m / \rho_{\rm cr,0}$")
    plt.plot(a, 0.5*y**2. + U(phi, params), 'k-', lw=1.8, 
             label=r"$\rho_\phi / \rho_{\rm cr,0}$")
    plt.axhline(np.log10(1. - params['omega_m'] - params['omega_k']), 
                color='k', ls='dotted')
    
    plt.xlabel('a')
    plt.yscale('log')
    plt.legend(loc='upper right')
    
    
    # Plot eqn. of state vs. scale factor
    plt.subplot(224)
    plt.plot(a, w, 'y--', lw=1.8)
    plt.ylabel("$w_\phi(a)$")
    plt.xlabel('a')
    
    plt.tight_layout()
    plt.show()
    
    
    
    
    
    
