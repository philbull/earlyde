#!/usr/bin/env python
"""
Plot results of MCMC chains.
"""
import numpy as np
import pylab as P
import wztanh as model
from load_data_files import load_chain

#MODE = 'scatter'
#MODE = 'wz'
MODE = 'percentiles'
#MODE = 'hist'

SAMPLE_SEED = 88 #44
BURNIN = 1500000 #1200000
NSAMP = 5000
np.random.seed(SAMPLE_SEED)

#fnames = ["wcdm_cmbonly.dat", "wcdm_cmb_lss.dat", "wcdm_cmb_lss_lya.dat",
#          "wztanh_cmb_lss.dat", "wztanh_zc100_cmb_lss.dat"]
#labels = ['wCDM CMB only', 'wCDM CMB + BAO', 'wCDM CMB + BAO + Lya', 
#          'tanh CMB + BAO', 'tanh CMB + BAO, zc<100']
#fnames = ["wztanh_cmb_lss.dat", "wztanh_zc100_cmb_lss.dat", "wztanh_cmb_lss_MADEUP.dat"]
#labels = ['tanh CMB + BAO', 'tanh CMB + BAO, zc<100', 'tanh CMB + madeup']

#fnames = ["wztanh_cmb_lss.dat", "wztanh_cmb_lss_desi.dat"]
#labels = ['tanh CMB + LSS', 'tanh CMB + LSS + DESI']

"""
fnames = [ "new_wztanh_cmb_lss.dat", 
           "new_wztanh_cmb_lss_desi.dat", 
           "new_wztanh_cmb_lss_desi_cosvis.dat",
           "new_wztanh_cmb_lss_desi_cosvis_nowedge.dat" ]
labels = [ 'tanh CMB + LSS', 
           'tanh CMB + LSS + DESI', 
           'tanh CMB + LSS + DESI + Stage 2',
           'tanh CMB + LSS + DESI + Stage 2 no wedge' ]
"""

#fnames = [ "final2_wztanh_cmb_lss.dat", 
#           "final_wztanh_cmb_lss_desi.dat", 
#           "final2_wztanh_cmb_lss_hirax_hw.dat",
#           "final2_wztanh_cmb_lss_hirax_nw.dat" ]
#labels = [ 'tanh CMB + LSS', 
#           'tanh CMB + LSS + DESI', 
#           'tanh CMB + LSS + HIRAX (horiz. wedge)',
#           'tanh CMB + LSS + HIRAX (no wedge)' ]

fnames = [ "chains/final_wztanh_seed30_cmb_lss_desi", ]
labels = [ "DESI" ]


#def load_chain(fname):
#    """
#    Load MCMC chain output from file
#    """
#    # Load header
#    f = open(fname, 'r')
#    hdr = f.readline()
#    f.close()
#    names = hdr[2:-1].split(' ')
#    
#    # Load data
#    dat = np.genfromtxt(fname).T
#    data = {}
#    for i, n in enumerate(names):
#        data[n] = dat[i]
#    return data


if MODE == 'scatter':
    P.subplot(111)
    for i, fname in enumerate(fnames):
        dat = load_chain(fname)
        print(dat.keys())
        print(dat['w0'].shape)
        
        # Select random subsample
        sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
        
        # Steppyness
        #x = (dat['w0'][10000:] - dat['winf'][10000:]) / (1. + dat['zc'][10000:])
        #y = dat['logl'][10000:] #dat['zc'][10000:] # logL
        
        # w0, winf
        x = dat['w0'][sample_idxs] # FIXME: 100k!
        y = dat['winf'][sample_idxs]
        
        P.plot(x, y, ls='none', marker='.', alpha=0.15, label=labels[i])

    #P.xlabel("Steppyness $(w_0 - w_\infty) / (1 + z_c)$", fontsize=15)
    #P.ylabel("$\log\mathcal{L}$", fontsize=15)
    P.xlabel("$w_0$", fontsize=15)
    P.ylabel("$w_\infty$", fontsize=15)

    #P.plot(0.67, 0.318, 'ko')
    #P.plot(-1., 0.318, 'ko')
    P.axvline(-1., color='k', ls='dashed')

    #P.yscale('log')
    P.legend(loc='upper left', frameon=False)
    P.tight_layout()

    P.show()


if MODE == 'hist':
    P.subplot(111)
    for i, fname in enumerate(fnames):
        dat = load_chain(fname)
        # Select random subsample
        sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)
        
        print(dat.keys())
        print(dat['w0'].shape)
        x = (dat['w0'][sample_idxs] - dat['winf'][sample_idxs]) \
          / (1. + dat['zc'][sample_idxs])
        #dat['w0'][10000:] - dat['winf'][10000:]
        #y = dat['logl'][10000:] #dat['zc'][10000:]
        print("%s: std = %3.3f" % (labels[i], np.std(x)))
        P.hist(x, bins=100, range=(-0.2, 0.2), histtype='step', label=labels[i])

    P.xlabel("Steppyness $(w_0 - w_\infty) / (1 + z_c)$", fontsize=15)
    #P.ylabel("$\log\mathcal{L}$", fontsize=15)

    #P.plot(0.67, 0.318, 'ko')
    #P.plot(-1., 0.318, 'ko')
    P.axvline(0., color='k', ls='dashed')

    #P.yscale('log')
    P.legend(loc='upper left', frameon=False)
    P.tight_layout()

    P.show()


if MODE == 'wz':
    # Load data and select random samples
    #dat = load_chain("wztanh_cmb_lss.dat")
    #dat = load_chain("wztanh_cmb_lss_desi.dat")
    dat = load_chain(fnames[0])
    
    # Select random subsample
    sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)

    #idxs = np.where(np.logical_and(dat['logl'] > -3., 
    #                              (dat['w0'] - dat['winf'])/(1. + dat['zc']) > -0.8))
    #idxs = idxs[0]

    #idxs = np.where(dat['logl'] > -3.2)
    #idxs = idxs[0]
    #print idxs.size

    # Scale factor array
    z = np.logspace(-4., 3., 500)
    a = 1./(1.+z)

    P.subplot(111)
    #P.text(6., -0.25, "$\log\mathcal{L}>-3$,\n$(w_0 - w_\infty)/(1+z_c) > 0.2$")
    #P.text(6., -0.25, "$\log\mathcal{L}>-4$")

    for j, i in enumerate(sample_idxs):
        if j % 100 == 0: print("Calculating sample %d / %d" % (j, sample_idxs.size))
        p = {key: dat[key][i] for key in dat.keys()}
        #p['mocker'] = True
        P.plot(z, model.wz(a, p), 'b-', lw=1.8, alpha=0.1)
        omegaK = 0.
        OmegaDE0 = 1. - p['omegaM'] - model.omegaR(p) - omegaK
        #P.plot(z, OmegaDE0*model.omegaDE(a, p), 'b-', lw=1.8, alpha=0.1)

    #P.ylim((-2., 0.))
    
    plcdm = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
                        w0=-1., winf=-1., zc=1e5, deltaz=0.5)
    OmegaDE0 = 1. - plcdm['omegaM'] - model.omegaR(plcdm) - omegaK
    #P.plot(z, OmegaDE0*model.omegaDE(a, plcdm), 'k-', lw=1.8, alpha=1.)
    P.plot(z, model.wz(a, plcdm), 'k-', lw=1.8, alpha=1.)
    
    P.xlabel("z", fontsize=15)
    P.ylabel("w(z)", fontsize=15)
    
    P.ylim((-2., 0.))

    #P.yscale('log')
    P.xscale('log')
    P.tight_layout()
    P.show()



if MODE == 'percentiles':
    # Load data and select random samples
    #fnames = ["wztanh_cmb_lss.dat", "wztanh_cmb_lss_desi.dat",
    #          "wzmocker_cmb_lss.dat", "wzmocker_cmb_lss_desi.dat",]
    fn = fnames[0]
    dat = load_chain(fn)
    
    # Select random subsample
    sample_idxs = np.random.randint(low=BURNIN, high=dat['h'].size, size=NSAMP)

    z = np.linspace(0., 10., 100)
    a = 1./(1.+z)
    
    ode = []; wz = []
    for j, i in enumerate(sample_idxs):
        if j % 100 == 0: print("Calculating sample %d / %d" % (j, sample_idxs.size))
        pp = {key: dat[key][i] for key in dat.keys()}
        #pp['mocker'] = True # FIXME
        wz.append( model.wz(a, pp) )
        #OmegaDE0 = 1. - pp['omegaM'] - model.omegaR(pp) - 0.
        #ode.append( OmegaDE0 * model.omegaDE(a, pp) )
    
    # Get percentiles of w(z) in each z bin
    #wz = [model.wz(a, {key: dat[key][i] for key in dat.keys()}) 
    #      for i in range(10000, dat['h'].size)] # dat['h'].size
    pcts = [2.5, 16., 50., 84., 97.5]
    print(dat['h'].size)

    #pct = np.percentile(ode, pcts, axis=0)
    pct = np.percentile(wz, pcts, axis=0)
    for i in range(pct.shape[0]):
        P.plot(z, pct[i], 'k-', lw=1.8, alpha=1. - np.abs(pcts[i] - 50.)/70.)

    #P.plot(z, np.mean(ode, axis=0), 'r-', lw=1.8)
    P.plot(z, np.mean(wz, axis=0), 'r-', lw=1.8)
    
    # LCDM curves
    plcdm = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
                        w0=-1., winf=-1., zc=1e5, deltaz=0.5)
    #OmegaDE0 = 1. - plcdm['omegaM'] - model.omegaR(plcdm) - 0.
    #P.plot(z, OmegaDE0*model.omegaDE(a, plcdm), 'y-', lw=1.8, alpha=1.)
    
    #P.plot(z, plcdm['omegaM'] * a**-3., 'g-', lw=1.8)
    
    #lbls = ["z", "OmegaDE_LCDM",]
    #lbls += ["OmegaDE[%2.1f pct]" % _pct for _pct in pcts]
    #np.savetxt("%s.percentiles" % fn,
    #           np.column_stack((z, OmegaDE0*model.omegaDE(a, plcdm), pct.T)), 
    #           header=" ".join(lbls))
    
    #P.ylim((-2., 0.))
    #P.xlim((0., 10.))
    #P.ylim((0., 5.))
    P.ylim((-2., 1.))

    P.xlabel("$z$", fontsize=15)
    #P.ylabel(r"$\Omega_{\rm DE}(z)$", fontsize=15)
    P.ylabel(r"$w(z)$", fontsize=15)

    P.tight_layout()
    P.show()

