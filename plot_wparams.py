#!/usr/bin/python
"""
Plot results of MCMC chains.
"""
import numpy as np
import pylab as P
import wztanh as model

#fnames = ["wcdm_cmbonly.dat", "wcdm_cmb_lss.dat", "wcdm_cmb_lss_lya.dat",
#          "wztanh_cmb_lss.dat", "wztanh_zc100_cmb_lss.dat"]
#labels = ['wCDM CMB only', 'wCDM CMB + BAO', 'wCDM CMB + BAO + Lya', 
#          'tanh CMB + BAO', 'tanh CMB + BAO, zc<100']
fnames = ["wztanh_cmb_lss.dat", "wztanh_zc100_cmb_lss.dat", "wztanh_cmb_lss_MADEUP.dat"]
labels = ['tanh CMB + BAO', 'tanh CMB + BAO, zc<100', 'tanh CMB + madeup']




def load_chain(fname):
    """
    Load MCMC chain output from file
    """
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


"""
P.subplot(111)
for i, fname in enumerate(fnames):
    dat = load_chain(fname)
    print dat.keys()
    print dat['w0'].shape
    x = (dat['w0'][10000:] - dat['winf'][10000:]) / (1. + dat['zc'][10000:])
    #dat['w0'][10000:] - dat['winf'][10000:]
    y = dat['logl'][10000:] #dat['zc'][10000:]
    P.plot(x, y, ls='none', marker='.', alpha=0.15, label=labels[i])

P.xlabel("Steppyness $(w_0 - w_\infty) / (1 + z_c)$", fontsize=15)
P.ylabel("$\log\mathcal{L}$", fontsize=15)

#P.plot(0.67, 0.318, 'ko')
#P.plot(-1., 0.318, 'ko')
P.axvline(0., color='k', ls='dashed')

#P.yscale('log')
P.legend(loc='upper left', frameon=False)
P.tight_layout()

P.show()
exit()
"""


#fnames = ["wztanh_cmb_lss_MADEUP.dat",]

# Load data and select random samples
#dat = load_chain(fnames[0])
dat = load_chain("wzmocker_cmb_lss.dat")
np.random.seed(10)
print dat.keys()

idxs = np.random.randint(low=10000, high=dat['h'].size, size=1500)

#idxs = np.where(np.logical_and(dat['logl'] > -3., 
#                              (dat['w0'] - dat['winf'])/(1. + dat['zc']) > -0.8))
#idxs = idxs[0]

#idxs = np.where(dat['logl'] > -3.2)
#idxs = idxs[0]
#print idxs.size

# Scale factor array
z = np.linspace(0., 10., 500)
a = 1./(1.+z)

P.subplot(111)
#P.text(6., -0.25, "$\log\mathcal{L}>-3$,\n$(w_0 - w_\infty)/(1+z_c) > 0.2$")
#P.text(6., -0.25, "$\log\mathcal{L}>-4$")

for i in idxs:
    p = {key: dat[key][i] for key in dat.keys()}
    p['mocker'] = True
    P.plot(z, model.wz(a, p), 'b-', lw=1.8, alpha=0.1)

P.ylim((-2., 0.))

P.xlabel("z", fontsize=15)
P.ylabel("w(z)", fontsize=15)

#P.yscale('log')
P.tight_layout()
P.show()

"""
# Load data and select random samples
fnames = ["wztanh_cmb_lss_MADEUP.dat",]
dat = load_chain(fnames[0])

z = np.linspace(0., 10., 100)
a = 1./(1.+z)

# Get percentiles of w(z) in each z bin
wz = [model.wz(a, {key: dat[key][i] for key in dat.keys()}) 
      for i in range(10000, dat['h'].size)] # dat['h'].size
pcts = [2.5, 16., 50., 84., 97.5]
print dat['h'].size

pct = np.percentile(wz, pcts, axis=0)
for i in range(pct.shape[0]):
    P.plot(z, pct[i], 'k-', lw=1.8, alpha=1. - np.abs(pcts[i] - 50.)/70.)

P.plot(z, np.mean(wz, axis=0), 'r-', lw=1.8)

P.ylim((-2., 0.))

P.xlabel("z", fontsize=15)
P.ylabel("w(z)", fontsize=15)

P.tight_layout()
P.show()
"""
