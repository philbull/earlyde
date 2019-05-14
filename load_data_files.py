import numpy as np
import wztanh as model
import os

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


def load_fisher_data(froot):
    """
    Load Fisher forecast data.
    """
    f_fisher = "%s.dat" % froot
    f_zbins = "%s.zbins" % froot
    
    # Load header
    f = open(f_fisher, 'r')
    hdr = f.readline()
    f.close()
    names = hdr[2:-1].split(' ')
    
    # Load Fisher matrix and bin centre redshifts
    F = np.genfromtxt(f_fisher)
    zbinc = np.genfromtxt(f_zbins)
    
    # Calculate mean values of H and D_A in LCDM model
    abinc = 1. / (1. + zbinc)
    #plcdm = model.pdict(h=0.6731, omegaM=0.315, omegaK=0.0, omegaB=0.045, 
    #                    w0=-1., winf=-1., zc=1e5, deltaz=0.5)
    plcdm = model.pdict(h=0.6727, omegaM=0.3166, omegaK=0.0, omegaB=0.04941, 
                        w0=-1., winf=-1., zc=1e5, deltaz=0.5)
                        
    # Convert units to same as Fisher matrix: H ~ (100 km/s/Mpc); D_A ~ Gpc
    Hc = model.Hz(abinc, plcdm) / 1e2 # 100 km/s/Mpc
    DAc = model.DAz(abinc, plcdm) / 1e3 # Gpc
    
    # Construct mean vector and list of parameter names
    mean_vec = np.concatenate((Hc, DAc))
    lbls = ["H%d" % i for i in range(zbinc.size)] \
         + ["DA%d" % i for i in range(zbinc.size)]
    
    # Invert Fisher matrix to get the covariance
    cov = np.linalg.inv(F)
    
    # Select only certain parameters from the covariance
    new_cov = np.zeros((len(lbls), len(lbls)))
    for i, pni in enumerate(lbls):
        for j, pnj in enumerate(lbls):
            new_cov[i,j] = cov[names.index(pni), names.index(pnj)]
    
    # Return the inverse covariance for this subset of parameters
    icov = np.linalg.inv(new_cov)
    return zbinc, mean_vec, icov


def load_chain(fname, burnin=0):
    """
    Load MCMC chain output from file
    """
    def load_header(hdrname):
        f = open(hdrname, 'r')
        hdr = f.readline()
        f.close()
        names = hdr[2:-1].split(' ')
        return names
    
    # Load from cache if one exists; otherwise, load from dat file
    if os.path.isfile("%s.npy" % fname):
        # Load from cache
        dat = np.load("%s.npy" % fname)
        names = load_header("%s.params" % fname)
        print("Loaded from cache, shape", dat.shape)
    else:
        # Load data from dat file and store in cache
        dat = np.genfromtxt("%s.dat" % fname).T
        np.save("%s" % fname, dat)
        
        # Load parameter names and store in cache
        names = load_header("%s.dat" % fname)
        f = open("%s.params" % fname, 'w')
        f.write("# %s\n" % " ".join(names))
        f.close()
        print("Loaded from dat file, shape", dat.shape)
        
    # Check that data array and list of param names have right shape/length
    assert dat.shape[0] == len(names), \
            "Number of data columns does not match number of column names"
    
    # Reformat data array into dict
    data = {}
    for i, n in enumerate(names):
        data[n] = dat[i][burnin:] # trim burn-in
    
    print(data.keys())
    
    return data
