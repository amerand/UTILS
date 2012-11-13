"""
Bayesian Inference first experimentation
"""

import numpy as np
from scipy import special
from matplotlib import pyplot
import dpfit

def normDistr(x,mu,s):
    """
    normal distribution, or probability density to have x +/- s
    knowing the model predicts mu
    """
    return 1./(s*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*s**2))

def v2(bd_l):
    """
    simple UD v2 function. 
    bd_l in m*mas/um
    """
    c = np.pi*np.pi/(180*3600*1000.0)*1e6

    bd_l = np.array(c*bd_l)
    bd_l += (bd_l==0)*1e-6
    return (2*special.jv(1,bd_l)/bd_l)**2

def dpV2(b_l, p):
    c = np.pi*np.pi/(180*3600*1000.0)*1e6
    bd_l = np.array(c*b_l*p['diam'])
    bd_l += (bd_l==0)*1e-6
    return (2*special.jv(1,bd_l)/bd_l)**2

def test1(N=10):
    diamMin = 1; diamMax = 2; P=1000
    bMin = 50; bMax = 150 
    err = 0.01 # fractional error
    
    # actual diameters
    diam = diamMin + np.random.rand()*(diamMax-diamMin)

    # observations:
    B = bMin + np.random.rand(N)*(bMax-bMin)
    obs = v2(B*diam/2.2)

    # error bars
    obs += np.random.randn(N)*err

    
    # prior:
    diams = np.linspace(diamMin, diamMax, P)
    prior = np.ones(P)/float(P)

    # inference
    for i,o in enumerate(obs):
        newPrior = []
        probs = np.array([normDistr(o, v2(B[i]*d/2.2), err*o)
                          for d in diams])
        norma = (prior*probs).sum()
        for k,p in enumerate(prior):
            newPrior.append(probs[k]*prior[k]/norma)
        # update
        prior = np.array(newPrior)

    # plot
    pyplot.clf()
    pyplot.vlines(diam, 0, 1.1*prior.max(), color='k',
                  linewidth=3, alpha=0.2, label='real value')
    pyplot.plot(diams, prior, label='Bayesian')
    
    print 'Least Square:'
    guess = {'diam':0.5*diamMin+0.5*diamMax}
    pfix, uncer, reducedChi2, model = \
          dpfit.leastsqFit(dpV2, B/2.2, guess, obs,
                           err=err*obs)
    print 'CHI2 reduced=', reducedChi2
    pyplot.vlines(pfix['diam'],0,prior.max(), 'r')
    pyplot.hlines(prior.max()/2,
                  pfix['diam']-uncer['diam'],
                  pfix['diam']+uncer['diam'], color='r',
                  label='least square')
    pyplot.legend()
    return
    
