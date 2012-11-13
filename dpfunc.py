"""
collection of functions parametrized by dictionnaries, to allow easy fitting
with dpfit.py

author: amerand@eso.org
"""

import numpy as np

def polyN(x, params):
    """
    Polynomial function. e.g. params={'A0':1.0, 'A2':2.0} returns
    x->1+2*x**2. The coefficients do not have to start at 0 and do not
    need to be continuous. Actually, the way the function is written,
    it will accept any 'Ax' where x i a float.
    """
    res = 0
    for k in params.keys():
        res += params[k]*np.array(x)**float(k[1:])
    return res

def gaussian(x,params):
    """
    1D gaussian function of moment 'MU' and variance 'SIGMA'**2
    
    params: {'MU':, 'SIGMA':, ('OFFSET':,) ('AMP':)} 'AMP' and
    'OFFSET' are optional: if AMP is not given, the amplitude is set
    to 1/(sigma*sqrt(2*pi)).
    """
    res = np.exp(-(x-params['MU'])**2/(2*params['SIGMA']**2))
    if params.has_key('AMP'):
        res *= params['AMP']
    else:
        res *= 1./(params['SIGMA']*np.sqrt(2*np.pi))
        
    if params.has_key('OFFSET'):
        res += params['OFFSET']
    return res

def lorentzian(x,params):
    """
    1D lorentzian function of moment 'MU' and parameter 'GAMMA'
    
    params: {'MU':, 'GAMMA':, ('OFFSET':,) ('AMP':)}. 'AMP' and
    'OFFSET' are optional: if AMP is not given, the amplitude is set to
    1/(sigma*sqrt(2*pi)). """
    res = 1/(1+(x-params['MU'])**2/(params['GAMMA']**2))
    if params.has_key('AMP'):
        res *= params['AMP']
    else:
        res *= 1./(params['GAMMA']*np.pi)
    if params.has_key('OFFSET'):
        res += params['OFFSET']
    return res

def sin(x, params):
    """
    sinusoidal wave
    params: {'AMP':amplitude, 'WAV':wavelength, ('OFFSET':offset), ('PHI':phase in radian)}
    
    returns AMP*sin(2*pi*x/WAV + PHI) + OFFSET
    """
    xx = 2*np.pi*x/params['WAV']
    if params.has_key('PHI'):
        xx += params['PHI']
    res = params['AMP']*np.sin(xx)
    if params.has_key('OFFSET'):
        res += params['OFFSET']
    return res

def cos(x, params):
    """
    cosinusoidal wave
    params: {'AMP':amplitude, 'WAV':wavelength, ('OFFSET':offset), ('PHI':phase in radian)}
    
    returns AMP*cos(2*pi*x/WAV + PHI) + OFFSET
    """
    xx = 2*np.pi*x/params['WAV']
    if params.has_key('PHI'):
        xx += params['PHI']
    res = params['AMP']*np.cos(xx)
    if params.has_key('OFFSET'):
        res += params['OFFSET']
    return res
    
def fourier(x, params):
    """
    fourier serie:
    (A0) + sum(Ak*sin(k*2*pi*x/WAV +PHIk))_k=1...
    
    A0 is optional. Parameters should contain at least WAV and one (Ak,PHIk)   
    """
    phi = filter(lambda k: k[:3]=='PHI', params.keys())
    res = np.copy(x)
    res *=0
    for f in phi:
        res += params['A'+f[3:]]*np.sin(float(f[3:])*2*np.pi*x/params['WAV']
                                                + params[f])
    if params.has_key('A0'):
        res += params['A0']
    return res

