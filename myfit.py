"""
The idea is to implement a flexible wrapup of minimization found in
scipy.optimize which lacks some features: in particular only fit a
subset of the parameters.

The module also contains some useful functions that can be fitted.

see myfit.example for a simple example how to use the module

"""

import numpy
import scipy
from scipy.optimize import leastsq
from scipy.interpolate import interp1d

def example():
    """
    simple example
    """
    # define data
    x = numpy.linspace(0,2,100)
    y = x + 0.01*numpy.random.randn(len(x))
    err = 0.01*numpy.ones(len(x))

    # first case
    print '-- fit everything --'
    f = fit(PolyN, x, [0,1,0], y, err=err)
    f.leastsqfit() 
    print 'reduced chi2=', f.leastsq_chi2red
    print 'best params =', f.leastsq_best_param
    print 'error bars  =', f.leastsq_errbars
        
    # second case
    print '-- fit only first order --'
    f = fit(PolyN, x, [0,1,0], y, err=err, fit=[1])
    f.leastsqfit()   
    print 'reduced chi2=', f.leastsq_chi2red
    print 'chi2 + 1sigma=', f.chi2red(f.leastsq_best_param[1]+
                                      f.leastsq_errbars[1])
    print 'chi2 - 1sigma=', f.chi2red(f.leastsq_best_param[1] -
                                      f.leastsq_errbars[1])
    print 'best params =', f.leastsq_best_param
    print 'error bars  =', f.leastsq_errbars

    return 

def PolyN(x,a):
    """
    returns the polynomial sum_n(an*x**n)
    a = [a0, a1, a2,...]
    results is a numpy array type
    """
    n = numpy.array(range(len(a)))    
    xx = numpy.array(x).flatten()
    res = numpy.sum(numpy.array(a)[numpy.newaxis,:]*\
                    xx[:,numpy.newaxis]**\
                    n[numpy.newaxis,:], axis=1)
    res = res.reshape(numpy.array(x).shape)
    return res

def Gaussian1D(x,a):
    """
    returns a[0] + a[1]*numpy.exp(-(x-a[2])**2/a[3]**2)
    """
    xx = numpy.array(x).flatten()
    res = a[0] + a[1]*numpy.exp(-(x-a[2])**2/a[3]**2)
    res = res.reshape(numpy.array(x).shape)
    return res

def Gaussian2D(xy,a):
    """
    returns a[0] + a[1]*numpy.exp(-(x'-a[2])**2/a[3]**2 
                                  -(y'-a[4])**2/a[5]**2)
    if a[6] is set, x' and y' (originally xy[0] and xy[1]) will be rotated

    """
    xx = numpy.array(xy).flatten()[:len(numpy.array(xy).ravel())/2]
    yy = numpy.array(xy).flatten()[len(numpy.array(xy).ravel())/2:]
    if len(a)>6:
        xx_ = numpy.cos(a[6])*xx + numpy.sin(a[6])*yy
        yy_ =-numpy.sin(a[6])*xx + numpy.cos(a[6])*yy
        xx = xx_
        yy = yy_
        
    res = a[0] + a[1]*numpy.exp(-(xx-a[2])**2/a[3]**2 -(yy-a[4])**2/a[5]**2)
    if len(numpy.array(xy).shape)>1:
        res = res.reshape(numpy.array(xy[0]).shape)
    return res

def erf(x,a):
    """
    error function for normaly distributed data 
    """
    xx = numpy.array(x).flatten()
    res = a[0] + (1-a[1])*(1+scipy.special.erf((xx-a[2])/a[3]))/2.
    res = res.reshape(numpy.array(x).shape)
    return res

def splineXY(x, a):
    """
    a = [x0,y0, x1,y1, x2,y2...]
    """
    xx = numpy.array(x).flatten()
    res = interp1d(a[0::2], a[1::2], kind='quadratic', \
                   bounds_error=False, fill_value=0.0)(xx)
    res = res.reshape(numpy.array(x).shape)
    return res

def PsplineXY(x,a):
    """
    periodic spline (P=1)

    a = [x0,y0, x1,y1, x2,y2...]
    x in [0,1]
    """
    coef = numpy.array(a)
    xp = numpy.zeros(3*len(coef)/2)
    yp = numpy.zeros(3*len(coef)/2)
    x0 = numpy.mod(coef[::2], 1.0)
    s = numpy.array(x0).argsort()
    xp[0:len(coef)//2]      = numpy.mod(x0[s], 1.0)-1
    xp[len(coef)//2:len(coef)] = xp[0:len(coef)//2]+1
    xp[-len(coef)//2:]      = xp[0:len(coef)//2]+2
    yp[0:len(coef)//2]      = coef[1::2][s]
    yp[len(coef)//2:len(coef)] = yp[0:len(coef)//2]
    yp[-len(coef)//2:]      = yp[0:len(coef)//2]
        
    xx = numpy.array(x).flatten()
    res = interp1d(xp, yp, kind='quadratic', \
                   bounds_error=False, fill_value=0.0)\
                   (numpy.mod(xx, 1))
    res = res.reshape(numpy.array(x).shape)
    return res

def PsplineIntegXY(xt,a):
    """
    periodic spline (P=1)

    a = [c0, c1, x0,y0, x1,y1, x2,y2...]
    t in [0,1]
    x = (t, type)
    type = 0  -> direct
    type = 1 -> c0 + c1*integ
    """ 
    x = xt[0]
    t = xt[1]
    res  = PsplineXY(x,a[2:])
    if isinstance(t, float) or isinstance(t, int):
        w = numpy.where(numpy.array([t])>0)
        x = numpy.array([x])
        res = numpy.array([res])
    else:
        w = numpy.where(t>0)
        
    if len(w[0])>0:
        xx = numpy.linspace(0,1,1000)
        ix = PsplineXY(xx,a[2:])
        ix -= ix.mean()
        ix = numpy.cumsum(ix)
        res[w] = a[0] + a[1]*interp1d(xx,ix)(numpy.mod(x[w],1.0))*\
                 numpy.diff(xx).mean()
    return res

def V_UD(bl, a):
    """
    uniform disk diameter interferometric visibility btl is
    Baseline(m)/wavelength(um), a is UD diameter in mas
    """
    if isinstance(a, float):
        a = [a]
    c = numpy.pi*a[0]*numpy.pi/(180*3600*1000.0)*1e6
    return 2*(scipy.special.jv(1,c*bl)+\
              numpy.float_(bl==0)*1e-6)/\
           (c*bl +numpy.float_(bl==0)*2e-6)

def exp(x, a):
    """
    simple exponentional
    """
    return a[0]+a[1]*numpy.exp(a[2]*(x-a[3]))
    
class fit():
    def __init__(self, func, x, first_guess, data,\
                 err=None, fit=None, **args):
        """    
        optimizes data=func(x, first_guess, **args) 
        
        - 'err' can be a scalar or an array of dims of 'data'        

        - 'fit' is an array of integers defining the parameters to fit
          in 'first_guess'
        """
        self.x    = numpy.array(x)
        self.data = numpy.array(data)
        self.first_guess = numpy.array(first_guess)
        self.n = len(first_guess)
        self.args=args

        # errors
        if err==None:
            self.err = numpy.ones(self.data.shape)
        else:
            self.err = err
                
        # where to fit and not:
        if fit==None:
            self.fit = range(len(first_guess))
        else:         
            tmp_fit = numpy.array(fit)
            tmp_fit.sort()
            self.fit = tmp_fit
            no_fit = []
            for k in range(self.n):
                if not (k==self.fit).max():
                    no_fit.append(k)
            self.no_fit=numpy.array(no_fit)

        self.func = func        
        self.nfree = len(self.x)-len(self.fit)+1
        return
    
    def fit_func(self,param):
        """
        proxy to the actual function. required to handle the fit to a
        subset of parameters
        """
        if len(self.fit)<self.n:
            xparam = numpy.zeros(self.n)
            xparam[self.fit] = param
            xparam[self.no_fit] = self.first_guess[self.no_fit]            
        else:
            xparam = param
        return self.func(self.x,xparam,**self.args)

    def delta(self, param):
        return (self.data-self.fit_func(param))/self.err

    def chi2red(self, param):
        return (self.delta(param)**2).sum()/self.nfree    
    
    def leastsqfit(self, verbose=False, epsfcn=1e-6):
        """
        perform a least square fit
        .leastsq_best_param has the best params (including non fitted)
        .leastsq_errbars has the error bars
        """
        tmp_fg = self.first_guess[self.fit]
        
        tmp_best = leastsq(self.delta, numpy.array(tmp_fg), epsfcn=epsfcn)[0]
        cov = None

        best = numpy.zeros(self.n)
        self.leastsq_errbars = numpy.zeros(self.n)

        ### chi2s
        self.leastsq_chi2 = self.chi2red(tmp_best)*self.nfree
        self.leastsq_chi2red = self.chi2red(tmp_best)        

        ### fitted param
        best[self.fit] = tmp_best
        if not cov is None:
            self.leastsq_errbars[self.fit] = numpy.sqrt(numpy.diag(cov))*\
                                             numpy.sqrt(self.leastsq_chi2)
                                         
        ### unfitted parameters
        if len(self.fit)<self.n:
            best[self.no_fit]=self.first_guess[self.no_fit]
            
        ### store the best solution:
        self.leastsq_best_param = best
                
        return best
