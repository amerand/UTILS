import scipy.optimize
import numpy as np
import time
import dpfunc
from dpfunc import *
"""
IDEA: fit Y = F(X,A) where A is a dictionnary describing the
parameters of the function.

note that the items in the dictionnary should all be scalar!

author: amerand@eso.org
"""

verboseTime=time.time()

def example():
    """
    very simple example
    """
    X = [ 0,   1,   2,   3  ]
    Y = [-0.1, 1.1, 4.1, 8.9]
    E = [ 0.1, 0.1, 0.1, 0.1]
    #best, unc, chi2, model =\
    #      leastsqFit(dpfunc.polyN, X, 
    #                 {'A0':0., 'A1':0.,'A2':0.1},
    #                 Y, err=E, fitOnly=['A2', 'A0'])
    fit=leastsqFit(dpfunc.polyN, X, 
                     {'A0':0., 'A1':0.,'A2':0.1},
                     Y, err=E, doNotFit=['A1'])
    
    print 'CHI2=', fit['chi2']
    for k in fit['best'].keys():
        print k, '=', fit['best'][k],
        if fit['uncer'][k]>0:
            print '+/-', fit['uncer'][k]
        else:
            print ''
    print 'Y=', Y
    print 'MODEL=', fit['model']
    return

def meta(x, params):
    """
    allows to call any combination of function defines inside dpfunc:
    
    params={'funcA;1:p1':, 'funcA;1:p2':,
            'funcA;2:p1':, 'funcA;2:p2':,
            'funcB:p1':, etc}
    
    funcA and funcB should be defined in dpfunc.py. Allows to call many
    instances of the same function (here funcA) and combine different functions.
    Outputs of the difference functions will be sumed usinf operator '+'. """
    
    # -- list of functions:
    funcs = set([k.strip().split(':')[0].strip() for k in params.keys()])
    #print funcs

    res = 0
    for f in funcs: # for each function
        # -- keep only relevant keywords
        kz = filter(lambda k: k.strip().split(':')[0].strip()==f, params.keys())
        tmp = {}
        for k in kz:
            # -- build temporary dict pf parameters
            tmp[k.split(':')[1].strip()]=params[k]
        ff = f.split(';')[0].strip() # actual function name
        if not dpfunc.__dict__.has_key(ff):
            raise NameError(ff+' not defined in dpfunc')
        # -- add to result the function result
        res += dpfunc.__dict__[ff](x, tmp)
    return res

def leastsqFit(func, x, params, y, err=None, fitOnly=None,
               verbose=False, doNotFit=[], epsfcn=1e-7,
               ftol=1e-5, fullOutput=True, normalizedUncer=True):
    """
    - params is a Dict containing the first guess.

    - fits 'y +- err = func(x,params)'. errors are optionnal.

    - fitOnly is a LIST of keywords to fit. By default, it fits all
      parameters in 'params'. Alternatively, one can give a list of
      parameters not to be fitted, as 'doNotFit='

    - doNotFit has a similar purpose: for example if params={'a0':,
      'a1': 'b1':, 'b2':}, doNotFit=['a'] will result in fitting only
    the 'b1' and 'b2'. WARNING: if you name parameter 'A' and another one 'AA',
    you cannot use doNotFit to exclude only 'A' since 'AA' will be excluded as
    well... 
    
    returns bestparam, uncertainties, chi2_reduced, func(x, bestparam)
    """
    # fit all parameters by default
    if fitOnly is None:
        if len(doNotFit)>0:
            fitOnly = filter(lambda x: x not in doNotFit, params.keys())
        else:
            fitOnly = params.keys()

    # build fitted parameters vector:
    pfit = [params[k] for k in fitOnly]

    # built fixed parameters dict:
    pfix = {}
    for k in params.keys():
        if k not in fitOnly:
            pfix[k]=params[k]

    if verbose:
        print '[dpfit] %d FITTED parameters:'%(len(fitOnly)), fitOnly
        
    # actual fit
    plsq, cov, info, mesg, ier = \
              scipy.optimize.leastsq(fitFunc, pfit,
                    args=(fitOnly,x,y,err,func,pfix, verbose),
                    full_output=True, epsfcn=epsfcn, ftol=ftol)
    
    # best fit -> agregate to pfix
    for i,k in enumerate(fitOnly):
        pfix[k] = plsq[i]

    # reduced chi2
    model = func(x,pfix)
    tmp = fitFunc(plsq, fitOnly, x, y, err, func, pfix)
    try:
        chi2 = (np.array(tmp)**2).sum()
    except:
        chi2=0
        for x in tmp:
            chi2+=np.sum(x**2)
    reducedChi2 = chi2/float(reduce(lambda x,y: x+y,
                  [1 if np.isscalar(i) else len(i) for i in y])-len(pfit)+1)

    # uncertainties:
    uncer = {}
    for k in pfix.keys():
        if not k in fitOnly:
            uncer[k]=0 # not fitted, uncertatinties to 0
        else:
            i = fitOnly.index(k)
            if cov is None:
                uncer[k]=-1
            else:
                uncer[k]= np.sqrt(np.abs(np.diag(cov)[i]))
                if normalizedUncer:
                    uncer[k]*= np.sqrt(reducedChi2)

    if verbose:
        print '-'*20
        print 'REDUCED CHI2=', reducedChi2
        tmp = pfix.keys(); tmp.sort()
        for k in tmp:
            print k, '=', pfix[k],
            if uncer[k]!=0:
                print '+/-', uncer[k]
            else:
                print ''
    # result:
    if fullOutput:
        pfix={'best':pfix, 'uncer':uncer,
              'chi2':reducedChi2, 'model':model}
    return pfix

def fitFunc(pfit, pfitKeys, x, y, err=None, func=None,
            pfix=None, verbose=False):
    """
    interface to leastsq from scipy:
    - x,y,err are the data
    - pfit is a list of the paramters
    - pfitsKeys are the keys to build the dict
    pfit and pfix (optional) and combines the two
    in 'A', in order to call F(X,A)
    """
    global verboseTime
    params = {}
    # build dic from parameters to fit and their values:
    for i,k in enumerate(pfitKeys):
        params[k]=pfit[i]
    # complete with the non fitted parameters:
    for k in pfix:
        params[k]=pfix[k]
    if err is None:
        err = np.ones(np.array(y).shape)
    # return residuals
    try:
        # assumes y is a numpy array
        y = np.array(y)
        res= ((func(x,params)-y)/err).flatten()
    except:
        # much slower: this time assumes y (and the result from func) is
        # a list of stuff, each convertible in np.array
        res = []
        tmp = func(x,params)
        for k in range(len(y)):
            df = (np.array(tmp[k])-np.array(y[k]))/np.array(err[k])
            try:
                res.extend(list(df))
            except:
                res.append(df)
        res= np.array(res)
        
    if verbose and time.time()>(verboseTime+1):
        verboseTime = time.time()
        print time.asctime(),
        try:
            chi2=(res**2).sum/(len(res)-len(pfit)+1.0)
            print 'CHI2:', chi2
        except:
            # list of elements
            chi2 = 0
            N = 0
            res2 = []
            for r in res:
                try:
                    chi2 += np.sum(r**2)
                    N += len(r)
                    res2.extend(list(r))
                except:
                    chi2 += r**2
                    N += 1
                    res2.append(r)   
            res = res2
            print 'CHI2:', chi2/float(N-len(pfit)+1)
    return  res
        

