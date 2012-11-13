import pyfits
import numpy as np
import scipy.signal
from matplotlib import pyplot
import dpfit # not a standard Python library!
import os
import time

# -- mouse button+key to fit
keyCont='q'
keyLorent='z'
keyGauss='a'

class Open():
    """
    read spectrum from FITS files (works for Elodie at least). Fit spectral
    lines and save the results in a binary table. FITS files having the same
    binary table will be reckognized.
    """
    def __init__(self, filename, chi2threshold=1.0):
        """
        open file and read it
        """
        self.chi2threshold = chi2threshold
        self.filename = filename
        # -- check if file exist
        if not os.path.isfile(self.filename):
            raise NameError(self.filename+' does not exist')
        
        self.fits = pyfits.open(filename) # open FITS file
        # -- extract spectrum: data of the first HDU
        sp = self.fits[0].data 

        # -- reconstruct wavelength table
        wl = (np.arange(1,len(sp)+1)-self.fits[0].header['CRPIX1'])*\
             self.fits[0].header['CDELT1']+self.fits[0].header['CRVAL1']

        # -- try to guess the unit in meters:
        self.unitStr = self.fits[0].header['CUNIT1'].strip() # from header
        unitDic = {'nm':1e-9, 'um':1e-6, 'm':1, '0.1 nm':1e-10,
                   'A':1e-10, 'cm':1e-2}
        self.unit = unitDic[self.unitStr] # in meters
        print 'WAVELENGTHS: ', '"'+self.unitStr+'"',\
            unichr(0x21d2), self.unit, 'm'
        
        # -- keep only finite values in spectrum:
        wl, sp = wl[np.isfinite(sp)], sp[np.isfinite(sp)]
        
        # -- set public variables
        self.spectrum = sp # raw spectrum
        self.noise = (sp-scipy.signal.medfilt(self.spectrum, 5)).std()
        self.wlTable = wl # wavelength table
        self.model={'polyN:A0':sp.mean()} # stores the model
        self.whereCont=[] # user defined continuum 

    def __fitContinuum__(self):
        """
        fit polynomial continuum and adjust the model (adds parameters)
        """
        fitOnly = filter(lambda k: k[:5]=='polyN', self.model.keys())
        #print self.model, fitOnly
        best, uncer, reducedChi2, model = dpfit.leastsqFit(
            dpfit.meta, self.wlTable[self.whereCont],
            self.model, self.spectrum[self.whereCont], err=self.noise,
            fitOnly=fitOnly, verbose=False)
        #print 'Continuum:  CHI2=', reducedChi2
        while reducedChi2>self.chi2threshold:
            N = len(fitOnly)
            best['polyN:A'+str(N)] = 0.0
            fitOnly = filter(lambda k: k[:5]=='polyN', best.keys())
            best, uncer, reducedChi2, model = dpfit.leastsqFit(
                dpfit.meta, self.wlTable[self.whereCont],
                best, self.spectrum[self.whereCont], err=self.noise,
                fitOnly=fitOnly, verbose=False)
            #print '  continuum -> adding order', N, 'CHI2=', reducedChi2
        self.model = best
        return

    def __fitAll__(self, noCont=True, restrict=True):
        if noCont:
            doNotFit=filter(lambda k: 'polyN' in k, self.model.keys())
        else:
            doNotFit=[]
    
        if restrict:
            w = np.zeros(len(self.spectrum))
            for k in self.model.keys():
                if 'SIGMA' in k or 'GAMMA' in k:
                    w += (np.abs(self.wlTable-self.model[k.split(':')[0]+':MU'])<=
                          np.abs(self.model[k]))
            w = np.where(w)
        else:
            w = np.where(np.isfinite(self.spectrum))
        
        best, uncer, reducedChi2, model = dpfit.leastsqFit(
            dpfit.meta, self.wlTable[w], self.model, self.spectrum[w],
            err=self.noise, verbose=False, doNotFit=doNotFit)
        self.model = best

        return
        
    def __replot__(self):
        xlim = pyplot.xlim()
        ylim = pyplot.ylim()
        pyplot.clf()
        # -- spectrum
        pyplot.plot(self.wlTable, self.spectrum, 'k')
        # -- continuum
        pyplot.plot(self.wlTable[self.whereCont],
                    self.spectrum[self.whereCont],
                    '.', color='orange', alpha=0.5)
        # -- model
        #print self.model
        pyplot.plot(self.wlTable,
                    dpfit.meta(self.wlTable, self.model),
                    color='red', linewidth=3, alpha=0.5)
        pyplot.xlim(xlim[0], xlim[1])
        pyplot.ylim(ylim[0], ylim[1])
        return

    def __handler__(self, event):
        """
        React to a click of the mouse: will fit a line (Lorentzian) on
        the line. Private method, not to be called by the user.
        """
        if event.name=='button_press_event' and \
               event.button==1:
            # -- set first limit for fit
            self.wlMIN = event.xdata
           
        if event.name=='button_release_event' and \
               event.button==1: 
            # -- set second limit
            self.wlMAX = event.xdata
            # -- invert min and max if needed
            try:
                if self.wlMAX<self.wlMIN:
                    self.wlMIN, self.wlMAX = self.wlMAX, self.wlMIN
            except:
                pass
            if event.key==keyCont:
                # -- add range to continuum 
                self.whereCont.extend(np.where((self.wlTable<=self.wlMAX)*
                    (self.wlTable>=self.wlMIN))[0])
                self.__fitContinuum__()
                self.__replot__()
            elif event.key==keyGauss:
                # --    
                g = set(filter(lambda k: k.split(';')[0]=='gaussian',
                               self.model.keys()))
                N = str(len(g)/3+1)
                if self.wlMIN == self.wlMAX:
                    cont = dpfit.meta(event.xdata, self.model)
                    self.model['gaussian;'+N+':MU']=event.xdata
                    self.model['gaussian;'+N+':AMP']=event.ydata-cont
                    self.model['gaussian;'+N+':SIGMA']=1e-3*event.xdata
                else:
                    w = np.where((self.wlTable>=self.wlMIN)*
                                (self.wlTable<=self.wlMAX))
                    cont = dpfit.meta(self.wlTable[w].mean(), self.model)
                    self.model['gaussian;'+N+':MU']=self.wlTable[w].mean()
                    self.model['gaussian;'+N+':AMP']=self.spectrum[w].min()-cont
                    self.model['gaussian;'+N+':SIGMA']=(self.wlMAX-self.wlMIN)/2   
                self.__fitAll__()
                self.__replot__()
            elif event.key==keyLorent:
                # --    
                g = set(filter(lambda k: k.split(';')[0]=='lorentzian',
                               self.model.keys()))
                N = str(len(g)/3+1)
                if self.wlMIN == self.wlMAX:
                    cont = dpfit.meta(event.xdata, self.model)
                    self.model['lorentzian;'+N+':MU']=event.xdata
                    self.model['lorentzian;'+N+':AMP']=event.ydata-cont
                    self.model['lorentzian;'+N+':GAMMA']=1e-3*event.xdata
                else:
                    w = np.where((self.wlTable>=self.wlMIN)*
                                (self.wlTable<=self.wlMAX))
                    cont = dpfit.meta(self.wlTable[w].mean(), self.model)
                    self.model['lorentzian;'+N+':MU']=self.wlTable[w].mean()
                    self.model['lorentzian;'+N+':AMP']=self.spectrum[w].min()-cont
                    self.model['lorentzian;'+N+':GAMMA']=(self.wlMAX-self.wlMIN)/2
                self.__fitAll__()
                self.__replot__()
        return
    
    def Plot(self, win=0):
        """
        Initialize the plotting window and plot the spectrum.
        """
        # -- init figure
        pyplot.close(win)
        self.fig = pyplot.figure(win)
        #pyplot.title('draw a box with %s click%s%s to fit a line'%
        #             (clickToFit, ' +'if not keyToFit is None else '',
        #              keyToFit))
        pyplot.xlabel(r'wavelength in vacuum ('+ self.unitStr+')')
        # -- plot spectrum
        pyplot.plot(self.wlTable, self.spectrum, 'k')
            
        # -- tell matplotlib what to do in case of mouse click
        self.fig.canvas.mpl_connect('button_press_event', self.__handler__)
        self.fig.canvas.mpl_connect('button_release_event', self.__handler__)
        return
   
    def reduceContinuumOrder(self, n=1, refit=True):
        """
        reduce continuum's polynomial order by n (default 1), redo the fit and
        replot.
        """
        for k in range(1):
            fitOnly = filter(lambda k: k[:5]=='polyN', self.model.keys())
            N = len(fitOnly)
            self.model.pop('polyN:A'+str(N-1))
        if refit:
            self.__fitContinuum__()
        self.__replot__()
        return
        
    def Save(self, filename=None, overwrite=True):
        """
        save files with the fitted line in a FITS file. if no name
        given, will come up with a default one. 'overwrite=True'
        (default) will overwrite the file if it already exists.
        """
        # -- generate default filename
        if filename is None:
            spl = os.path.split(self.filename) # (dir, filename)
            # -- prefer os.path.join(dir, filename)
            # -- to dir+'/'+filename
            filename = os.path.join(spl[0],
                                    os.path.splitext(spl[1])[0]+'_RED'+
                                    os.path.splitext(spl[1])[1])
        
        # -- check if file already exists
        if os.path.isfile(filename):
            if overwrite:
                print 'file', filename, 'already exists: overwriting.'
                os.remove(filename)
            else:
                raise NameError(filename+
                    ' already exists: use another one of use overwrite=True')
        print 'saving to', filename

        self.fits.writeto(filename) # write new FITS file
        return