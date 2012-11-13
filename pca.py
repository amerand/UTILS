import numpy

def test():
    """
    """
    from matplotlib import pyplot
    p = 4
    n = 100
    
    x1 = numpy.array([1,-1,1,-1])    
    x2 = numpy.array([1,1, -1,-1])    
    
    sig = 10*x1[numpy.newaxis,:]*numpy.random.rand(n)[:,numpy.newaxis]+\
          1*x2[numpy.newaxis,:]*numpy.random.rand(n)[:,numpy.newaxis]+\
          numpy.random.randn(n,p)
    p = pca(sig)
    print p.base[:,0]
    print p.base[:,1]
    sigf = p.comp(0)
    
    pyplot.figure(0)
    pyplot.clf()
    pyplot.plot(sig[:,0], '+r')
    pyplot.plot(sigf[:,0], 'r-')
    pyplot.plot(sig[:,1], '+b')
    pyplot.plot(sigf[:,1], 'b-')
    pyplot.plot(sig[:,2], '+g')
    pyplot.plot(sigf[:,2], 'g-')
    pyplot.plot(sig[:,3], '+y')
    pyplot.plot(sigf[:,3], 'y-')
    return

class pca():
    def __init__(self,data):
        """
        data[i,:] is ith data

        creates 'var' (variances), 'base' (base[i,:] are the base,
        sorted from largest constributor to the least important. .coef
        """
        self.data     = data
        self.data_std = numpy.std(self.data, axis=0)
        self.data_avg = numpy.mean(self.data,axis=0)
     
        # variance co variance matrix
        self.varcovar  = (self.data-self.data_avg[numpy.newaxis,:])                      
        self.varcovar = numpy.dot(numpy.transpose(self.varcovar), self.varcovar)

        # diagonalization
        e,l,v= numpy.linalg.svd(self.varcovar)
       
        # projection
        coef = numpy.dot(data, e)
        self.var  = l/numpy.sum(l)
        self.base = e
        self.coef = coef
        return
    def comp(self, i=0):
        """
        returns component projected along one vector
        """
        return self.coef[:,i][:,numpy.newaxis]*\
               self.base[:,i][numpy.newaxis,:]

