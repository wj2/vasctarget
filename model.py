"""
This module contains classes and methods related to data modeling.
"""

import numpy as np
try:
    import rpy
    from rpy import r as R
    R.library("CircStats")
except:
    print "could not load R-project"

__all__ = ['phasedist','mises','mises_params','cvar','kappa_to_cvar','fit_gamma','fit_exp_nonlinearity','likelihood','rate_information']

def phasedist(phases):
    """
    Returns parameters of Von Mises distribution fitted to phase data
    """
    direction = np.mean(np.exp(1j*phases))
    n = len(phases)
    (kappa,mu,p) = mises_params(direction,n)
    return (kappa,mu,p,n)

def mises(phi,kappa,mu):
    """
    Returns the Von Mises distribution with mean mu and concentration kappa
    """
    from scipy.special import i0
    return (1./(2.*np.pi*i0(kappa)))*np.exp(kappa*np.cos(phi-mu))

def mises_params(direction,n=1):
    """
    Computes parameters of Von Mises distribution from direction vector
    """
    from scipy.optimize import fmin
    from scipy.special import i0,i1
    def bess(x,r):
        return (i1(x)/i0(x)-r)**2
    try: # using R (faster by a factor of 10)
        kappa = R.A1inv(np.abs(direction))
    except:
        kappa = float(fmin(bess,np.array([1.]),(np.abs(direction),),disp=0));
    mu = np.angle(direction)
    z = float(n)*np.abs(direction)**2
    p = np.exp(-z)*(1.+(2.*z-z**2)/(4.*n)-
                 (24.*z-132.*z**2+76.*z**3-9.*z**4)/(288.*n**2))
    return kappa,mu,p

def fit_gamma(data):
    """
    fit gamma distribution to data and returns mean and shape parameter
    """
    xmean = np.array(data,'d').mean()
    k = xmean**2/np.array(data,'d').var()
    print "gamma k = %1.1f, mean = %f" %(k,xmean)
    return k,xmean

def gamma(t,k,xmean):
    from scipy import special
    kappa = k/xmean
    return (kappa**k)*(t**(k-1))*np.exp(-kappa*t)/special.gamma(k)

def cvar(phases):
    """
    returns the circular variance of phase distribution
    """
    direction = np.mean(np.exp(1j*phases))
    return 1-np.abs(direction)**2

def kappa_to_cvar(kappa):
    """
    computes circular variance from phase concentration parameter kappa
    """
    from scipy.special import i0,i1
    return 1-(i1(kappa)/i0(kappa))**2

def prinComp2(x,k):
    """
    PRINCOMP Principal Component Analysis (centered and scaled data).
    [PC, SCORE, LATENT, TSQUARE] = PRINCOMP(X) takes a data matrix X and
    returns the principal components in PC, the so-called Z-scores in SCORES,
    the eigenvalues of the covariance matrix of X in LATENT, and Hotelling's
    T-squared statistic for each data point in TSQUARE.
    """
    pass


def normalDist(x,params):
    """
    returns a normal distribution with params "u" and "sig2"
    given by: p(r) = sqrt(1/2/pi/sig^2)*exp(-((r-)^2)/(2*sig^2))
    """
    pass


def nCr(n,r):
    """
    result = nCr(n,r)

    Returns n choose R.  n can be an array or matrix.
    """
    pass


def nPr(n,r):
    """
    result = nPr(n,r)

    Returns n pick R. n can be an array or matrix.
    """
    pass


def fitMLnormal( x,hAx ):
    """
    fit_ML_normal - maximum likelihood fit of the normal distribution of i.i.d.
    samples. Given the samples of a normal distribution, the PDF parameters are
    returned. Fits data to the probability of the form:
        p(r) = sqrt(1/2/pi/sig^2)*exp(-((r-u)^2)/(2*sig^2))
    ...with parameters: u,sig^2

    input:    x   - vector, samples with normal distribution to be parameterized
            hAx - handle of an axis, on which the fitted distribution is plotted
                  if h is given empty, a figure is created.

    output:   result  - structure with the fields
                     sig^2,u          - fitted parameters
                     CRB_sig2,CRB_u   - Cramer-Rao Bound for the estimator value
                     RMS              - RMS error of the estimation
                     type             - 'ML'
    """
    pass


def fit_exp_nonlinearity(spktrain,**signals):
    """
    Fit exponential nonlinearity given spktrain and feed forward signal
    """

    rpy.set_default_mode(rpy.NO_CONVERSION)
    datadict = dict(sp=spktrain)
    datadict.update(signals)
    data = R.as_data_frame(datadict)
    vars = signals.keys()
    obj = R.glm('sp ~ '+' + '.join(vars), data=data, family='poisson')
    rpy.set_default_mode(rpy.BASIC_CONVERSION)

    # print obj.as_py()['coefficients']
    print R.coef(obj)
    coefmatrix = R.summary(obj)['coefficients']
    for i,var in enumerate(['(Intercept)',]+vars):
        print '#%d %s = %f +/- %f \t (z = %f, p = %f)'%tuple([i,var,]+coefmatrix[i].tolist())

    return R.coef(obj)

def likelihood(rate,events):
    """
    likelihood of poisson spike events with varying density function
    """
    ntrials = events.has_trials()
    Lf = np.zeros(ntrials)
    dt = rate.time().max()-rate.time().min()
    samplingrate = rate.i.samplingrate
    for i in range(ntrials):
        if rate.has_trials():
            Lf[i] += np.log(rate[i][events.select_dim(i)]).sum()
            Lf[i] -= rate[i].sum()/samplingrate
        else:
            Lf[i] += np.log(rate[events.select_dim(i)]).sum()
            Lf[i] -= rate.sum()/samplingrate
        Lf[i] /= dt
    Lfmean = Lf.mean()
    Lfstderr = Lf.std()/np.sqrt(ntrials)
    print 'Likelihood = %3.2f +/- %2.2f bit/s' %(Lfmean,Lfstderr)
    return (Lfmean,Lfstderr)

def rate_information(rate,meanrate=None):
    """
    computes the shannon information carried by an event based on the 1d event rate in bit/event

    see
      Brenner N, Strong SP, Koberle R, Bialek W, de Ruyter van Steveninck RR.
      Synergy in a neural code. Neural Comput. 2000 Jul;12(7):1531-52.
    """
    if rate.ndim > 1: raise NotImplementedError
    if meanrate is None: meanrate = rate.mean()
    p = rate.astype('d')/meanrate
    p[p>0]=p[p>0]*np.log2(p[p>0])
    return float(p.mean())


