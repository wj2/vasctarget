"""
Library for functions related to digital signal processing
"""

import numpy as np
import utils
import data
import scipy.signal.filter_design as filter_design
from scipy.signal import lfilter, cheby1, firwin
import sys,os

try:
    import scipy.fftpack as F
except:
    import numpy.fft as F
# if NP.has_pystream:
#     from cufftpack import fft,ifft
# else:
try:
    from scipy.fftpack import fft,ifft
except:
    from numpy.fft import fft,ifft

try:
    import rpy
    from rpy import r as R
except:
    print "could not load R-project"

__all__ = ['denoise_time_series','analytic_signal','analytic_signal_fft','wmorlet','wfmorlet','wgauss','smoothg','smooth_timeseries',
           'islocmin','locmin','arglocmin','islocmax','locmax','arglocmax',
           'mtpsd','psd_sig_noise','psd_to_snr','correct_psd', 'iirfilter','downsample']

module_dir = os.path.split(__file__)[0]

def normalize(arr): 
    return (arr - arr.mean()) / arr.std()

def wmorlet(f0,s_f,samplingrate=10**3,ns=3,normed='area'):
    """
    returns a complex morlet wavelet in the time domain

    :Parameters:
        f0 : center frequency
        s_f : standard deviation
        samplingrate : samplingrate
        ns : window length in number of stanard deviations
    """
    s_t = 1./(2.*np.pi*s_f)
    w_sz = float(int(ns*s_t*samplingrate)) # half time window size
    t = np.arange(-w_sz,w_sz+1,dtype=float)/samplingrate
    if normed == 'area':
        w = np.exp(-t**2/(2.*s_t**2))*np.exp(
            2j*np.pi*f0*t)/np.sqrt(np.sqrt(np.pi)*s_t*samplingrate)
    elif normed == 'gain':
        w = np.exp(-t**2/(2.*s_t**2))*np.exp(
            2j*np.pi*f0*t)*2*s_f*np.sqrt(2*np.pi)/samplingrate
    else:
        assert 0, 'unknown norm %s'%normed
    return w

def wfmorlet(f0,s_f,samplingrate=10**3,ns=3,nt=None):
    """
    returns a complex morlet wavelet in the frequency domain

    :Parameters:
        f0 : center frequency
        s_f : standard deviation
        samplingrate : samplingrate
        ns : window length in number of stanard deviations
        nt : window length in number of sample points
    """
    if nt==None:
        s_t = 1./(2.*np.pi*s_f)
        nt = 2*int(ns*s_t*samplingrate)+1
    f = F.fftfreq(nt,1./samplingrate)
    wf = 2*np.exp(-(f-f0)**2/(2*s_f**2))*np.sqrt(samplingrate/(np.sqrt(np.pi)*s_f))
    wf[f<0] = 0
    wf[f==0] /= 2
    return wf

def wflogmorlet(f0,s_f,samplingrate=10**3,ns=3,nt=None):
    """
    returns a complex log morlet wavelet in the frequency domain

    :Parameters:
        f0 : center frequency
        s_f : standard deviation
        samplingrate : samplingrate
        ns : window length in number of stanard deviations
        nt : window length in number of sample points
    """
    if nt==None:
        s_t = 1./(2.*np.pi*s_f)
        nt = 2*int(ns*s_t*samplingrate)+1
    f = F.fftfreq(nt,1./samplingrate)

    sfl = np.log(1+1.*s_f/f0)
    wf = 2*np.exp(-(np.log(f)-np.log(f0))**2/(2*sfl**2))*np.sqrt(samplingrate/(np.sqrt(np.pi)*s_f))
    wf[f<0] = 0
    wf[f==0] /= 2
    return wf

def wlogmorlet(f0,s_f,samplingrate=10**3,ns=3,normed='area'):
    """
    returns a complex log morlet wavelet in the time domain

    :Parameters:
        f0 : center frequency
        s_f : standard deviation
        samplingrate : samplingrate
        ns : window length in number of stanard deviations
    """
    s_t = 1./(2.*np.pi*s_f)
    w_sz = int(ns*s_t*samplingrate) # half time window size
    wf = wflogmorlet(f0,s_f,samplingrate=samplingrate,nt=2*w_sz+1)
    w = F.fftshift(F.ifft(wf))
    if normed == 'area':
        w /= w.real.sum()
    elif normed == 'max':
        w /= w.real.max()
    elif normed == 'energy':
        w /= np.sqrt((w**2).sum())
    else:
        assert 0, 'unknown norm %s'%normed
    return w

def wgauss(nt=None,ns=3,sdt=None,normed='area'):
    """
    returns a gaussian window

    :Parameters:
        nt : window length in number of sample points
        ns : window length in number of stanard deviations
        sdt: standard deviation in number of sample points
    """
    if nt==None:
        nt = 2*int(ns*sdt)+1
    th = (nt-1)/2; # half window size
    s  = float(nt)/ns/2;  # standard deviation
    t = np.array(np.arange(-th,nt-th),'d')
    w = np.exp(-t**2/(2*s**2))
    if normed == 'area':
        w /= w.sum()
    elif normed == 'max':
        w /= w.max()
    elif normed == 'energy':
        w /= np.sqrt((w**2).sum())
    else:
        assert 0, 'unknown norm %s'%normed
    return w

def analytic_signal(f,sd,signal,win_type=0,neg=0,pos=0,ns=5,normed='area',logfreq=False):
    sig = signal.astype(float)
    if logfreq:
        w = wlogmorlet(f,sd,samplingrate=signal.i.samplingrate,ns=ns,normed=normed)
    else:
        w = wmorlet(f,sd,samplingrate=signal.i.samplingrate,ns=ns,normed=normed)
    nd = (w.shape[0]-1)/2
    if neg: w[:nd+1] = np.zeros(nd+1)     # only use negative times
    if pos: w[nd:2*nd+1] = np.zeros(nd+1) # only use positive times
    if win_type < 0:                   # don't use the first n time points
                                       # where n = -win_type
        w[nd:nd+1-win_type] = np.zeros(-win_type+1)
    elif (win_type == 'z')|(win_type == 'Z'):              # don't use time = t0
        w[nd] = 0
    a_signal = data.TimeSeries(np.zeros(signal.shape),dtype='D',
                            samplingrate=signal.i.samplingrate,i=signal.i)
    if signal.has_trials():
        for i in range(signal.has_trials()):
            a_signal[i,:] += (np.convolve(sig[i,:], np.real(w), mode='same') + 1j*
                              np.convolve(sig[i,:], np.imag(w), mode='same'))
    else:
        a_signal[:] += (np.convolve(sig, np.real(w), mode='same') + 1j*
                        np.convolve(sig, np.imag(w), mode='same'))
    return a_signal

def analytic_signal_fft(f,sd,signalf):
    t0 = utils.now()
    sigf = signalf.astype(float)
    wf = wfmorlet(f,sd,samplingrate=signalf.i.samplingrate,nt=len(sigf))
    asigf = sigf*wf
    asig = ifft(asigf)
    dt = utils.now()-t0
    print 'time needed: %s s' %dt
    return data.TimeSeries(asig,samplingrate=signalf.i.samplingrate,info=signalf.i)


def denoise_time_series(obj,**args):
    # default frequencies:
    noise = data.Info(freqs_line=[60,120,180,240,300,360],sf_line=.075,
                   freqs_refresh=[],sf_refresh=.01,
                   freqs_frame=[],sf_frame=.025)
    noise.update(obj.i)
    noise.update(args)
    
    ndata = len(obj)
    samplingrate = obj.i.samplingrate
    df = float(ndata)/samplingrate
    fft_data = fft(obj.astype('d'))
    filt_notch = np.ones(ndata,'d')
    
    # remove dc
    filt_notch[0] = 0
    # remove line noise +/- 0.075 Hz (incl. harmonics)
    sf = noise.sf_line
    for freq in noise.freqs_line:
        filt_notch[int((freq-sf)*df):int((freq+sf)*df)+1] = 0
        filt_notch[int((-freq-sf)*df):int((-freq+sf)*df)+1] = 0
    # remove refresh artifact  +/- 0.01
    sf = noise.sf_refresh
    for freq in noise.freqs_refresh:
        filt_notch[int((freq-sf)*df):int((freq+sf)*df)+1] = 0
        filt_notch[int((-freq-sf)*df):int((-freq+sf)*df)+1] = 0
    # remove frame artifact +/- 0.02
    sf = noise.sf_frame
    for freq in noise.freqs_frame:
        filt_notch[int((freq-sf)*df):int((freq+sf)*df)+1] = 0
        filt_notch[int((-freq-sf)*df):int((-freq+sf)*df)+1] = 0
    return data.TimeSeries(np.real(ifft(fft_data*filt_notch)),samplingrate,i=noise)


def smoothc(mI, Nr, Nc):
    """
    SMOOTHC.M: Smooths matrix data, cosine taper. MO=SMOOTHC(MI,Nr,Nc)
    smooths the data in MI using a cosine taper over 2*N+1 successive
    points, Nr, Nc points on each side of the current point.

           Inputs: mI - original matrix
                   Nr - number of points used to smooth rows
                   Nc - number of points to smooth columns
           Outputs:mO - smoothed version of original matrix

    """
    pass


def smoothg(sig,win_sz=None,ns=3,std=None,window=wgauss,**kargs):
    """
    smoothes signal using a gaussian window

    either win_sz (window size in bins) or sigma (std in s) has to be given
    """
    if std != None:
        # compute window size; sig has to be of class TimeSeries
        win_sz = int(np.floor(2*ns*std*sig.i.samplingrate))

    if isinstance(sig,data.TimeSeries):
        ret = data.TimeSeries(np.empty(sig.shape,dtype=sig.dtype),samplingrate=sig.i.samplingrate,i=sig.i.copy())
    else:
        ret = np.empty(sig.shape,dtype=sig.dtype)
    if sig.ndim == 1:
        if sig.dtype == np.dtype(complex):
            ret[:] = np.convolve(np.real(sig),window(win_sz,ns=ns,**kargs),mode='same'
                                )+1j*np.convolve(np.imag(sig),window(win_sz,ns=ns,**kargs),mode='same')
        elif sig.dtype == np.dtype(float):
            ret[:] = np.convolve(sig,window(win_sz,ns=ns,**kargs),mode='same')
        else:
            assert 0, 'bad typecode'
    elif sig.ndim == 2:
        for i in range(sig.shape[0]):
            ret[i,:] = smoothg(sig[i,:],win_sz,ns=ns,window=window,**kargs)
    else:
        assert 0, 'only 1D and 2D data supported'
    return ret


def smooth_timeseries(sig,width=1,**kargs):
        """
        smoothes timeseries using a gaussian window with standard deviation width (in s).

        optional parameters:
        ns : number of standard deviations used for window (on each side of the window center)
        """
        args = dict(ns=3,window=wgauss)
        args.update(kargs)
        args['win_sz'] = 2*int(args['ns']*sig.i.samplingrate*width)+1 # time points in gauss
        return smoothg(sig,**args)


def islocmin(sig):
    """
    Returns a boolean array of local minima.
    True if sig[n] < sig[n-1] <= sig[n+1].
    Array boundaries (undefined) are False.

    The algorithm used here is ~35 times faster than
    using an element-wise comparision loop in Python.
    """
    neg_slope = np.diff(sig) < 0
    ret = np.zeros(sig.shape,'bool')
    ret[1:-1] = neg_slope[:-1] # t-1 > t, neg slope at t-1
    ret[:-1] &= ~neg_slope # t <= t+1, pos slope at t
    return ret


def locmin(sig):
    """
    Returns array of locally minimal values.
    Calls islocmin.
    """
    return sig[islocmin(sig)]


def arglocmin(sig):
    """
    Returns indices of locally minimal values.
    Calls islocmin.
    """
    return np.where(islocmin(sig))[0]


def islocmax(sig):
    """
    Returns a boolean array of local maxima.
    True if sig[n] > sig[n-1] >= sig[n+1].
    Array boundaries (undefined) are False.

    The algorithm used here is ~35 times faster than
    using an element-wise comparision loop in Python.
    """
    pos_slope = np.diff(sig) > 0
    ret = np.zeros(sig.shape,'bool')
    ret[1:-1] = pos_slope[:-1] # t-1 < t, pos slope at t-1
    ret[:-1] &= ~pos_slope # t >= t+1, pos slope at t
    return ret


def locmax(sig):
    """
    Returns array of locally maximal values.
    Calls islocmax.
    """
    return sig[islocmax(sig)]


def arglocmax(sig):
    """
    Returns indices of locally maximal values.
    Calls islocmax.
    """
    return np.where(islocmax(sig))[0]

def multitapers(nt,tw):
    """
    Returns the slepian window for multitaper spectral power estimation

    (at this point, the function loads pre-computed tapers)
    """
    k = np.floor(2*tw-1)
    try:
        rpy.set_default_mode(rpy.BASIC_CONVERSION)
        R.library("waveslim")
        tapers = R.dpss_taper(nt,k,tw)
    except:
        # load tapers produced by matlab
        from scipy.io import loadmat
        fn = os.path.join(module_dir,'examples','sample_data','dpss-'+str(nt)+'-'+str(tw))
        mydict = loadmat(fn)
        tapers = np.array(mydict['tapers'])
    return tapers

def mtpsd(sig,tw=3,nwin=1024,remove_mean=1):
    """
    Spectral power estimation using multitaper method
    """
    # compute frequencies
    samplingrate = sig.i.samplingrate
    nfft = int(2**(np.floor(np.log(nwin)/np.log(2))+1))
    f = F.fftfreq(nfft,1./samplingrate)
    ind = (f>=0)&(f<=samplingrate/2.)
    f = f[ind]

    # prepare sig
    mydata = sig.view(np.ndarray).copy()
    if sig.ndim == 1:
        mydata.shape = (1,len(mydata)) # make data 2-d
    (repeats,ndata) = mydata.shape
    nt = ndata/nwin
    if remove_mean:
        mydata -= mydata.mean(axis=1)[:,np.newaxis]

    # padding to integer multiple of window size
    if np.mod(ndata,nwin):
        pad = nwin-np.mod(ndata,nwin)
        mydata = np.hstack((mydata,np.zeros((repeats,pad),'d')))
        nt+=1

    # compute multitapers
    tapers = multitapers(nwin,tw)
    k = tapers.shape[1]
    tapersn = np.transpose(tapers*np.sqrt(samplingrate))
    tapersn = np.transpose(np.array([tapersn]*repeats),axes=[1,0,2])

    # compute psd
    psd = np.zeros((repeats,len(f)),'d')
    for i in range(nt):
        dat = mydata[:,i*nwin:(i+1)*nwin]
        dat_proj = tapersn*np.array([dat]*k,'d')
        ft_tmp = F.fft(dat_proj,n=nfft)
        ft = ft_tmp[...,ind]
        psd += np.mean(abs(ft)**2,axis=0)
    psd /= nt
    if sig.ndim == 1: psd = np.squeeze(psd)
    return f,psd

def psd_sig_noise(sig,signal=None,nwin=1024,correct_repeats=False,**args):
    """
    Computes spectral power estimate for signal- and noise for repeated trials.

    If signal is not given, the trial mean is assumed to be the signal.
    If correct_repeats is given, the signal power is corrected down for finite trial number.
    """

    # compute mean sig power
    (f,psd_data) = mtpsd(sig,nwin=nwin,**args)
    (f,psdmean) = mtpsd(sig.mean(axis=0),nwin=nwin,**args)

    psd_tot = psd_data.mean(axis=0)
    if correct_repeats:
        psd_sig,psd_noise = correct_psd(psdmean,psd_tot,repeats=correct_repeats)
    else:
        psd_sig,psd_noise = (psdmean,psd_tot-psdmean)

    if signal!=None: # compute noise as deviation from signal
        if signal.ndim == 1:
            signal = signal[np.newaxis,:]   # make signal 2-D
            noise = signal-sig
            f,psd_noise = mtpsd(noise,nwin=nwin,**args)
            psd_noise = psd_noise.mean(axis=0)

    return f,psd_sig,psd_noise

def correct_psd(psdmean,psd_tot,repeats=None,noisemin=np.finfo(float).eps):
    """
    Corrects the spectral estimate for finite number of trials

    psdmean = psd(mean) is the spectral power of the average trial,
              considered to be the signal power (for repeats->infinity)
    psd_tot = mean(psd) is the average power of a trial

    noise_min is a lower bound for the noise power in each freqency band
         to ensure that the signal-to-noise ratio is finite

    psd_noise is the estimate for the noise power, computed as
      psd_noise = cfactor*(psd_tot-psdmean)

    cfactor is a correction factor, correcting for the fact that for finite number
    of repeats, psdmean contains a contribution of 1./repeats * psd_noise.
    """
    psd_noise = (psd_tot-psdmean)
    if repeats: psd_noise /= (1-1./repeats) # correct for finite repeats
    psd_noise = np.where(psd_noise<psd_tot,psd_noise,psd_tot)   # psd_sig < psd_tot
    psd_noise = np.where(psd_noise>noisemin,psd_noise,noisemin) # psd_noise > eps
    psd_sig = psd_tot-psd_noise
    return psd_sig,psd_noise


def psd_to_snr(f,psd_sig,psd_noise,f_cutoff=np.inf):
    """
    Returns signal-to-noise ration and information of Gaussian chanel
    """
    df = f[1]
    imax = np.argmax(f[f<f_cutoff])+1
    snr = psd_sig[1:imax]/psd_noise[1:imax]
    inf = np.log(snr+1).sum()/np.log(2)*df
    return f[1:imax],snr,inf

def iirfilter(signal, N_ord=3,Wn=[.5, 150.0],samplingrate=None,):
    """
    Filter signal with iirfilter specified by its order and corner frequencies
    
    :Parameters:
        signal : the signal to be filtered
        N : order of the filter
        Wn : corner frequency (or pair of corner frequencies)
        samplingrate : the sampling rate of the signal

    Returns an ndarray.
    """
    if samplingrate is None:
        samplingrate = signal.i.samplingrate
    nyq = samplingrate / 2.0
    Wn = np.array(Wn)
    b,a = filter_design.iirfilter(N_ord,Wn/nyq,signal)
    return lfilter(b,a,signal);
    
def downsample(signal, q, n=None, ftype='iir', axis=-1):
    """
    Downsample signal x by an integer factor q, using an order n filter

    By default, an order 8 Chebyshev type I filter is used or a 30 point FIR 
    filter with hamming window if ftype is 'fir'.

    (port to Python of the GNU Octave function decimate.)

    :Parameters:
        signal -- the signal to be downsampled (N-dimensional array)
        q -- the downsampling factor
        n -- order of the filter (1 less than the length of the filter for a
             'fir' filter)
        ftype -- type of the filter; can be 'iir' or 'fir'
        axis -- the axis along which the filter should be applied
    """
    if type(q) != int: raise Error, "q should be an integer"
    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 8
    if ftype == 'fir':
        b = firwin(n+1, 1./q, window='hamming')
        y = lfilter(b, 1., signal, axis=axis)
    else:
        (b, a) = cheby1(n, 0.05, 0.8/q)
        y = lfilter(b, a, signal, axis=axis)
    return y.swapaxes(0,axis)[::q].swapaxes(0,axis)
