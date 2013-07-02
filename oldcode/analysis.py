
import numpy as np

def normalize(arr): 
    return (arr - arr.mean()) / arr.std()

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

def smoothg(sig,win_sz=None,ns=3,std=None,window=wgauss,**kargs):
    """
    smoothes signal using a gaussian window

    either win_sz (window size in bins) or sigma (std in s) has to be given
    """
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
