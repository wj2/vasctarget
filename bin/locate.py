
import numpy as np

from numpy.fft import fftn, ifftn, fft2, ifft2
from scipy.signal import fftconvolve
from proc import normalize, pad_to_size_of, pad_by

def nccfft(s, p, room, fact=1):
    """Used for all Patterns that do not fall other categories.

    Cross correlates normalized Source and Pattern images while
    taking advantage of FFTs for the convolution step. 
    
    ----------
    s, p : Image
       Pattern and Source images for comparison.
       
    fact : int, optional
       Factor by which both Source and Pattern are
       scaled down.

    ----------
    out1 : ndarray[float]
       Confidence matrix for matches.

    out2 : float
       Threshold for deciding if a match has been found.

    out3 : float
       Mean of the confidence matrix.
    """
    
    # subtract mean from Pattern
    pmm = p - p.mean()
    pstd = p.std()
    n = p.size
    
    # make matrix of ones the same size as pattern
    u = np.ones(p.shape)

    # pad matrices (necessary for convolution)
    s = pad_by(s, room)

    upad = pad_to_size_of(u, s)
    pmmpad = pad_to_size_of(pmm, s)
    
    # compute neccessary ffts
    fftppad = fftn(pmmpad)
    ffts = fftn(s)
    fftss = fftn(s**2)
    fftu = fftn(upad)

    # compute conjugates
    cfppad = np.conj(fftppad)
    cfu = np.conj(fftu)

    # do multiplications and ifft's
    top = ifftn(cfppad * ffts)
    bot1 = n * ifftn(cfu * fftss)
    bot2 = ifftn(cfu * ffts) ** 2

    # finish it off!
    bottom = pstd * np.sqrt(bot1 - bot2)
    full = top / bottom

    return np.where(full.real.max() == full.real)

def nccfft2(s, p):
    
    sn = normalize(s)
    pn = normalize(p)

    prob = fftconvolve(sn, pn)

    return np.where(prob == prob.max())

