
import numpy as np

from numpy.fft import fftn, ifftn, fft2, ifft2
from scipy.signal import fftconvolve
from proc import normalize

def nccfft(s, p, fact=1):
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
    s = pad_one_half(s)

    upad = pad_to_size_of(u, s)
    pmmpad = pad_to_size_of(pmm, s)
    
    print s, pmmpad

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

    print full.real.max(), full.real.min(), full.real.shape

    return np.where(full.real.max() == full.real)

def nccfft2(s, p):
    
    sn = normalize(s)
    pn = normalize(p)

    prob = fftconvolve(sn, pn)

    print prob.min(), prob.max(), prob.shape
    print np.where(prob == prob.max())

    return np.where(prob == prob.max())

def pad_zeroes(m, axis, thick, side='post'):
    """Pads given matrix with zeros. 

    Pads m with thick zeros on axis. 

    ----------
    m : ndarray
       Matrix to be padded.

    axis : int
       Axis of m to be padded.

    thick : int
       Number of zeros to pad with.

    ----------
    out1 : ndarray
       Simply m, padded with given number of zeros on given axis.
    """

    if axis == 1:
        rows = m.shape[0]
        cols = thick
    elif axis == 0:
        rows = thick
        cols = m.shape[1]

    if side == 'pre':
        con = (np.zeros((rows, cols)), m)
    elif side == 'post':
        con = (m, np.zeros((rows, cols)))
    elif side == 'both':
        con = (np.zeros((rows, cols)), m, np.zeros((rows, cols)))
        
    return np.concatenate(con, axis)

# pad a to size of b, with zeroes
def pad_to_size_of(a, b):
    """Pads first matrix to the size of second matrix.

    Uses zeros to pad first matrix to size of second matrix.
    Will fail if first matrix is larger than second. 

    ----------
    a, b : ndarray

    ----------
    out : ndarray
       Original first array padded with zeros s.t it is
       now the size of b. 
    """

    a2 = pad_zeroes(a, 0, b.shape[0] - a.shape[0])
    return pad_zeroes(a2, 1, b.shape[1] - a.shape[1])

def pad_one_half(a):
    a2 = pad_zeroes(a, 0, a.shape[0] / 4, 'both')
    return pad_zeroes(a2, 1, a.shape[1] / 4, 'both')
                 
