
import math

import numpy as np

from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate

def rect_mask(l, w):
    mask = np.zeros((w, l))
    mask[0] = 1; mask[-1] = 1
    mask[:, 0] = 1; mask[:, -1] = 1
    
    return mask

def create_masks(psize, inter, n):
    masks = [rect_mask(psize[0], 1)]
    for i in xrange(n):
        add = i * inter
        masks.append(rect_mask(psize[0] + add, psize[1] + add))

    return masks

def rotate_mask_horizontal(arr, rot, loc):
    dtr = math.pi / 180.0 # conv radians

    # put origin of loc at center
    cent = np.array(arr.shape) / 2
    y1, x1 = loc - cent
    theta1 = np.arctan(y1 / x1)
    theta2 = theta1 + (dtr * -rot)
    h = np.sqrt(y1**2 + x1**2)
    
    rotted = rotate(arr, rot)
    newcent = np.array(rotted.shape) / 2
    
    # print x1, '->', h * np.cos(theta2)
    # print y1, '->', h * np.sin(theta2)

    y2 = h * np.sin(theta2) + newcent[0]
    x2 = h * np.cos(theta2) + newcent[1]
    
    return rotted, np.array((y2, x2))
    
def collapse_rect_mask(mask):
    if min(mask.shape) < 2 or len(mask.shape) < 2:
        return mask.T[mask.T.nonzero()]

    else:
        return np.hstack((mask[0,1:-1], mask[:, -1], 
                          mask[-1, :-1][::-1], 
                          mask[:-1, 0][::-1]))

def make_probe(psize):
    l = psize[0]; w = psize[1]
    probe = np.zeros((l, l))
    probe[(l/2) - (w/2):(l/2) + (w/2), :] = 1
    return probe

def collapse_stack(stack, delta, factor=50):
    factor = factor / delta

    numSlices = int(np.ceil(stack.shape[0] / factor))

    if numSlices != stack.shape[0]:
        holster = np.empty((numSlices, stack.shape[1], stack.shape[2]))
        for i in xrange(numSlices):
            holster[i] = np.mean(stack[i*factor:i*factor + factor],
                                 axis=0)
        return holster
    else:
        return stack

def clean_rotate(probe, t, interp='nearest'):
    rotated = imrotate(probe, t, interp='nearest')
    # get rid of artifacts
    rotated[rotated < rotated.max()] = 0
    rotated[rotated == rotated.max()] = 1
    return rotated

def normalize(arr): 
    return (arr - arr.mean()) / arr.std()

def get_thetas(num):
    inter = 180 / num
    return xrange(0, 180, inter)


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
