
import sys

import numpy as np

from proc import normalize, smoothg, clean_rotate
from scipy.misc import imrotate
from scipy.signal import argrelmax
from scipy.ndimage import gaussian_filter


def profile(section, probe, thresh=-0.5):
    dam = 0
    p = probe * section
    p = p.T[p.T.nonzero()]
    # if avg luminance is high (eg, on a vessel)
    pmean = p.mean()
    if pmean > 0.5: 
        dam += 1
        p[p > 1.5] = 1
    
    p = smoothg(p, 9) # smooth 
    p[p < thresh] = thresh # elim small peaks

    # find number of peaks (ie, local maxima)
    dam += argrelmax(p, order=2)[0].size
    
    return dam, p
        
def line_profiles(views, stack, thetas, down, gauss, probesize):
    # represent probe as single line
    l = probesize[0]
    probe = np.zeros((l, l))
    probe[l / 2] = 1
    
    # create container for damage map
    damage = np.zeros((len(thetas), 
                       (stack.shape[1] - l) / down + 1,
                       (stack.shape[2] - l) / down + 1))
    # container for individual line profiles
    lineprofiles = np.empty((stack.shape[0],
                             len(thetas), 
                             (stack.shape[1] - l) / down + 1,
                             (stack.shape[2] - l) / down + 1), 
                            dtype=object)

    # begin loop procedure
    for lnum, layer in enumerate(stack):
        sys.stdout.write('layer '+str(lnum+1)+' of '+str(stack.shape[0])+'\n')
        layer = normalize(layer)
        for i,t in enumerate(thetas):
            sys.stdout.write('at '+str(t)+' degrees\n')
            rotated = clean_rotate(probe, t, interp='nearest')
            ys = 0
            for y in xrange(0, layer.shape[0] - l + 1, down):
                xs = 0
                for x in xrange(0, layer.shape[1] - l + 1, down):
                    # djo line profile
                    dam, line = profile(layer[y:y+l, x:x+l],
                                        rotated)
                    damage[i, ys, xs] += dam
                    lineprofiles[lnum, i, ys, xs] = line
                    xs += 1
                ys += 1

    views['line'] = damage
    views['lineprof'] = lineprofiles
    if gauss:
        smoothdamage = gaussian_filter(damage, 1)
        views['linegauss'] = np.ceil(smoothdamage)
    return views

def luminance(views, stack, thetas, probesize):
    # construct probe representation
    probe = make_probe(probesize)
    # container for lumdamage
    damage = np.empty((stack.shape[0], len(thetas), stack[1], stack[2]))
    
    for lnum, layer in enumerate(stack):
        for i, t in enumerate(thetas):
            rotated = clean_rotate(probe, t, interp='nearest')
            damage[lnum, i] = ndimage.convolve(layer, rotated)

    views['luminance'] = damage.sum(axis=0)
    return views
    
