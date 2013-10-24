
import sys, csv, time, re
import gui

import numpy as np

from proc import *
from scipy.misc import imrotate
from scipy.signal import argrelmax
from scipy.ndimage import gaussian_filter
from scipy.stats.mstats import kurtosis
from scipy.stats import skew

def attempt_meaningful_name(pre, post):
    psize_re = '-[0-9]{1,4}um-'
    date_re = '-[0-9]{8}-'
    ppost_des_re = '-((pre)|(post))-[0-9]{3}'
    try: 
        psize = re.search(psize_re, pre).group(0)
        date = re.search(date_re, pre).group(0)
        pre_des = re.search(ppost_des_re, pre).group(0)
        post_des = re.search(ppost_des_re, post).group(0)

        name = 'ZSeries'+date+psize.strip('-')+pre_des+post_des+'.csv'

    except:

        t = time.strftime('%X%x%Z').replace(' ','').replace(':','.')
        t = t.replace('/','-')
        
        name = 'damage'+t+'.csv'

    return name
    

def print_damage_stats(pre, post, masks):
    filename = attempt_meaningful_name(pre[0], post[0])
    with open(filename, 'wb') as cf:
        writ = csv.writer(cf)
        
        writ.writerow([pre[0]]+['' for x in xrange(len(masks) - 1)]
                      +[post[0]]+['' for x in xrange(len(masks) - 1)])
        writ.writerow([m.shape for m in masks] + [m.shape for m in masks])
        for i,m in enumerate(pre[1]):
            writ.writerow(list(pre[1][i])+list(post[1][i]))
            
    return None

def print_misc_stats(path, alldata, key, write=True): 
    name = path.replace('/', '-').replace('.','l') + '_' + key + '_stats'
    want = [np.mean, np.median, np.min, np.max, skew, kurtosis]
    names, results = ['data'], [path]
    for x in want:
        names.append(x.__name__)
        results.append(x(alldata[key], axis=None))
        
    if write:
        with open(name+'.csv', 'wb') as cf:
            writ = csv.writer(cf)
            writ.writerow(names)
            writ.writerow(results)
        return None
    else:
        return names, results

def print_tddict(dic, title):
    t = time.strftime('%X%x%Z').replace(' ','').replace(':','.').replace('/','-')
    filename = 'stats_'+title+'_'+t+'.csv'
    
    with open(filename, 'wb') as cf:
        writ = csv.writer(cf)
        writ.writerow(dic['title'])
        for l in dic['data']:
            writ.writerow(l)
    

def profile(p, thresh=-0.7, upper=1.5, smooth=20, peaksize=5):
    dam = 0
    
    smooth_only = smoothg(p, smooth) # smooth 
    p = smooth_only
    # if whole luminance is high (eg, on a vessel)
    if p.min() >= upper: 
        dam += 1
    
    # count high peaks only once
    t = p > upper
    dam += np.diff(t).nonzero()[0][::2].size
    p[p > upper] = upper

    p[p < thresh] = thresh # elim small peaks

    # find number of peaks (ie, local maxima)
    dam += argrelmax(p, order=peaksize)[0].size

    return dam, p, smooth_only
        
def line_profiles(views, stack, thetas, down, gauss, probesize):
    # avoid big top level vessels, how to do it visually
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
                    # do line profile
                    section = layer[y:y+l, x:x+l]
                    p = rotated * section
                    p = p.T[p.T.nonzero()]

                    dam, line, onsm = profile(p)
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
    
def damage_profiles(stack, locs, rots, psize, args):
    thresh = -0.7; upper = 1.5
    # get items out of args namespace
    inter = args.interval
    n = args.n_profiles
    debug = args.debug
    if args.smooth:
        sm = 25
    else:
        sm = 20

    masks = create_masks(psize, inter, n)
    if debug or args.manual:
        import matplotlib.pyplot as plt
    count = np.zeros((stack.shape[0], n + 1))
    lines = np.empty((stack.shape[0], n + 1), dtype=object)
    for i, layer in enumerate(stack):
        # print rots[i], locs[i]
        rotlayer, yx = rotate_mask_horizontal(layer, rots[i], locs[i])
        rotlayer = normalize(rotlayer)
        # print rotlayer.shape, yx
        for j, mask in enumerate(masks):         
            # print mask.shape

            my_shape = mask.shape[0]/2.0
            mx_shape = mask.shape[1]/2.0
            myb_shape = mask.shape[0]/2.0 + args.buffer
            mxb_shape = mask.shape[1]/2.0 + args.buffer
            super_section = rotlayer[yx[0]-myb_shape:yx[0]+myb_shape,
                                     yx[1]-mxb_shape:yx[1]+mxb_shape]
            section = super_section[args.buffer:-args.buffer,
                                    args.buffer:-args.buffer]
            dat = mask * section
            prof = collapse_rect_mask(dat)
        
            count[i,j], lines[i,j], onsm = profile(prof, thresh=thresh, 
                                                   upper=upper, smooth=sm)
            if args.manual:
                results = gui.make_thresh_setter(prof, thresh, upper, mask, 
                                                 super_section, args.buffer, sm)
                count[i,j], lines[i,j], args.manual = results
            if debug:
                printfig = plt.figure()
                layax = printfig.add_subplot(211)
                layax.imshow(rotlayer, 'gray', interpolation='none')
                pim = layax.imshow(mask, extent=[0,0,0,0], interpolation='none')
                pim.set_extent([yx[1]-(mask.shape[1]/2.0),
                                yx[1]+(mask.shape[1]/2.0),
                                yx[0]-(mask.shape[0]/2.0),
                                yx[0]+(mask.shape[0]/2.0)])
                layax.set_ylim(rotlayer.shape[1], 0)
                layax.set_xlim(0, rotlayer.shape[0])
                layax.set_title(count[i,j])

                lineax = printfig.add_subplot(212)
                lineax.plot(lines[i, j])

                printfig.savefig(str(i)+str(mask.shape)+'.pdf', format='pdf')
                plt.close(printfig)


    return count, lines, masks
