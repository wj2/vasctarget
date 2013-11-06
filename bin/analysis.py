
import sys, csv, time, re, os
import gui

import numpy as np

from proc import *
from scipy.misc import imrotate
from scipy.signal import argrelmax, find_peaks_cwt
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
    
def save_damage(data_name, data, args, psize):
    name = data_name.split('/')
    if name[-1] == '':
        name = name[-2]
    else:
        name = name[-1]
    name = os.path.splitext(name)[0]+'-damage.npz'
    file_ = open(name, 'wb')
    np.savez(file_, damage=data, lw=np.array(psize),
             context=np.array((args.downsample, 
                               args.z_thickness,
                               args.channel)))

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
    
def decide_thresh(region):
    
    lower = 0; upper = 0

    return lower, upper

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

def profile_cwt(p, sizes=[10]):
    
    peaks = find_peaks_cwt(p, np.array(sizes))

    print 'peaks: '+str(peaks)
    return len(peaks), p, p
    
    
        
def line_profiles(views, stack, thetas, down, gauss, probesize):
    # avoid big top level vessels, how to do it visually
    # represent probe as single line
    l = probesize[0]
    probe = np.zeros((l, l))
    probe[l / 2] = 1
    
    # create container for damage map
    damage = np.zeros((stack.shape[0],
                       len(thetas), 
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
                    damage[lnum, i, ys, xs] = dam
                    lineprofiles[lnum, i, ys, xs] = line
                    xs += 1
                ys += 1

    views['line'] = damage.sum(axis=0)
    views['line_notreduced'] = damage
    views['lineprof'] = lineprofiles
    # save damage data 
    if gauss:
        smoothdamage = gaussian_filter(views['line'], 1)
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
            # these give negative overshoot
            pre_my_zeros = abs(min(yx[0] - my_shape, 0))
            post_my_zeros = max((yx[0] + my_shape) - layer.shape[0], 0)
            mx_shape = mask.shape[1]/2.0
            pre_mx_zeros = abs(min(yx[1] - mx_shape, 0))
            post_mx_zeros = max((yx[1] + mx_shape) - layer.shape[1], 0)
            myb_shape = mask.shape[0]/2.0 + args.buffer
            mxb_shape = mask.shape[1]/2.0 + args.buffer
            super_section = rotlayer[max(0, yx[0]-myb_shape):yx[0]+myb_shape,
                                     max(0, yx[1]-mxb_shape):yx[1]+mxb_shape]
            section = rotlayer[max(yx[0]-my_shape, 0):yx[0]+my_shape,
                               max(yx[1]-mx_shape, 0):yx[1]+mx_shape]

            print yx[0]-myb_shape, yx[0]+myb_shape
            print yx[1]-mxb_shape, yx[1]+mxb_shape
            print pre_my_zeros, post_my_zeros, pre_mx_zeros, post_mx_zeros
            print section.shape, super_section.shape
            dat = apply_mask(mask, section, [pre_my_zeros, post_my_zeros],
                             [pre_mx_zeros, post_mx_zeros])
            prof = collapse_rect_mask(dat)
            thresh = super_section.mean()
            upper = thresh + super_section.std()
        
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
