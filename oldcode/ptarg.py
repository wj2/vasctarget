#!/usr/bin/python

import numpy as np, dsp.dsplib as dsp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
from scipy import ndimage
from scipy import misc
from scipy.signal import argrelmax, find_peaks_cwt
from random import choice

# replace this section with a few lines
import vasc.vastifarr as va

from config import path, probesize

a = va.tiffile(path)
tiflist = a.gettiflist(path)
tiflist2, xml = a.get_tiffs_and_xml(tiflist)
scales = a.get_scale(xml)
tiflist2 = tiflist2[400:450]
nfiles = a.listsize(tiflist2)
array3d = a.initarray(nfiles)
tif3darray = a.tif2array(path,tiflist2,array3d)
###########################################

regionSize = 50.0

meanvasdam = np.mean(tif3darray,axis=0)

cross = (40, 10)
def mouseclick(event):
    if event.inaxes == ax1:
        xdat = event.xdata * downsample
        ydat = event.ydata * downsample
        outfig = plt.figure()
        ax = outfig.add_subplot(111)
        axim = ax.imshow(subz[0], 'gray', interpolation='none')
        cross1 = plt.imshow(np.ones((cross[1], cross[0])), 
                           extent=[0,0,0,0])
        cross2 = plt.imshow(np.ones((cross[0], cross[1])),
                           extent=[0,0,0,0])
        extent1 = [xdat, xdat + cross[0],
                   ydat + cross[0]/2 - cross[1]/2, 
                   ydat + cross[0]/2 + cross[1]/2]
        extent2 = [xdat + cross[0]/2 - cross[1]/2,
                   xdat + cross[0]/2 + cross[1]/2,
                   ydat, ydat + cross[0]]
        cross1.set_extent(extent1)
        cross2.set_extent(extent2)
        ax.set_ylim(subz[0].shape[0], 0)
        ax.set_xlim(0, subz[0].shape[1])
        # ax2im.set_data(subz[0])
        # plt.draw()
        # extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        # fig.savefig('out.pdf', bbox_inches=extent, format='pdf')
        outfig.savefig('out.pdf', format='pdf')
        print 'pdf saved'
    elif event.inaxes == ax2:
        pass
# 1 count and compare to compute (3 or 4 times)
# 2 output, arrow key through best possible 
# 3 compare where suggested with where she went
def mousemove(event):
    # plot best orientation at a given location
    if event.inaxes == ax1 or event.inaxes == ax2:
        
        if event.inaxes == ax1:
            xdat = event.xdata - kxy/(2*downsample)
            ydat = event.ydata - kxy/(2*downsample)
    
            i = np.argmin(damage[:,event.ydata, 
                                   event.xdata])
            dam = np.min(damage[:,event.ydata,event.xdata])
            extent = [event.xdata * downsample, (event.xdata * downsample) + kxy, 
                      event.ydata * downsample, (event.ydata * downsample) + kxy]
            linep = lineprof[0, i, event.ydata, event.xdata]
        elif event.inaxes == ax2:
            xdat = event.xdata - kxy/2
            ydat = event.ydata - kxy/2

            i = np.argmin(damage[:,ydat / downsample, 
                                   xdat / downsample])
            dam = np.min(damage[:,ydat / downsample,
                                  xdat / downsample])
            extent = [xdat, xdat + kxy, 
                      ydat, ydat + kxy]
            linep = lineprof[0, i, ydat / downsample, xdat / downsample]

        img.set_data(probe[i])
        img.set_extent(extent)
        ax2.set_ylim(meanvasdam.shape[0],0)
        ax2.set_xlim(0,meanvasdam.shape[1])
        print('damage=%1d'%dam)
        
        ax3.clear()
        ax3.plot(linep)

        plt.draw()
    elif event.inaxes == ax4:
        i = np.argmin(damage[:,event.ydata,event.xdata])
        img.set_data(probe[i])
        extent = [event.xdata - kxy/2, 
                  event.xdata + kxy/2, 
                  event.ydata - kxy/2, 
                  event.ydata + kxy/2]
        img.set_extent(extent)
        ax2.set_ylim(512,0)
        ax2.set_xlim(0,512)
        print('damage=%1d'%np.min(lumdam[:,event.ydata,event.xdata]))

        plt.draw()



def micron2pix(probsize, reverse=False):       
    scales = a.get_scale(xml)
    if reverse:
        return (probsize[0] * scales[0], probsize[1] * scales[1])
    s = (1 / scales[0], 1 / scales[1])
    newprobsize = [int(probsize[0]*s[0]),int(probsize[1]*s[1])]
    return newprobsize

def build_kernel(probesizemicron, theta = [0]):
    # make mask of probe cross-section
    probesizepix = micron2pix(probesizemicron)
    l = max(probesizepix)
    w = 1# min(probesizepix)
    k = np.zeros((len(theta),l,l))    
    k[0,(l/2.-w/2.):(l/2.+w/2.),:] = 1 # horizontal probe
    for i, t in enumerate(theta):
        k[i] = ndimage.rotate(k[0],t,reshape=False)
    return k
    
def build_histogram(probesizemicron, theta = [0]):
    # make mask of probe cross-section
    probesizepix = micron2pix(probesizemicron)
    l = max(probesizepix)
    w = min(probesizepix)
    k = np.zeros((l,l))    
    k[(l/2.-w/2.):(l/2.+w/2.),:] = 1 # horizontal probe
    damage2d = np.zeros((len(theta),512,512))
    for i, t in enumerate(theta):
        k = ndimage.rotate(k,t,reshape=False)
        print('%s degrees'', all translations...'%t)
        damage2d[i] = ndimage.convolve(meanvasdam,k)
    return damage2d

def subz_from_3dtiff(tiffs, slicesize):
    slices = xrange(0, tiffs.shape[0], slicesize)
    holder = np.ones((len(slices), tiffs.shape[1], 
                      tiffs.shape[2]))
    i = 0
    while i < len(slices):
        print str(i*slicesize) +' to '+ str(i*slicesize+slicesize)
        holder[i] = np.mean(tiffs[i*slicesize:i*slicesize + slicesize], 
                            axis=0)
        i += 1
    
    return holder

def find_events(psth, thresh):
    """
    Finds events in a PSTH. 

    Given a 'psth', returns 'peaks' and 'epochs' corresponding to 
    maximal firing events and the boundaries of these events. Event
    boundaries are defined as significant (p < 0.05) minima between
    adjacent peaks as described in several Warland and Reinegal
    papers [refs needed]. Minima are defined relative to the
    geometric mean of adjacent maxima, that is: sqrt(Mi*M(i+1)) <= 3
    
    """
    # find PSTH peaks
    psth = np.array(psth) - thresh
    maxidx = argrelmax(psth, order=5)[0]
    print maxidx
    if maxidx.size < 2:
        return maxidx
    pkval  = psth[maxidx]
    # pktime = edges[maxidx]

    # find PSTH epochs
    mint=[0]
    minv=[0]
    peaksv=[]
    for j in xrange(len(maxidx)):
        if j == 0:
            thresh = np.sqrt(pkval[j]*pkval[j+1])
            mi = psth[maxidx[j]:maxidx[j+1]].argmin()+maxidx[j]
        elif j == len(maxidx) - 1:
            thresh = np.sqrt(pkval[j]*pkval[j-1])
            mi = psth[maxidx[j-1]:maxidx[j]].argmin()+maxidx[j-1]
        else:
            thresh = (pkval[j-1]*pkval[j]*pkval[j+1]) ** (1/3.0)
            mi = psth[maxidx[j-1]:maxidx[j+1]].argmin()+maxidx[j-1]
        m  = psth[mi]
        print m, thresh
        if m <= thresh / 3.0:
            minv.append(m)
            # mint.append(edges[mi])
            peaksv.append(pkval[j])
            # peakst.append(pktime[j])

    return peaksv

# downsample factor
downsample = 3
def count_vessels(probesizemicron, theta = [0]):
    # algorithm pseudocode:
    # normalize image slice
    # rotate probe
    # slice rotated image to generate line profiles 
    # low pass filter line profiles to remove noise
    # find peaks of line profile (local maxima)
    # discard peaks below empirical threshold  
    # count the number of peaks --> number of crossed vessels at  location
    # no need to adjust discard theshold with depth
    # repeat above for all substacks in the whole imaging volume
    # p.s. tested on 50um substacks
    thresh = -0.6
    # create sub-z tiff stack
    subzstack = subz_from_3dtiff(tif3darray, 50)
    print tif3darray.shape
    print subzstack.shape
    
    # make mask of probe cross-section
    probesizepix = micron2pix(probesizemicron)
    l = max(probesizepix)
    k = np.zeros((l,l))    
    k[l/2] = 1 # horizontal probe
    damage2d = np.zeros((len(theta),
                         ((subzstack.shape[1] - l) / downsample + 1),
                         ((subzstack.shape[2] - l) / downsample + 1)))
    lineprofiles = np.zeros((subzstack.shape[0], len(theta), 
                             ((subzstack.shape[1] - l) / downsample + 1),
                             ((subzstack.shape[1] - l) / downsample + 1)), dtype=object)
    print damage2d.shape
    for lnum,layer in enumerate(subzstack):
        print 'layer '+str(lnum+1)+' of '+str(subzstack.shape[0])
        layer = dsp.normalize(layer)
        for i, t in enumerate(theta): 
            print('%s degrees'', all translations...'%t)
            rotated = misc.imrotate(k, t, interp='nearest')
            # get rid of rotation artifacts
            rotated[rotated < rotated.max()] = 0
            rotated[rotated == rotated.max()] = 1
            # for y in xrange(layer.shape[0]+l): # add step here...
            ys = 0
            for y in xrange(0, layer.shape[0]-l, downsample):
                xs = 0 
                # for x in xrange(layer.shape[1]-l): # add step here...
                for x in xrange(0, layer.shape[1]-l, downsample): 
                    a=rotated*layer[y:y+l,x:x+l]
                    b=a.T[a.T.nonzero()]
                    # luminance est of collision

                    # profl = b
                    profl = dsp.smoothg(b,20) # smoothed profile for this location
                    # if b.mean() > .5:
                    c = profl > 1.5
                    damage2d[i,ys,xs] += np.diff(c).nonzero()[0][::2].size
                    profl[profl > 1.5] = 1.5
                    
                    profl[profl < thresh] = thresh # discard small peaks
                    
                    lineprofiles[lnum, i, ys, xs] = profl
                    
                    # count peaks --> vessels
                    # damage2d[i,y,x] += sum(dsp.islocmax(profl))
                    damage2d[i,ys,xs] += argrelmax(profl, order=5)[0].size 
                    # damage2d[i, ys, xs] += len(find_events(profl, thresh))
                    # damage2d[i, ys, xs] += len(find_peaks_cwt(b, np.arange(5, 20)))
                    xs += 1
                ys += 1

    return damage2d, lineprofiles, subzstack

# ori = np.array([0,30,60,90,120,150])
ori = np.array([0])
damage, lineprof, subz = count_vessels(probesize, ori)
lumdam = build_histogram(probesize, ori)
probe = build_kernel(probesize, ori)
kxy = probe.shape[1]

# plot vascular pattern with optimal locations overlaid
fig = plt.figure(figsize=(10,6))
fig.canvas.mpl_connect('motion_notify_event', mousemove)
fig.canvas.mpl_connect('button_press_event', mouseclick)

# find best points in each region
yregions = int(np.ceil(damage.shape[1] / regionSize))
xregions = int(np.ceil(damage.shape[2] / regionSize))
print damage.shape, xregions, yregions

best = {}
for y in [y*regionSize for y in xrange(yregions)]:
    for x in [x*regionSize for x in xrange(xregions)]:
        print str(x) +' to '+str(x+regionSize)
        print str(y) +' to '+str(y+regionSize)
        region = damage[:, y:y+regionSize, x:x+regionSize]
        rmin = region.min()
        regionmin = np.where(region == rmin)
        # get absolute index and convert to microns
        ys = regionmin[1] + y; xs = regionmin[2] + x
        i = choice(xrange(ys.size))
        print xs[i], ys[i], region.min()
        xm, ym = (xs[i] * downsample * scales[0], 
                  ys[i] * downsample * scales[1])

        theta = ori[regionmin[0][i]]

        best[(x*downsample, y*downsample)] = (rmin, theta, xm, ym)
        
print 'best locations: '

def pretty_best(dic):
    for k in dic.keys():
        print 'region  : '+str(k)
        print 'damage  : '+str(dic[k][0])
        print 'rotation: '+str(dic[k][1])
        print 'locaton : '+str(dic[k][2:])

pretty_best(best)
                              
# # plot the histogram of damage at all translations and rotations
# plt.subplot(131)
# plt.hist(damage.flatten(), bins = 50, alpha=0.25, label='probe size = '+str(probesz))
# plt.ylabel('Count')
# plt.xlabel('Z-collapsed integrated luminance')
# plt.title('Vascular damage vs. probesize: translation and rotation')
# plt.legend()
# 



# plot the estimated damage for the least damaging orientation
ax1 = fig.add_subplot(221)
ax1.imshow(damage.min(0), interpolation='bicubic')
ax1.set_title('Vascular damage at best orientation')

# plot the interactive GUI
ax2 = plt.subplot(222)
ax2im = ax2.imshow(meanvasdam,'gray',interpolation='none')

# plot probe
img = plt.imshow(probe[0], extent = [0,0,0,0], alpha = 0.5)

# plot dynamic timecourse
# ax3 = fig.add_subplot(223)
# ax3.plot([1,2,3,4,5])

# plot histogram with slider
ax3 = fig.add_subplot(223)
ax3.hist(damage.flatten())
slide = plt.Slider(ax3, "", 0, 10)

# plot luminance
ax4 = fig.add_subplot(224)
ax4.imshow(lumdam.min(0))

plt.draw()
plt.show()
plt.ioff()


