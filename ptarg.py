#!/usr/bin/python

import numpy as np, dsp.dsplib as dsp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import misc
from scipy.signal import argrelmax

# replace this section with a few lines
import vasc.vastifarr as va

from config import path, probesize

a = va.tiffile(path)
tiflist = a.gettiflist(path)
tiflist2, xml = a.get_tiffs_and_xml(tiflist)
tiflist2 = tiflist2[50:75]
nfiles = a.listsize(tiflist2)
array3d = a.initarray(nfiles)
tif3darray = a.tif2array(path,tiflist2,array3d)
###########################################

meanvasdam = np.mean(tif3darray,axis=0)

def mousemove(event):
    # plot best orientation at a given location
    if event.inaxes == ax1 or event.inaxes == ax2:
        i = np.argmin(damage[:,event.ydata,event.xdata])
        img.set_data(probe[i])
        extent = [event.xdata - kxy/2, 
                  event.xdata + kxy/2, 
                  event.ydata - kxy/2, 
                  event.ydata + kxy/2]
        img.set_extent(extent)
        ax2.set_ylim(512,0)
        ax2.set_xlim(0,512)
        print('damage=%1d'%np.min(damage[:,event.ydata,event.xdata]))
        
        ax3.clear()
        ax3.plot(lineprof[0, i, event.ydata, event.xdata])

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


def micron2pix(probsize):       
    scales = a.get_scale(xml)
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
        holder[i] = np.mean(tiffs[i*slicesize:i*slicesize + slicesize], 
                            axis=0)
        i += 1
    
    return holder

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
    thresh = -0.9
    # create sub-z tiff stack
    subzstack = subz_from_3dtiff(tif3darray, 25)
    print tif3darray.shape
    print subzstack.shape
    
    # make mask of probe cross-section
    probesizepix = micron2pix(probesizemicron)
    l = max(probesizepix)
    k = np.zeros((l,l))    
    k[l/2] = 1 # horizontal probe
    damage2d = np.zeros((len(theta),meanvasdam.shape[0],
                         meanvasdam.shape[1]))
    lineprofiles = np.zeros((subzstack.shape[0], len(theta), meanvasdam.shape[0],
                             meanvasdam.shape[1]), dtype=object)
    for lnum,layer in enumerate(subzstack):
        print 'layer '+str(lnum+1)+' of '+str(subzstack.shape[0])
        layer = dsp.normalize(layer)
        for i, t in enumerate(theta): 
            print('%s degrees'', all translations...'%t)
            # TODO : find better rotate function
            # rotated = ndimage.rotate(k,t,reshape=False)
            rotated = misc.imrotate(k, t, interp='nearest')
            # get rid of rotation artifacts
            rotated[rotated < rotated.max()] = 0
            rotated[rotated == rotated.max()] = 1
            # for y in xrange(layer.shape[0]+l): # add step here...
            for y in xrange(layer.shape[0]/2-100, layer.shape[0]/2+100):
                # for x in xrange(layer.shape[1]-l): # add step here...
                for x in xrange(layer.shape[1]/2-100, layer.shape[1]/2+100): 
                    # print lnum, x, y
                    a=rotated*layer[y:y+l,x:x+l]
                    b=a.T[a.T.nonzero()]
                    # luminance est of collision
                    if b.mean() > 0:
                        damage2d[i,y,x] += 1
                    # profl = b
                    profl = dsp.smoothg(b,9) # smoothed profile for this location
                    profl[profl < thresh] = thresh # discard small peaks
                    
                    lineprofiles[lnum, i, y, x] = profl

                    # count peaks --> vessels
                    # damage2d[i,y,x] += sum(dsp.islocmax(profl))
                    damage2d[i,y,x] += argrelmax(profl, order=1)[0].size 

    return damage2d, lineprofiles

ori = np.array([0,30,60,90,120,150])
damage, lineprof = count_vessels(probesize, ori)
lumdam = build_histogram(probesize, ori)
probe = build_kernel(probesize, ori)
kxy = probe.shape[1]

# plot vascular pattern with optimal locations overlaid
fig = plt.figure(figsize=(10,6))
fig.canvas.mpl_connect('motion_notify_event',mousemove)

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
ax1.imshow(damage.min(0))
ax1.set_title('Vascular damage at best orientation')

# plot the interactive GUI
ax2 = plt.subplot(222)
ax2.imshow(meanvasdam,'gray',interpolation='none')

# plot probe
img = plt.imshow(probe[0], extent = [0,0,0,0], alpha = 0.5)

# plot dynamic timecourse
ax3 = fig.add_subplot(223)
ax3.plot([1,2,3,4,5])

# plot luminance
ax4 = fig.add_subplot(224)
ax4.imshow(lumdam.min(0))

plt.draw()
plt.show()
plt.ioff()


