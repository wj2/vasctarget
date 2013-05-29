import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import vastifarr as va

path = '/Users/MacbookAir/Dropbox/Blanche_Lab/vasculature/Sample data/Z-Stack_PRE-000'

a = va.tiffile(path)
tiflist = a.gettiflist(path)
tiflist2 = a.puretiflist(tiflist)#[100:]
num_list = a.listsize(tiflist2)
array3d = a.initarray(num_list)
tif3darray = a.tif2array(path,tiflist2,array3d)

meanvasdam= np.mean(tif3darray,axis=0)

def micron2pix(probsize):       # conversion from microns to pixels
    s  = 1/0.909090909090908
    newprobsize = [int(probsize[0]*s),int(probsize[1]*s)]
    return newprobsize

micron0 = [7,7]
micron1 = [23,25]
micron2 = [23,75]
micron3 = [15,205]

def vashist(probsizemicron):
    probsize = micron2pix(probsizemicron)
    b = va.vasdam(probsize)

    damage = b.damcalc(probsize,meanvasdam)
    damage1d = damage.reshape(-1,order='F')
    damagesize = damage.size

    meandam = np.mean(damage1d)
    stddam = np.std(damage1d)  # how do I display mean and std as a second line in title?

    # the histogram of the data
    # probarea = probsize[0]*probsize[1]
    plt.hist(damage1d,np.arange(0,3e6,5e4), alpha=0.25, label=str(probsizemicron))

    plt.ylabel('Count')
    plt.xlabel('Z-collapsed integrated luminance')
    plt.title('Vascular damage vs. probesize: translation only')
    plt.legend()

    plt.grid(True)

    plt.show()

# vashist(micron0)

vashist(micron1)

# vashist(micron2)
# 
# vashist(micron3)


