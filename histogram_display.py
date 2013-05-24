import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import vastifarr as va

path = '/users/april/desktop/Research/Z-Stack_PRE-000/'

a = va.tiffile(path)
tiflist = a.gettiflist(path)
tiflist2 = a.puretiflist(tiflist)
num_list = a.listsize(tiflist2)
array3d = a.initarray(num_list)
tif3darray = a.tif2array(path,tiflist2,array3d)

meanvasdam= np.mean(tif3darray,axis=0)

def micron2pix(probsize):       # conversion from microns to pixels
    s  = 1/0.909090909090908
    newprobsize = [int(probsize[0]*s),int(probsize[1]*s)]
    return newprobsize

micron1 = [15,25]
micron2 = [15,50]
micron3 = [15,100]
probsize1 = micron2pix(micron1)
probsize2 = micron2pix(micron2)
probsize3 = micron2pix(micron3)

def vashist(probsize):
    b = va.vasdam(probsize)

    damage = b.damcalc(probsize,meanvasdam)
    damage1d = damage.reshape(-1,order='F')
    damagesize = damage.size

    meandam = np.mean(damage1d)
    stddam = np.std(damage1d)  # how do I display mean and std as a second line in title?

    # the histogram of the data
    n, bins, patches = plt.hist(damage1d,150, facecolor='green')

    plt.ylabel('Count')
    plt.xlabel('Summated Vascular Damage over a Given Probesize')
    plt.title("Vascular Damage with a probe size ="+ str(probsize))

    plt.grid(True)

    plt.show()


vashist(probsize1)


vashist(probsize2)


vashist(probsize3)


