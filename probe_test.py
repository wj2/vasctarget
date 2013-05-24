import vastifarr as va
import numpy as np

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

probsize1 = micron2pix([15,25])
probsize2 = micron2pix([15,50])
probsize3 = micron2pix([15,100])


b = va.vasdam(probsize1)
c = va.vasdam(probsize2)
d = va.vasdam(probsize3)

damage1 = b.damcalc(probsize1,meanvasdam)
damage2 = c.damcalc(probsize2,meanvasdam)
damage3 = d.damcalc(probsize3,meanvasdam)

location1 = b.maxmindam(probsize1,meanvasdam)
location2 = c.maxmindam(probsize2,meanvasdam)
location3 = d.maxmindam(probsize3,meanvasdam)

print location1, location2, location3

import matplotlib.pyplot as plt
import matplotlib.cm as cm


def disp_maxmin(location):      #display of maximum and minimum damage on z-collapsed 2d array
    ymaxlow = location[0][0][0]
    ymaxhi = location[0][0][1]
    xmaxlow = location[0][1][0]
    xmaxhi = location[0][1][1]

    yminlow = location[2][0][0]
    yminhi = location[2][0][1]
    xminlow = location[2][1][0]
    xminhi = location[2][1][1]

    annotmax = 'maximum ='+ str(location[1])
    annotmin = 'minimum ='+ str(location[3])

    plt.imshow(meanvasdam, cmap = cm.Greys_r)
    plt.axhspan(yminlow, yminhi,xmin=1.*xminlow/512, xmax=1.*xminhi/512, facecolor='g', edgecolor='w',alpha=0.5)
    plt.axhspan(ymaxlow, ymaxhi,xmin=1.*xmaxlow/512, xmax=1.*xmaxhi/512, facecolor='r', alpha=0.5)
    plt.xlabel('pixels(x)')
    plt.ylabel('pixels(y)')
    plt.ylim([512,0])
    
    plt.annotate('maximum ='+ str(location[1]),xy= (xmaxhi,ymaxhi),  xycoords='data')
    plt.annotate('minimum ='+ str(location[3]),xy= (xminhi,yminhi),  xycoords='data')
    plt.show()

disp_maxmin(location1)
disp_maxmin(location2)
disp_maxmin(location3)

    
