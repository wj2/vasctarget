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

import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.imshow(meanvasdam, cmap = cm.Greys_r)
plt.show()



        
