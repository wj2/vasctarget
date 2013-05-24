import os
import glob

import numpy as np
import Image

import xml_parser as xp


class tiffile:

    def __init__(self,path, xmlpath):
        self.path = path
        self.xmlpath = xmlpath

    def gettiflist(self,path):
        tiflist = os.listdir(path)  # get the list of all the file names in the directory
        
        return tiflist        

    def get_tiffs_and_xml(self,tiflist):
        
        for idx,infile in enumerate(tiflist):# purify the list only with tiff files
            name,ext = os.path.splitext(infile)
            if ext == '.xml':
                xml = infile
            elif ext != '.tif':
                tiflist.remove(infile)

        for idx,infile in enumerate(tiflist):# WHY???
            name,ext = os.path.splitext(infile)
            if ext != '.tif':
                tiflist.remove(infile)

        return tiflist, xml


    def listsize(self,xmlpath):    # Z-size of array
        
        num_list = xp.countframe(xmlpath)

        return num_list


    def initarray(self,num_list,xmlpath):

        pixels = xp.getpixels(xmlpath)
        
        x = pixels[0]
        y = pixels[1]

        array3d = np.empty([num_list, x, y]) #create an empty array

        return array3d
    

    def tif2array(self,path,xmlpath,tiflist,array3d):

        pixels = xp.getpixels(xmlpath)
        
        for idx,infile in enumerate(tiflist):
            im = Image.open(os.path.join(path,infile))
            A = list(im.getdata())
            tifarray = np.array(A,np.uint16).reshape(pixels[0],pixels[1])
            array3d[idx,:,:] = tifarray # fill the empty array with values
            
        return array3d



class vasdam:

    def  __init__(self,probsize):
        self.probsize = probsize


    def damcalc(self,probsize,meanvasdam):
        m = probsize[0]
        n = probsize[1]

        S = meanvasdam.shape
        p = S[0] +1 - n
        q = S[1] + 1 - m
        damage = np.zeros([p,q])

        for i in range(p):
            for j in range(q):
                partif = meanvasdam[i:n+i,j:m+j]
                damage[i][j] = partif.sum()

        return damage




    def maxmindam(self,probsize,meanvasdam):

        m = probsize[0]
        n = probsize[1]

        damage = vasdam.damcalc(self,probsize,meanvasdam)

        maxdam = np.argmax(damage)
        mindam = np.argmin(damage)

        maxdamval = np.max(damage)
        mindamval = np.min(damage)
        
        dimsdam = damage.shape
        
        maxidx = np.unravel_index(maxdam,dimsdam)
        minidx = np.unravel_index(mindam,dimsdam)

        maxvasdam = np.array([maxidx[0],m+maxidx[0],maxidx[1],n+maxidx[1]]).reshape([2,2])
        minvasdam = np.array([minidx[0],m+minidx[0],minidx[1],n+minidx[1]]).reshape([2,2])

        location = np.array([maxvasdam,int(maxdamval),minvasdam,int(mindamval)])
        

        return  location





