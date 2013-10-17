
import os, sys

import numpy as np
import BeautifulSoup as bs
# import sys

from scipy.misc import imread
from PIL import Image

from tifffile import TIFFfile
from proc import collapse_stack


def find_tifs_and_xml_folder(path, channel='Ch1'):
    files = os.listdir(path)
    tifseries = []
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.xml' or ext == '.cfg':
            xml = os.path.join(path, f)
        elif ext == '.tif' and channel in f:
            tifseries.append(os.path.join(path, f))
    return xml, sorted(tifseries)

def find_tifs_and_xml_file(path):
    data = TIFFfile(path).asarray()
    name = os.path.splitext(path)[0]
    xml = name + '.xml'

    return xml, data

def get_config_info(xml):
    soup = bs.BeautifulStoneSoup(open(xml, 'rb'))
    info = {}
    info['xsize'] = int(soup.find('key', 
                                    key='pixelsPerLine')['value'])
    info['ysize'] = int(soup.find('key', 
                                    key='linesPerFrame')['value'])
    info['xmpp'] = float(soup.find('key', 
                                   key='micronsPerPixel_XAxis')['value'])
    info['ympp'] = float(soup.find('key', 
                                   key='micronsPerPixel_YAxis')['value'])
    if os.path.splitext(xml)[1] == '.xml':
        widths = soup('key', key='positionCurrent_ZAxis')
        info['z_width'] = abs(float(widths[1]['value']) 
                              - float(widths[0]['value']))
    # elif os.path.splitext(xml)[1] == '.cfg':
    #     width = float(soup('key', key='motorStepSize_ZAxis')[0]['value'])
    #     info['z_width'] = abs(width)
    else:
        info['z_width'] = 1

    return info

def stack_tifs(tifseries, xsize, ysize):
    """
    Assumes tif series is in order
    """
    holster = np.empty((len(tifseries), ysize, xsize))
    for i, layer in enumerate(tifseries):
        lay = Image.open(layer)
        lis = list(lay.getdata())
        holster[i] = np.array(lis, np.uint16).reshape(512, 512)

    return holster
              
def get_data(path, zthick, chan):
    if os.path.isdir(path):
        xml, series = find_tifs_and_xml_folder(path, chan)
        sys.stdout.write('in stack: '+str(len(series))+'\n')
        info = get_config_info(xml)        
        data = stack_tifs(series, info['xsize'], info['ysize'])
    elif os.path.isfile(path):
        xml, data = find_tifs_and_xml_file(path)
        info = get_config_info(xml)    
    else:
        sys.stderr.write('InputError: path given does not exist\n')
        sys.exit(1)

    data = collapse_stack(data, info['z_width'], zthick)
    return info, data

