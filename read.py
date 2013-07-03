
import os

import numpy as np
import BeautifulSoup as bs

from scipy.misc import imread
from PIL import Image

from tifffile import TIFFfile
from proc import collapse_stack


def find_tifs_and_xml_folder(path, channel='Ch1'):
    files = os.listdir(path)
    tifseries = []
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.xml':
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
    widths = soup('key', key='positionCurrent_ZAxis')
    info['z_width'] = abs(float(widths[1]['value']) 
                          - float(widths[0]['value']))

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
              
def get_data(path, zthick):
    if os.path.isdir(path):
        xml, series = find_tifs_and_xml_folder(path)
        data = stack_tifs(series, info['xsize'], info['ysize'])
    elif os.path.isfile(path):
        xml, data = find_tifs_and_xml_file(path)
    else:
        sys.stderr.write('InputError: path given does not exist')
        sys.exit(1)

    info = get_config_info(xml)    
    data = collapse_stack(data, info['z_width'], zthick)
    return info, data

