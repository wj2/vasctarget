import xml.etree.ElementTree as ET
import BeautifulSoup as bs
        

def getsequence(xmlpath):

    tree = ET.parse(xmlpath)
    root = tree.getroot()
    sequence = root[1]

    return sequence



def getpvss(xmlpath):
    
    tree = ET.parse(xmlpath)
    root = tree.getroot()
    sequence = root[1]
    pvss = sequence[0][1]

    return pvss



def countframe(xmlpath):
    count = 0

    tree = ET.parse(xmlpath)
    root = tree.getroot()
    sequence = root[1]
    
    for frame in sequence.iter('Frame'):
        count = count + 1

    return count


def getpixels(xmlpath):

    tree = ET.parse(xmlpath)
    root = tree.getroot()
    sequence = root[1]
    pvss = sequence[0][1]

    for Key in pvss.findall('Key'):
        if Key.get('key') == 'pixelsPerLine':
            ppl = Key.get('value')
            x = int(ppl)

    for Key in pvss.findall('Key'):
        if Key.get('key') == 'linesPerFrame':
            lpf = Key.get('value')
            y = int(lpf)

    return [x,y]


def getscale(xmlpath):
        
    soup = bs.BeautifulStoneSoup(open(xmlpath))
    xscale = float(soup.find('key',key="micronsPerPixel_XAxis")['value'])
    yscale = float(soup.find('key',key="micronsPerPixel_YAxis")['value'])

    return (xscale, yscale)

		


