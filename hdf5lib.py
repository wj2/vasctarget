"""
Library of functions to save/load in hdf5 binary format
"""
from __future__ import with_statement

import numpy as np
import data as D
import custom_data_tim as DT
import custom_data_hirschlab as DH

try:
    import tables
except:
    print "could not load pytables"

__all__ = ['load_hdf5','save_hdf5']

_recording_classes = [D.Recording,DT.CatRecording,DT.MonkeyRecording,DH.WholeCellRecording]
_recording_class_names = [c.__name__ for c in _recording_classes]
_recording_class = dict(zip(_recording_class_names,_recording_classes))
def empty_instance(cls):
    class Empty(cls):
        def __init__(self): pass
    x = Empty()
    x.__class__ = cls
    return x


def save_hdf5(data,filename):
    """
    saves data into hdf5 binary file
    """
    with tables.openFile(filename,'w') as h5:
        save_data_hdf5(data,h5)

def save_data_hdf5(data,fileh,node=None):
    """
    saves data into node in hdf5 binary file
    """
    if node == None: node = fileh.root
    node._v_attrs.TITLE = data.__class__.__name__
    if data == None:
        save_none_hdf5(data,fileh,node)
    elif type(data) in [float,int,str,np.double]:
        save_scalar_hdf5(data,fileh,node)
    elif type(data) in [list,tuple,set]:
        save_iterable_hdf5(data,fileh,node)
    elif type(data) in [dict,D.Info,D.Container]:
        save_dict_hdf5(data,fileh,node)
    elif type(data) in [np.ndarray,np.matrix]:
        save_array_hdf5(data,fileh,node)
    elif isinstance(data,D.Events) or isinstance(data,D.Stimulus):
        save_events_hdf5(data,fileh,node)
    elif isinstance(data,D.TimeSeries) or isinstance(data,D.Image):
        save_info_array_hdf5(data,fileh,node)
    elif isinstance(data,D.Epochs):
        save_epochs_hdf5(data,fileh,node)
    elif isinstance(data,D.Recording):
        save_recording_hdf5(data,fileh,node)
    else:
        raise NotImplementedError


def load_hdf5(filename_or_node):
    """
    loads data classes from hdf5 binary file
    """
    if type(filename_or_node) == str:
        with tables.openFile(filename_or_node,'r') as h5:
            return load_hdf5(h5.root)
    elif isinstance(filename_or_node,tables.Node):
        node = filename_or_node
        title = node._v_title
        if title == 'NoneType':
            return None
        elif title in ['int','float','str','float64']:
            return load_scalar_hdf5(node)
        elif title in ['list','tuple','set']:
            return load_iterable_hdf5(node)
        elif title in ['dict','Info','Container']:
            return load_dict_hdf5(node)
        elif title in ['ndarray','matrix','alist','atuple','aset']:
            return load_array_hdf5(node)
        elif title in ['Events','Stimulus']:
            return load_events_hdf5(node)
        elif title == 'TimeSeries':
            return load_time_series_hdf5(node)
        elif title in ['Epoch','Epochs']:
            return load_epochs_hdf5(node)
        elif title  == 'Image':
            return load_info_array_hdf5(node)
        elif title  in _recording_class_names:
            return load_recording_hdf5(node)
        else:
            raise NotImplementedError
        return obj
    else:
        assert 0, 'hdf5 file not recognized'


def save_none_hdf5(data,fileh,node=None):
    """
    saves None into node in hdf5 file
    """
    if node == None: node = fileh.root
    node._v_attrs.none = None

def save_scalar_hdf5(data,fileh,node=None):
    """
    saves scalar into node in hdf5 file
    """
    if node == None: node = fileh.root
    node._v_attrs.scalar = data

def save_iterable_hdf5(data,fileh,node=None):
    """
    saves iterable into node in hdf5 file
    """
    if node == None: node = fileh.root
    typeset = set(map(type,data))
    if len(typeset) == 1 and typeset.pop().__name__ in ['int','float','float64']:
        node._v_attrs.TITLE = 'a'+node._v_attrs.TITLE
        save_array_hdf5(list(data),fileh,node=node)
    else:
        for i,val in enumerate(data):
            vnode = fileh.createGroup(node, "value%d"%i)
            save_data_hdf5(val,fileh,node=vnode)

def save_dict_hdf5(data,fileh,node=None):
    """
    saves dictionary into node in hdf5 file
    """
    if node == None: node = fileh.root
    for i,(key,val) in enumerate(data.items()):
        vnode = fileh.createGroup(node, "value%d"%i)
        knode = fileh.createGroup(node, "key%d"%i)
        save_data_hdf5(key,fileh,node=knode)
        save_data_hdf5(val,fileh,node=vnode)

def save_info_as_attr_hdf5(data,fileh,node=None):
    """
    saves info dictionary into attributes of node in hdf5 file
    """
    if node == None: node = fileh.root
    for key,value in data.iteritems():
        node._v_attrs.__setattr__(key,value)

def save_events_hdf5(data,fileh,node=None):
    """
    saves events into node of hdf5 file
    """
    if node == None: node = fileh.root
    fileh.createArray(node, "ndarray", data)
    infnode = fileh.createGroup(node, "info")
    save_data_hdf5(data.i,fileh,infnode)
    valnode = fileh.createGroup(node, "values")
    for key in data.keys():
        fileh.createArray(valnode, key, getattr(data,key))

def save_info_array_hdf5(data,fileh,node=None):
    """
    saves info array into node in hdf5 file
    """
    if node == None: node = fileh.root
    fileh.createArray(node, "ndarray", data)
    infnode = fileh.createGroup(node, "info")
    save_data_hdf5(data.i,fileh,infnode)

def save_array_hdf5(data,fileh,node=None):
    """
    saves info array into node in hdf5 file
    """
    title = node._v_attrs.TITLE
    if node == None: node = fileh.root
    fileh.createArray(node, "ndarray", data)

def save_epochs_hdf5(data,fileh,node=None):
    """
    saves epochs into node in hdf5 file
    """
    if node == None: node = fileh.root
    if np.iterable(data):
        node._v_attrs.TITLE = 'Epochs'
        fileh.createArray(node, "ndarray", np.array(zip(data['tstart'],data['duration'],data['t0_offset'])))
    else:
        node._v_attrs.TITLE = 'Epoch'
        fileh.createArray(node, "ndarray", np.array((data['tstart'],data['duration'],data['t0_offset'])))

def save_recording_hdf5(data,fileh,node=None):
    """
    saves recording into node of hdf5 file
    """
    if node == None: node = fileh.root
    for key,val in data.__dict__.iteritems():
        valnode = fileh.createGroup(node, key)
        save_data_hdf5(val,fileh,valnode)

def load_scalar_hdf5(node):
    """
    loads scalar from node in hdf5 file
    """
    return node._v_attrs.scalar

def load_iterable_hdf5(node):
    """
    loads iterable from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    nodedict = node._v_children
    iterable = (load_hdf5(nodedict["value%d"%i]) for i in range(node._v_nchildren))
    if title == 'list':
        return list(iterable)
    elif title == 'tuple':
        return tuple(iterable)
    elif title == 'set':
        return set(iterable)
    else:
        raise NotImplementedError

def load_dict_hdf5(node):
    """
    loads dictionary from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    nodedict = node._v_children
    iterkey = (load_hdf5(nodedict["key%d"%i]) for i in range(node._v_nchildren/2))
    iterval = (load_hdf5(nodedict["value%d"%i]) for i in range(node._v_nchildren/2))
    if title == 'dict':
        return dict(zip(iterkey,iterval))
    elif title == 'Info':
        return D.Info(zip(iterkey,iterval))
    elif title == 'Container':
        return D.Container(zip(iterkey,iterval))
    else:
        raise NotImplementedError

def load_info_as_attr_hdf5(node):
    """
    loads info dictionary into attributes of node in hdf5 file
    """
    attrs = node._v_attrs
    info = D.Info()
    for key in attrs._f_list('user'):
        info[key] = attrs.__getattr__(key)
    return info

def load_events_hdf5(node):
    """
    loads events from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    info = load_hdf5(node.info)
    events = D.Events(node.ndarray.read(),i=info)
    for val in node.values._f_iterNodes():
        setattr(events,val._v_name,val.read())
    if title == 'Events':
        return events
    elif title == 'Stimulus':
        return D.Stimulus(events)
    else:
        raise NotImplementedError

def load_info_array_hdf5(node):
    """
    loads info array from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    info = load_hdf5(node.info)
    infoarray = node.ndarray.read()
    if title == 'Image':
        return D.Image(infoarray,i=info)
    else:
        raise NotImplementedError

def load_array_hdf5(node):
    """
    loads info array from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    arr = node.ndarray.read()
    if title == 'ndarray':
        return arr
    if title == 'matrix':
        return np.matrix(arr)
    if title == 'alist':
        return list(arr)
    if title == 'atuple':
        return tuple(arr)
    if title == 'aset':
        return set(arr)
    else:
        raise NotImplementedError

def load_time_series_hdf5(node):
    """
    loads time series from node in hdf5 file
    """
    info = load_hdf5(node.info)
    timeseries = D.TimeSeries(node.ndarray.read(),samplingrate=info.samplingrate,i=info)
    return timeseries

def load_epochs_hdf5(node):
    """
    loads epochs from node in hdf5 file
    """
    epochs = node.ndarray.read()
    title = node._v_attrs.TITLE
    if title == 'Epoch':
        return D.Epochs(tuple(epochs))
    elif title == 'Epochs':
        return D.Epochs([tuple(epoch) for epoch in epochs])
    else:
        raise NotImplementedError

def load_recording_hdf5(node):
    """
    loads recording from node in hdf5 file
    """
    title = node._v_attrs.TITLE
    rec = empty_instance(_recording_class[title])
    for vnode in node._f_iterNodes():
        setattr(rec,vnode._v_name,load_hdf5(vnode))
    return rec
