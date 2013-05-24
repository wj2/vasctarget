"""
Library of functions related to data.py
"""
import numpy as np
import data
import utils

try:
    import scipy.fftpack as F
except:
    import numpy.fft as F
import sys
from itertools import izip

__all__ = ['load_events','load_event_timings','load_time_series','load_movie','crop',
           'cut_time_series','cut_events','events_to_time_series','buildSweepTable',
           'events_from_trials','epochs_from_onsets','rescale_event_time']
_eps = np.MachAr().eps


def raw_to_SI(dat,scalefactor,offset,overwrite=True):
    """
    converts raw binary-offset data into SI units
    """
    if overwrite:
        result = dat
    else:
        result = dat.copy()
    result -= offset
    result *= scalefactor
    return result


def has_trials(obj):
    "returns the number of trials, or zero if none"
    if type(obj) in [list,data.Epochs,data.Events]:
        if hasattr(obj, 'trials'):
            return obj.i.ntrials
        else:
            return 0
    elif type(obj) in [np.array,data.TimeSeries]:
        if obj.ndim == 1:
            return 0
        elif obj.ndim == 2:
            return obj.shape[0]
        else: assert 0, 'has_trials is not defined for %dd array'%obj.ndim


def load_time_series(filename,samplingrate,dtype,nsamples=None,byteorder='little'):
    """
    Load time-series data from file.

    output data type is float ('d')
    """
    dat = np.fromfile(filename,dtype=dtype)
    if byteorder != sys.byteorder:
        dat = dat.byteswap()
    if (nsamples != None) & (dat.shape[0] != nsamples):
        assert 0, 'wrong number of samples in %s: expected %d, read %d'%(filename,nsamples,dat.shape[0])
    return data.TimeSeries(dat.astype('d'),samplingrate)


def load_event_timings(filename,dtype,byteorder='little'):
    """
    Load time-stamps for events from a file.
    """
    dat = np.fromfile(filename,dtype=dtype)
    if byteorder != sys.byteorder:
        dat = dat.byteswap()
    return data.Events(dat.astype('d'))


def load_events(filename,dtype_t,dtype_v,byteorder='little',cast_value=None):
    """
    Load time-stamps and values from file.
    """
    if dtype_t != dtype_v:
        assert 0, 'Different data type for time stamps and values not implemented.'
    dat = np.fromfile(filename,dtype=dtype_t)
    if byteorder != sys.byteorder:
        dat = dat.byteswap()
    time = dat[::2].astype('d')
    if cast_value:
        value = dat[1::2].astype(cast_value)
    else:
        value = dat[1::2]
    return data.Events(time,v=value)


def load_movie(filename,dtype,nx,ny,nframes,hdrlen=0,byteorder='little'):
    """
    Loads a non-compressed movie from file.
    """
    dat = np.fromfile(filename,dtype=dtype)
    if byteorder != sys.byteorder:
        dat = dat.byteswap()
    dat = dat[hdrlen:] # skip over movie header, if present
    try:  # reshape into a movie
        dat = dat.reshape(nframes, nx, ny)
    except ValueError: # movie is not the expected size
        assert 0, 'Movie file incorrect size for dimensions provided.'

    return data.Image(dat)

def events_from_trials(dat, i=None, dtype=None, copy=False, **valuedict):
    cdata = np.concatenate(dat)
    cdict = {}
    for key,val in valuedict.iteritems():
        cdict[key] = np.concatenate(val)
    cdict['trials'] = np.concatenate([ii*np.ones(l,int) for ii,l in enumerate([len(d) for d in dat])])
    obj = data.Events(cdata,i=i,dtype=dtype,copy=copy,**cdict)
    obj.i.ntrials = len(dat)
    return obj

def epochs_from_onsets(onsets,duration,offset=0):
    """
    Returns epochs of fixed duration with start times given by onsets minus offset.
    """
    return data.Epochs([(float(t)-offset,duration,offset) for t in onsets], 'trials')


def cut_time_series(dat,epochs):
    """
    Partitions input array 'timeseries' according to 'epochs'.

    Partitions an array of continuously sampled 'data' according to 'epochs'.
    Epochs specifies start/end time pairs or start/duration/offset triples
    that 'cut' the timeseries (1d) into a 2d array of timeseries.

    If the timeseries are already cut (ie. the dimensionality of data is 2d)
    then only a single epoch can be specified, and the timeseries returned
    is a subset of the original timeseries array of the same length.
    """
    samplingrate = dat.i.samplingrate
    # if dat.has_trials():
    #     assert 0, 'not implemented yet'

    n = dat.shape[-1] # number of sample points
    if type(epochs) != data.Epochs:
        epochs = data.Epochs(epochs)

    if epochs.ndim == 0: # single epoch, don't add another dimension
        obj = dat[...,epochs.asslice(samplingrate,nmax=n,t0=dat.i.t0,include_last_bin=True)]
        t0 = -epochs['t0_offset']
    else:
        if len(epochs) > 0: # multiple epochs, add another dimension to data
            t0 = -epochs[0]['t0_offset']
            slices = epochs.asslice(samplingrate,nmax=n,t0=dat.i.t0,include_last_bin=True)
            nt = dat[...,slices[0]].shape[-1]
            ne = len(slices) # number of epochs
            newshape = (ne,)+dat.shape[:-1]+(nt,)
            obj = np.empty(newshape,dat.dtype) # allocate memory
            for i in xrange(ne):
                try:
                    obj[i,...,:] = dat[...,slices[i]]
                except:
                    assert (dat[...,slices[i]].shape == newshape[1:]), 'all epochs have to be of the same duration, got %s instead of %s'%(
                        dat[...,slices[i]].shape, newshape[1:])
                    assert 0, 'error cutting TimeSeries'
        else: # no epochs
            t0 = 0
            obj = np.array([])
    return data.TimeSeries(obj,samplingrate,i=dat.i,t0=t0)


def cut_events(events, epochs): #align=True
    """
    Partitions input array 'events' according to 'epochs'.

    Partitions an nd array of 'events' according to 'epochs'. Epochs
    specifies start/end time pairs or start/duration/offset triples
    that 'cut' the events (1d) into a 2d array of events.

    If the events are already cut (ie. the dimensionality of events is 2d)
    then only a single epoch can be specified, and the event array returned
    is a subset of the original event array of the same length.

    TO DO: if align = True (default), epoch start times are subtracted from
    event time-stamps, if False, the event time-stamps are returned unchanged.

    Optional argument: align = True.
    """
    if type(epochs) != data.Epochs:
        epochs = data.Epochs(epochs)

    name = getattr(epochs, 'name', 'idx0')
    # use the first available idx that is not used yet
    if name == 'idx0':
        idx = 0
        while name in events.keys():
            idx += 1
            name = 'idx%d'%i

    if epochs.ndim == 0: # single epoch
        inepoch = ((events >= epochs['tstart']) & (events < epochs['tstop']-10*_eps)) # 2*_eps might be enough
        obj = events[inepoch] - (epochs['tstart']+epochs['t0_offset'])
    else: # multiple epochs, create new idx values
        idx = []
        timestamps = []
        valuedict = dict.fromkeys(events.keys(),[])
        for i,epoch in enumerate(epochs):
            ev = events[epoch]
            timestamps.extend(ev.view(np.ndarray))
            idx.extend([i]*len(ev))
            for k in events.keys():
                valuedict[k].extend(ev.__dict__[k])
        valuedict[name] = idx
        obj = data.Events(timestamps,i=events.i,**valuedict)
        if name == 'trials': obj.i.ntrials = i+1 # new way to store ntrials for has_trials() method
    return obj

def events_to_time_series(evt,samplingrate,nbins=None,binary=False,weights=None,width=None,dtype='d',t0=None):
    """
    Returns a time series sampled with samplingrate that contains the number of events

    If binary == True, multiple events are represented by 1 instead of the number of events.
    If weights (time series) are given, the events are multiplied by the corresponding values of weights.
    If nbins is given, it determines the length of the raster.
    If t0<0 (t0>0) is given, sample points are added (removed) at the beginning of the raster in order
    to make the first sample point to coincide with t0. If no t0 is given, it is assumed to be 0 if all
    events have positive times and the time of the first event otherwise.
    If width is given, it determines the width of an event (stdev of gaussian smoothing window in s).
    """
    t_shift = -min(0,evt.min()) # time shift needed to make events positive
    if t0 is None: t0 = min(0,evt.min())
    t_adjust = t0-min(0,evt.min()) # time shift needed to make the first bin t0
    addn = int(max(0,-int((t_adjust-utils.eps)*samplingrate)))
    cutn = int(max(0,int((t_adjust+utils.eps)*samplingrate)))
    if nbins == None:
        # if evt.has_trials(): nbins = np.array(evt.asidx(samplingrate)[1],'i').max()+addn-cutn+1
        # else: nbins = np.array(evt.asidx(samplingrate),'i').max()+addn-cutn+1
        nbins = int((evt.max()+t_shift+utils.eps)*samplingrate)+addn-cutn+1
    # tmax = float(nbins+cutn-addn)/samplingrate
    tmax = t0+float(nbins)/samplingrate
    idx = (evt[evt<tmax]+t_shift).asidx(samplingrate)
    if not evt.has_trials():
        out = data.TimeSeries(np.zeros(nbins+cutn,dtype),samplingrate=samplingrate,t0=t0)
        if binary:
            if weights == None:
                out[addn:][idx] = 1
            else:
                out[addn:][idx] = weights
        else:
            if weights == None:
                for i in idx: out[addn+i] += 1
            else:
                for wi,i in enumerate(idx): out[addn+i] += weights[wi]
        out = out[cutn:]
    else:
        out = data.TimeSeries(np.zeros((evt.has_trials(),nbins+cutn),dtype),samplingrate=samplingrate,t0=t0)
        if binary:
            out[:,addn:][idx] = 1
        else:
            if weights==None:
                for i in apply(izip,idx): out[:,addn:][i] += 1
            else:
                for wi,i in enumerate(apply(izip,idx)): out[:,addn:][i] += weights[wi]
        out = out[:,cutn:]
    if width: out.smooth(width)
    return out

def rescale_event_time(evt,scale):
    """
    Rescale the timings of events using the scale time series.
    Produces constant event rate if event rate is used as the scale time
    series.
    """

    samplingrate = scale.i.samplingrate
    eigentime = np.zeros_like(scale)
    if scale.has_trials():
        dt = data.TimeSeries(np.cumsum(scale,1)/samplingrate/scale.mean()-scale.time(),samplingrate)
    else:
        dt = data.TimeSeries(np.cumsum(scale)/samplingrate/scale.mean()-scale.time(),samplingrate)

    tmax = scale.time().max()
    evt = evt[evt<tmax]

    tmg = []
    if evt.has_trials():
        for tr in range(evt.has_trials() or 1):
            e_tr = evt.select_dim(tr)
            ind = np.array(e_tr.asidx(samplingrate))
            if scale.has_trials():
                tmg.append((e_tr+dt[tr,:][e_tr]).tolist())
            else:
                tmg.append((e_tr+dt[e_tr]).tolist())
    else:
        tmg.extend(evt+dt[evt])
    return data.Events(tmg)


def buildSweepTable(stim):
    """
    TODO: handle co-varying dimensions (ie. multiple stimulus attributes in one dimension)
          add orioff here or later? add equivalent offsets for other parameters here?
    """
    n=[] # list of number of conditions per dimension
    dnames=[] # list of stimulus dimension names
    for i, dim in enumerate(stim.i.dimlist):
        dnames.append(dim['vars'][0])
        n.append(len(stim.i[dnames[i]]))
    ncond = np.multiply.reduce(n) # number of unique conditions

    # add trials attribute (assumes an equal number of trials per condition)
    if stim.i.nrepeats > 1:
        frms_per_sweep = round(stim.i.sweeptimeSec*stim.i.refreshrate)
        tr = (np.arange(stim.i.nrepeats)+1).repeat((ncond+1)*frms_per_sweep)
        tr = tr[:len(stim)] # crop needed as blank sweeps missing from end of last repeat
        setattr(stim, 'trial', tr)

    # hide blank stimuli
    blank_mask = stim.v==65535 # store indices of blank stimuli in mask
    stim.v[blank_mask] = 0 # temporary, so mapping below works

    # map stim.v to unique stimulus condition
    for i, dname in enumerate(dnames):
        sweeptable = np.array(stim.i[dname]*np.multiply.reduce(n[:i])).repeat(np.multiply.reduce(n[i+1:]))
        setattr(stim, dname, sweeptable[stim.v].astype(float)) # add stimulus dimension attribute
        stim.__getattribute__(dname)[blank_mask] = np.nan # blank stimulus frames not a number
    stim.v[blank_mask] = 65535 # restore blank v values

    return stim


def loadMultiSpikeTimes(varargins):
    """
    Batch load all .spk files from specified or current directory.
    Uses loadSpikeTimes.

    Optional varargins: 'dirname', 'filter'.
    """
    pass

def loadLFPs(varargins):
    """
    Batch load all LFP time-series from specified or current directory.
    Parses LFP channel numbers from filenames. Inherited from loadLFP.

    Optional varargins: 'dirname', 'filter'.
    """
    pass

def loadStim(filename,offset=0,units=1e-6,dtype='uint64',byteorder='little',vtype='i'):
    """
    Load time-stamps and values for stimuli (eg. movie frame
    onset/index) from file and convert to seconds. Takes input string
    argument 'filename'. Returns an nd Event array structure with .x,
    .t and .info [optional] ???sweepranges fields. sweepranges is a
    cell array containing the start and stop times of each sweep
    (indexed into with 1 base).

    Optional varargins: 'header';
    """
    pass


def loadStims(varargins):
    """
    stims = loadstims(varargins);

    Loads all .din files in the specified or current directory as
    stimulus structures.  Also loads all sweeptable files in the
    specified or current directory with the same base file name, and
    stores them in a sweeptable field in the appropriate  stimulus
    structure. Stimulus indexes are zero based.  Uses loadStim.

    Optional varargins: 'dirname', 'filter', 'header';
    """
    pass


def loadComments(filename,varargin):
    """
    Load from file annotations describing recording, or computer
    generated messages (eg. warnings) pertaining to recording.
    """
    pass


def saveTimeSeries(obj,filename,dtype='uint16',byteorder='little'):
    """
    Save time series data (e.g. LFP) to file.
    """
    pass


def saveEventTimings(obj,filename,units=1e-6,dtype='uint64',byteorder='little'):
    """
    Save time stamps for events (e.g. spikes) to file.
    """
    pass


def saveEvents(obj,filename,units=1e-6,dtype='uint64',byteorder='little'):
    """
    Save time stamps and values for events (e.g. frame onsets) to file.
    """
    pass

def cutSweeps(sweeptable,varargin):
    """
    cSweeps = cutSweeps(sweeptable,varargin);

    Given a sweeptable (a cell array of structures, see
    loadsweeptable) and stimulus parameter-value pairs, returns a list
    of sweep indices for which the specified parameters have the
    specified values.

    eg.cutsweeps(sweeptable,'ori',90,'speedDegSec',6) returns all the
    sweep indices (ordered from lowest to highest) with these
    parameter values. If 'sweeptable' is the only arg, returns all the
    (0-based) sweep indices in the sweeptable

    """
    pass


def loadSweepTable(filename):
    """
    sweeptable = loadSweepTable(filename);

    Given a sweep-table filename, returns a dictionary containing
    the stimulus sweep table. eg. sweeptable[7].ori contains the
    orientation of the stimulus on the 8th sweep. The attribute names
    come from those listed in the first line of the file.

    """
    pass


def parseStim(filename,varargin):
    """
    stim=parseStim(filename,varargin)

    Splits a single stimulus DIN file into multiple stimulus DIN files
    input file should have 64 bit timestamps / 16 bit word pairs
    output files will have 64 bit timestamps / 16 bit word pairs

    Make redundant by always exporting separate stimulus files?
    """
    pass


def BuildSweepRanges(dinEvents,varargins):
    """
    sweepRanges = buildsweepranges(din,varargins)

    Returns a cell array 'sweepranges' with 1-based sweep indices. The
    contents of each entry in 'sweepranges' is an nx2 array of the
    start and stop times of that sweep index, one row per occurrence.
    This is meant to be run and the results stored before they're
    needed (like when running loadstim()), as a field in a stimulus
    structure. It can then be accessed much more quickly than
    searching through the entire din (as you do with cutdin and
    cutdins) every time you want to get the start and stop times of
    some sweep index.

    Optional argument 'mode' determines what timestamp to consider as
    the stop time. In 'short' mode, the stop times are the timestamps
    of the last consecutive occurence of 'sweepi'. 'normal' mode is
    the same as 'short' except that an extra 'refreshtime' is added to
    each stop time. Default 'refreshtime' is the difference between
    the first 2 timestamps in din. In 'long' mode, the stop times
    returned are the timestamps of the sweep immediately following the
    last consecutive occurence of 'sweepi'. If there is a postsweep
    delay, 'long' mode will include it in the returned stop times.

    Optional varargins: 'mode', 'refreshtime'

    """
    pass


def regroupSweepIdxs(oldstims, sweepgrouping):
    """
    newstims = regroupSweepIdxs(oldstims, sweepgrouping);

    Regroups the sweepis for stims according to sweepgrouping, which
    is an ngroupings long cell array,  each entry being an array of
    length 2 containing the range of sweepis to replace. sweepis are
    replaced by the 0-based index of the range in the cell array

    eg: sweepgrouping = {[0 249], [250 500]}; replaces all
    sweepindices from 0 to 249 with a 0, and all those from 250 to 500
    with a 1, and leaves the rest unchanged.
    """
    pass


def getParamVals(sweeptable,paramname):
    """
    vals = getParamVals(sweeptable,paramname);

    Returns all unique paramater values for the specified parameter in
    the specified sweeptable.
    """
    pass


def crop(image,roi):
    """
    Crop spatial extend of image to roi=[xmin, xmax, ymin, ymax]
    """
    [xmin, xmax, ymin, ymax] = roi
    # xmin = roi['left']
    # xmax = roi['right']
    # ymin = roi['top']
    # ymax = roi['bottom']
    return image[..., ymin:ymax, xmin:xmax]




