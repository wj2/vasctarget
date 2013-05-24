"""
This module contains the base classes for different data types
"""

import numpy as np
import bisect

import utils
import datalib as dat
import dsplib as dsp


__all__ = ['TimeSeries','Events','Epochs','Info','Image','Container']

class Info(dict):
    """
    Base class providing data-specific meta information as attributes.

    :Parameters:
        filename : string
          optional filename of info-file to parse
    """
    def __init__(self,i=None,filename=None,translation=None,**args):
        if i != None:
            self.update(i)
        self.update(args)
        if filename != None:
            self.load(filename,translation)
    def __getattr__(self,key):
        return self[key]
    def __setattr__(self,key,value):
        self[key] = value
    def copy(self):
        return Info(i=dict.copy(self))
    def translate(self,translation):
        ''' translate info key names '''
        for k in translation.keys():
            if k in self.keys():
                self[translation[k]] = self.pop(k)
    def load(self,filename,translation=None):
        # note: makes all keys lower case
        import re
        import string
        pcomment = re.compile('^#')
        pempty = re.compile('^\s*$')
        pvariable = re.compile('^([\w\.]+)\s*=\s*(\S+[^\r\n]*)\s*$')
        f = open(filename)
        ln = 0
        for line in f.readlines():
            ln += 1
            if not pcomment.match(line) == None: continue
            if not pempty.match(line) == None: continue
            match = pvariable.match(line)
            if match == None:
                assert 0, 'syntax error in info file: %s\n[%d] %s'%(filename,ln,line)
            # varname = string.lower(re.sub('^[^\.]+\.','',match.group(1))) # use lower case keys
            varname = re.sub('^[^\.]+\.','',match.group(1))
            value = match.group(2)
            # print '%s = %s'%(varname,value)
            try:
                self[varname] = eval(value)
            except:
                self[varname] = str(value)
        if translation != None:
            self.translate(translation)

class TimeSeries(np.ndarray):
    """
    Base class for continuously sampled data.
    TimeSeries is either 1d (no trials) or 2d (several trials).
    """

    def __new__(subtype, data, samplingrate=1, i=None,
                dtype=None, copy=False, **more_info):

        # print "__new__ received %s" % type(data)
        # Make sure we are working with an array
        if hasattr(data, 'next'):
            # data comes from an iterable
            array = np.fromiter(data, dtype=dtype)
        else:
            array = np.array(data, dtype=dtype, copy=copy)

        # transform 'array' from an ndarray to our new subclass.
        array = array.view(subtype)

        # initialize info attribute
        array.i = Info()
        array.i.t0 = 0
        if i != None:
            array.i.update(i)

        # update info parameters
        array.i.samplingrate = samplingrate
        array.i.update(more_info)
        return array

    def __array_finalize__(self, obj):
        if not hasattr(self, 'i'):
            self.i = getattr(obj, 'i', None)

    def __getitem__(self, key):
        if type(key) == np.float:
            key = Events(key)
        if isinstance(key, Epochs):
            return dat.cut_time_series(self,key)
        elif isinstance(key, Events):
            return self[key.asidx(self.i.samplingrate)].view(np.ndarray)
        else:
            i = self.i # save info
            obj = self.view(np.ndarray)[key].view(type(self))
            if not np.isscalar(obj):
                obj.i = i
            return obj

    def has_trials(self):
        "returns the number of trials"
        return dat.has_trials(self)

    def len(self):
        return self.shape[-1]

    def time(self):
        "returns the time axis corresponding to the TimeSeries in s"
        dt = 1./self.i.samplingrate
        return Events(np.linspace(self.i.t0,self.i.t0+self.len()*dt,self.len(),endpoint=False))

    def duration(self):
        "returns the duration of the TimeSeries in s"
        dt = 1./self.i.samplingrate
        return self.len()*dt

    def tstop(self):
        "returns the time at the end of the TimeSeries in s (last timestamp + binwidth)"
        return self.i.t0+self.duration()

    def concatenate(self,data):
        obj = np.concatenate((self,data)).view(type(self))
        obj.i = self.i.copy()
        return obj

    def smooth(self,width=1,**kargs):
        """
        smoothes timeseries using a gaussian window with standard deviation width (in s).

        other parameters:
        ns : number of standard deviations used for window (on each side of the window center)
        """
        args = dict(width=width)
        args.update(kargs)
        np.ndarray.__setitem__(self,slice(None),dsp.smooth_timeseries(self,**args))
        return self

    def mean(self,*args,**kargs):
        return TimeSeries(self.view(np.ndarray).mean(*args,**kargs),samplingrate=self.i.samplingrate,i=self.i)

    def iirfilter(self,N_ord=3,Wn=[.5, 150.0],overwrite=True):
        """
        filters timeseries using an IIR filter of order N_ord and corner frequencies Wn (in Hz).

        Modifies the timeseries in place unless overwrite=False (alternatively use dsplib.iirfilter)
        """
        if overwrite:
            dat = self
        else:
            dat = self.copy()
            dat.i = self.i.copy() # XXX - where else? deep copy or explicit?
        np.ndarray.__setitem__(dat,slice(None),dsp.iirfilter(dat,N_ord, Wn, self.i.samplingrate))
        return dat
        
    def decimate(self,q):
        """
        Decimate signal x by an integer factor q (without interpolation)
        Updates info accordingly.

        :Parameters:
            signal -- the timeseries to be decimated
            q -- the decimation factor
        """
        if type(q) != int and q > 1: raise Error, "q should be a positive integer > 1"
        self = self[::q]
        self.i.samplingrate/=q
        self.i.nsamples/=q
        return self


class Events(np.ndarray):
    """
    Base class for discrete event data (time stamps and optional atttibutes)
    """
    def __new__(subtype, data, i=None, dtype=None, copy=False, **valuedict):
        # print "__new__ received %s" % type(data)
        if (type(data) == list) and (len(data) > 0):
            if not np.isscalar(data[0]):
                # print "received multiple trials", data
                return dat.events_from_trials(data,i=i,dtype=dtype,copy=copy,**valuedict)
        subarr = np.array(data, dtype=dtype, copy=copy).view(subtype)
        assert subarr.ndim <= 1, "Events have to be at most 1d!"

        # copy keys
        if hasattr(data, 'keys'):                   # Use data's 'values' attribute if it exists
            for key in data.keys():                 # update values from data
                setattr(subarr, key, getattr(data, key))

        # copy values
        for key,value in valuedict.items():         # Use the specified 'values' argument if given
            setattr(subarr,  key, np.array(value))   # update keys from data

        # copy info
        subarr.i = Info()                               # initialize info
        if hasattr(data, 'i'): subarr.i.update(data.i)  # update info from data
        if i != None:                                   # update info from arguments
            subarr.i.update(i)

        return subarr

    def __array_finalize__(self,obj):
        if not hasattr(self, 'i'): self.i = getattr(obj, 'i', None)
        if hasattr(obj,'keys'):
            for key in obj.keys():
                setattr(self,  key, getattr(obj, key))

    def valuedict(self):
        valuedict = self.__dict__.copy()
        valuedict.pop('i')
        return valuedict

    def keys(self):
        return self.valuedict().keys()

    def __getitem__(self, key):
        # print "__getitem__ received %s, type %s" % (str(key),type(key))
        if isinstance(key, Events):
            if key.dtype == bool:
                key = key.view(np.ndarray)
        if isinstance(key, Epochs):
            # The key is an epoch
            return dat.cut_events(self,key).view(type(self))
        elif (type(key) == str) and key in self.keys():
            return getattr(self,key)
        else:
            subarr = self.view(np.ndarray)[key]
            if np.isscalar(subarr):
                return subarr
            else:
                i = self.i # save info
                obj = subarr.view(type(self))
                obj.i = i
                for k in self.keys():
                    setattr(obj,k,getattr(self,k)[key])
                return obj

    def __getslice__(self,i,j):
        # print "__getslice__ received %s,%s" % (str(i),str(j))
        obj = np.array(np.ndarray.__getslice__(self,i,j)).view(type(self))
        obj.i = self.i
        for k in self.keys():
            setattr(obj,k,np.ndarray.__getslice__(getattr(self,k),i,j))
        return obj

    def __iter__(self,*args,**kargs):
        if self.has_trials():
            return self.iter_dim(*args,**kargs)
        else:
            return self.view(np.ndarray).__iter__(*args,**kargs)

    def concatenate(self,other):
        """
        Concatenates events with other events
        """
        obj = np.concatenate((self,other)).view(type(self))
        obj.i = self.i.copy()
        for k in self.keys():
            setattr(obj,k,np.concatenate((getattr(self,k),getattr(other,k))))
        return obj

    def has_trials(self):
        "returns the number of trials"
        return dat.has_trials(self)

    def has_events(self,dim_name='trials'):
        "returns the number of trials with one or more events"
        if hasattr(self, dim_name):
            return len(set(getattr(self,dim_name)))
        else:
            return 0

    def events_per_trial(self):
        "returns the number of events in each trial"
        if hasattr(self, 'trials'):
            return [(self.trials == t).sum() for t in range(self.has_trials())]
        else:
            return 0

    def select_dim(self,id_dim,dim_name='trials'):
        "returns the events of a specific trial"
        out = self[getattr(self,dim_name)==id_dim].copy()
        del out.__dict__[dim_name]
        out.i = self.i.copy()
        if dim_name == 'trials':
            out.i.pop('ntrials')
        return out

    def iter_dim(self,dim_name='trials'):
        for id_dim in set(getattr(self,dim_name)):
            yield self.select_dim(id_dim,dim_name=dim_name)

    def asidx(self,samplingrate):
        """
        returns the event timestamps as integer indices

        Allows time-stamps to be used as an indices of a TimeSeries
        with given samplingrate.
        """
        if self.has_trials():
            return [self.trials.tolist(),((self+utils.eps)*samplingrate).astype('i').tolist()]
        else:
            return ((self+utils.eps)*samplingrate).astype('i').tolist()


    def asraster(self,samplingrate,**kargs):
        """
        Returns event raster (binary timeseries) sampled with samplingrate

        optional arguments:

        nbins: if given, it determines the length of the raster.
        t0:  if t0<0 (t0>0) is given, sample points are added (removed) at the beginning of the raster in order
             to make the first sample point to coincide with t0.
        """
        args = dict(binary=True)
        args.update(kargs)
        return dat.events_to_time_series(self,samplingrate,**args)

    def asrate(self,samplingrate,**kargs):
        """
        Returns event rate (in 1/s) sampled with samplingrate, assuming given event width (in s)

        optional arguments:

        width: stdev of smoothing window
        nbins: if given, it determines the length of the raster.
        t0:  if t0<0 (t0>0) is given, sample points are added (removed) at the beginning of the raster in order
             to make the first sample point to coincide with t0.
        """
        args = dict(binary=False)
        args.update(kargs)
        return dat.events_to_time_series(self,samplingrate,**args)*samplingrate

#     def time2idx(self,time):
#         """
#         Returns indices corresponding to given time-stamps.
# 
#         If the exact event doesn't exist for a given time, the closest _previous_
#         event index is returned. Time-stamps that preceed the first event return -1, those
#         later than the last event return the last index of the events array.
#         """
#         return np.array(map(lambda x: bisect.bisect(self.view(np.ndarray),x)-1,time),'i')

    def evt2idx(self,events):
        """
        Returns corresponding time-stamp indices of another event array.

        If the exact time-stamp in 'events' doesn't exist for a given element,
        the closest _previous_ time-stamp index is returned. Time-stamps that
        preceed the first time-stamp of the 'events' array return -1, those
        later than the last last-stamp return the last index of the 'events' array.
        """
        return np.array(map(lambda x: bisect.bisect(events.view(np.ndarray),x)-1,self),'i')

    def epoch2slice(self,epoch):
        """
        Returns a slice that corresponds to given epoch

        Uses bisect method and assumes that events are sorted.
        """
        if type(epoch) == tuple:
            tstart,tstop = epoch
        else:
            tstart,tstop = (epoch['tstart'],epoch['tstop'])
        imin = bisect.bisect_left(self.view(np.ndarray),tstart)
        imax = bisect.bisect_right(self.view(np.ndarray),tstop)
        return slice(imin,imax)

    def sortedslice(self,epoch):
        """
        Slices events corresponding to given epoch

        Uses bisect method and assumes that events are sorted.
        """
        if type(epoch) == tuple:
            tstart = epoch[0]
            toffset = 0
        else:
            toffset = epoch['tstart']+epoch['t0_offset']
        return self[self.epoch2slice(epoch)] - toffset


    def sort(self):
        """
        Sorts timestamps and rearranges all other items accordingly
        """
        idx = self.view(np.ndarray).argsort()
        self.view(np.ndarray).sort()
        for k in self.keys():
            setattr(self,k,np.ndarray.__getitem__(getattr(self,k),idx))

    def shift(self,shifts):
        """
        shifts event timings by shifts.
        if evts have trials, shift has to be iterable.
        """
        ntrials = self.has_trials()
        if ntrials==None:
            assert np.isscalar(shifts), 'shifts has to be a scalar for single trials'
            self[:] += shifts
        else:
            assert len(shifts) == ntrials, 'shifts has to have length %d'%ntrials
            for i in range(ntrials):
                self[self.trials==i] += shifts[i]

    def __repr__(self):
        nr = self.has_trials()
        if nr:
            tr = "\n %(trials)s trials" % {'trials': str(nr)}
            return repr(self.select_dim(0))+tr

        s = repr(self.__array__()).replace('array', 'Events')
        # now, 'Events' has 6 letters, and 'array' 5, so the columns don't
        # line up anymore. We need to add a space.
        l = s.splitlines()
        for i in range(1, len(l)):
            if l[i]:
                l[i] = ' ' + l[i]
        s = '\n'.join(l)

        for key in self.keys():
            sv = repr(self.__dict__[key].__array__()).replace('array', key)
            l = sv.splitlines()
            for i in range(1, len(l)):
                if l[i]:
                    l[i] = ' ' + l[i]
            sv = '\n'+'\n'.join(l)
            s = s+sv
        return s


_epochtype = np.dtype({'names':['tstart','duration','t0_offset','tstop'],'formats':['d']*4})
class Epochs(np.ndarray):
    """
    Base class defines time intervals.

    Epochs can be single time intervals, or lists of time intervals.

    Time intervals are given by either
    - pairs of time points denoting tstart and tstop
    - triples of a tstart, duration, t0_offset

    """
    def __new__(subtype, data, data2=None, name=None, align=True):
        if data2 != None:
            if type(data2) == str:
                name = data2
            else: data = zip(data,data2) # two lists: tstart and tstop
        # print "__new__ received %s" % type(data)
        if type(data) == tuple:
            if len(data) == 2: # list of t_start, t_stop
                if align:
                    t0_offset = 0
                else:
                    t0_offset = -data[0]
                subarr = np.array((data[0],data[1]-data[0],t0_offset,data[1]), _epochtype)
            elif len(data) == 3: # list of t_start, duration, offset
                subarr = np.array(data + (data[0]+data[1],), _epochtype)
            else: assert 0, 'list of 2-tuples or 3-tuples required'
        elif type(data) == list:
            if type(data[0]) != tuple: assert 0, 'list of tuples required'
            if len(data[0]) == 2: # list of t_start, t_stop
                if align:
                    toff = 0
                else:
                    toff = -1
                subarr = np.array([(e[0],e[1]-e[0],toff*e[0],e[1])for e in data], _epochtype)
            elif len(data[0]) == 3: # list of t_start, duration, offset
                subarr = np.array([e + (e[0]+e[1],) for e in data], _epochtype)
            else: assert 0, 'list of 2-tuples or 3-tuples required'
        else: assert 0, "can't typecast data of %s" %type(data)

        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)
        if name != None:
            subarr.name = name
        return subarr

    def __array_finalize__(self, obj):
        if not hasattr(self, 'name'):
            name = getattr(obj, 'name', None)
            if name != None:
                self.name = name

    def __getitem__(self, key):
        # print "__getitem__ received %s, type %s" % (str(key),type(key))
        subarr = np.ndarray.__getitem__(self,key)
        if (isinstance(key,str) and ~np.isscalar(subarr)):
            # extract single column (tstart, etc), don't return Epochs
            return Events(subarr) # return Events
        elif (not isinstance(key,str) and np.isscalar(subarr)): # single epoch, return Epochs class
            obj = np.array(subarr,_epochtype).view(type(self))
            if hasattr(self,'name'): obj.name = self.name
            return obj
        else:
            return subarr

    def asslice(self,samplingrate,nmin=0,nmax=np.inf,t0=0,include_last_bin=False):
        # correct for time offset of timeseries data
        dn = int(np.ceil(t0*samplingrate))
        if self.ndim == 0:
            if include_last_bin:
                idx = (int((self['tstart'])*samplingrate),int(np.ceil((self['duration']-utils.eps)*samplingrate))) #incl. last bin
            else:
                idx = (int((self['tstart']+utils.eps)*samplingrate),int((self['duration']+utils.eps)*samplingrate)) #excl. last bin
            return slice(idx[0]-dn,idx[0]+idx[1]-dn)
        else:
            if include_last_bin:
                idx = [(int((e['tstart'])*samplingrate),
                        int(np.ceil((e['duration']-utils.eps)*samplingrate))) for e in self]  #incl. last bin
            else:
                idx = [(int((e['tstart']+utils.eps)*samplingrate),
                        int((e['duration']+utils.eps)*samplingrate)) for e in self]  #excl. last bin
            slices = [slice(i[0]-dn,i[0]+i[1]-dn) for i in idx
                      if not ((i[0]-dn<nmin) | (i[0]+i[1]-dn>nmax))]
            return np.array(slices)

    def __repr__(self):
        name = getattr(self,'name','Epochs')
        s = repr(self.__array__()).replace('array', name)
        return s

class Stimulus(Events):
    """
    Base class for visual stimulus information
    """
    def trials_epochs(self,vstart=0,vstop=0):
        "returns epochs corresponding to repeated presentations"
        # keep slightly modified version of trials_epochs and
        # make a new function "sweeps_epochs?" to deal with sweep-based stimuli?
        # some of the code in here should be in custom_tim_data
        # esp. stuff that uses DimStim specific attributes and fields
        # scalar or vector v specifies trial onset DIN index, default=0

        # check for crashes if vstart or vstop == 0 == false for all DIN

        frms_per_sweep = round(self.i.refreshrate/self.i.framerate)
        onsets = self[self.v==vstart][::frms_per_sweep]
        offsets = self[self.v==vstop][frms_per_sweep-1::frms_per_sweep]
        duration = min(offsets-onsets)+1/self.i.refreshrate

        # get duration of sweep
        if vstart==vstop==0:
            if len(onsets) > 1:
                duration = np.diff(onsets).min()
            else: # fix for non-repeating stimuli
                duration = self[-1]+self[1]-2*self[0]

        # get stimulus epochs
        epochs = dat.epochs_from_onsets(onsets,duration)
        return epochs

    def __repr__(self):
        nr = self.has_trials()
        if nr:
            tr = "\n %(trials)s trials" % {'trials': str(nr)}
            return repr(self[0])+tr

        s = repr(self.__array__()).replace('array', 'Stimulus')
        # now, 'Stimulus' has 8 letters, and 'array' 5, so the columns don't
        # line up anymore. We need to add a space.
        l = s.splitlines()
        for i in range(1, len(l)):
            if l[i]:
                l[i] = '   ' + l[i]
        s = '\n'.join(l)

        for key in self.keys():
            sv = repr(self.__dict__[key].__array__()).replace('array', key)
            l = sv.splitlines()
            for i in range(1, len(l)):
                if l[i]:
                    l[i] = ' ' + l[i]
            sv = '\n'+'\n'.join(l)
            s = s+sv
        return s


class Image(np.ndarray):
    """
    Base class for static images (2d) or movies (3d).
    """
    def __new__(subtype, data, i=None, dtype=None, copy=False):
        # transform 'data' from an ndarray to the new subclass
        subarr = np.array(data, dtype=dtype, copy=copy).view(subtype)
        # copy info
        subarr.i = Info()                               # initialize info
        if hasattr(data, 'i'): subarr.i.update(data.i)  # update info from data
        if i != None:                                   # update info from arguments
            subarr.i.update(i)

        return subarr

    def __array_finalize__(self,obj):
        if not hasattr(self, 'i'): self.i = getattr(obj, 'i', None) # update info

    def crop(self,roi):
        """
        Crop spatial extend of image to roi=[xmin, xmax, ymin, ymax]
        """
        return dat.crop(self,roi)

class Note:
    """
    Base class for text string with time stamp and date.
    """
    def __init__(self,msg,time,date):
        self.msg = msg    # message string
        self.t = time     # precision (float) time stamp [s]
        self.date = date  # ANSI date


class Container(dict):
    """
    Container class that groups multiple Events or Timeseries objects.

    A container groups multiple objects, such as spike timings
    corresponding to multiple neurons.
    Individual objects can be addressed using a key, such as the neuron ID.
    """

    def __getitem__(self,key):
        if type(key) in [list,tuple]:
            return Container([(k,self[k]) for k in key])
        else:
            return dict.__getitem__(self,key)

    def batch_process(self,function,overwrite=False,*params,**kparams):
        obj = Container()
        for k in self.iterkeys():
        # for k in np.sort(self.keys()):
            obj[k] = function(self[k],*params,**kparams)
            if overwrite:
                self[k] = obj[k]
        return obj


class Recording:
    """
    Container class that groups data sets related to a given
    recording session (e.g. spike timings, LFP, stimulus data).

    Inputs:
      recording_id - a tuple that specifies the recording
      info - a dictionary with meta information
    """
    def __init__(self):
        """
        Initialization method

        this method has to be customized
        """
        pass

    def load_cache(self,arrayname,samplingrate,dtype):
        """
        load pre-computed data from cache file
        """
        outfile = utils.load_cache(arrayname,self.i.cache_path,dtype)
        if outfile != None:
            return TimeSeries(outfile,samplingrate=samplingrate,dtype=dtype)
        else:
            return None

    def save_cache(self,data,arrayname):
        """
        save computed array data to cache file
        """
        utils.save_cache(data,arrayname,self.i.cache_path)

    def load_lfp(self,id_lfp):
        '''
        load LFP data

        Input:
          id_lfp - an id that identifies the LFP recording

        this function has to be customized
        '''
        assert 0, 'this function has to be customized'

    def load_spk(self,id_spk):
        '''
        load spike timing data

        Input:
          id_spk - an id that identifies the spike recording

        this function has to be customized
        '''
        assert 0, 'this function has to be customized'

    def load_stim(self):
        '''
        load stimulus timing data

        Input:
          id_spk - an id that identifies the spike recording

        this function has to be customized
        '''
        assert 0, 'this function has to be customized'

    def keys(self):
        ''' returns a list of attributes (except i)'''
        return [k for k in self.__dict__.keys() if k != 'i']


class Electrode:
    """
    Base class for recording electrode (or electrode array)
    """
    pass
