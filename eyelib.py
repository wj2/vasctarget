"""
Library of functions related to eye movement analysis
"""
import re
import numpy as np
from data import Epochs,Events,TimeSeries,Container
# from neuropy import dsplib
import time
from StringIO import StringIO

#try:
#    from Dimstim import pylink
#except ImportError,x:
#    print x

__all__ = [ 'Eyelink', 'detect_blinks', 'detect_saccades_trivial',
'smooth_boxcar', 'detect_saccades','reprocess_eyelink_msgs', 'findall_loadtxt',
'read_eyelink','read_eyelink_cached', 'ck_read_array', 'plot_eyetrace',
'EyelinkReplayer']

def analog_to_angles(delta_V):
    """ convert analog voltage trace to a rough approximation of degrees

    returns k , deg/V conversion factor
    """
    # 22 inch monitor, assuming 4:3 aspect ratio
    screen_height = 22. * 0.0254  * (3./5.)
    screen_height /= 2.
    #Thom Carney says viewing distance varies between .5 and .7 meteres
    viewing_distance = .5
    theta = np.arctan2(screen_height, viewing_distance)
    theta *=2.
    k = theta / delta_V
    k *= 10000 # empirical fudge factor so that voltage end up in the same ballpark as angles
    return k

def detect_blinks(p,x=None,y=None,offset=1):
    """
trivial blink detector, returns Epochs of blinks.

If you can, you should use Eyelink's own blink detector, this is a very poor
replacement.

:Parameters:
    p : TimeSeries
        pupil area time series

:Returns:
    : Epochs
        times when eyes are closed

:SeeAlso:
    - Eyelink : class whose ['discard'] field contains thorough blinks detected
    by the eye tracker (see Notes section of Eyelink's docstring)

Notes
-----

This code relies on p==0 when eyes are closed. 
    """
    # extract times from the last non-zero to the last zero
    blink_thresh = -1000
    #blinks = p.time()[np.logical_xor(
    #    (np.diff(p)<blink_thresh), (np.diff(p)>-blink_thresh))]
        #(np.diff(p)=p[1:]), (np.diff(p)==-1*p[:-1]))]
    boxcar = 301.*x.i.samplingrate/ 1000. #31 ms boxcar filter
    ps = smooth_boxcar(p,boxcar)
    boxcar = 101.*x.i.samplingrate/ 1000. #31 ms boxcar filter
    pss = smooth_boxcar(p,boxcar)
    #eye_closed = p < ps.mean() - ps.std()*2
    eye_closed = p < ps - ps.std()*3
    #eye_closed = np.logical_or(eye_closed, p < ps - ps.std()*1.5)
    eye_closed = np.logical_or(eye_closed, p < pss - pss.std()*3)
    eye_closed = np.logical_or(eye_closed, p < .95*p.min())
    blinks = p.time()[np.logical_xor(eye_closed[:-1],eye_closed[1:])]
    if len(blinks):
        if len(blinks) % 2:
            blinks = blinks.concatenate([p.time()[-1]])
        #print Epochs(blinks[0+offset::2],blinks[1+offset::2])
        return Epochs(blinks[0+offset::2],blinks[1+offset::2])
    else:
        return None

    blinks = p.time()[(np.diff(p)<blink_thresh)]
    blinke = p.time()[(np.diff(p)>-blink_thresh)]
    #if len(blinks) != len(blinke)
    #    for i in range(len(blinks):

    #use blinks[::2] and blinks[1::2] to get the onsets and ends of blinks
    if len(blinks):
        if len(blinks) != len(blinke):
            blinke = blinke.concatenate([p.time()[-1]])
        print Epochs(blinks,blinke)
        return Epochs(blinks,blinke)
        #print Epochs(blinks[::2],blinks[1::2])
        #return Epochs(blinks[::2],blinks[1::2])


def detect_saccades_trivial(e,x=None,y=None):
    """
    trivial saccade detector, returns Events of when saccades occured
    """
    # use diff 
    dx = np.diff(x) #dx = x[1:] - x[:-1]
    dy = np.diff(y) #dy = y[1:] - y[:-1]
    dx_s = smooth_boxcar(dx,131)
    dy_s = smooth_boxcar(dx,131)
    saccades = Events(dx.time()[abs(dx_s)>100000])
    print len(saccades), "sacs ", 
    # create a bitmask of when dx or dy cross some threshold?
    # make sure it's not a blink
    return saccades

def smooth_boxcar(x, length):
    """ 
    Simple boxcar smoothing (as I understand it -Paul)
    """
    boxcar = np.ones(length) #31 ms boxcar filter
    boxcar /= boxcar.size
    return np.convolve(x,boxcar,'same') #[boxcar.size/2-1:-boxcar.size/2]

def detect_saccades(e,x=None,y=None):
    """
Saccade detector, returns Events of when saccades occurred

Ref: Martinez-Conde, Macknik, Troncoso, and Dyar. Microsaccades counteract visual fading during fixation. Neuron 49, 297-305 (2006).

:Parameters:
    e : Eyelink or TimeSeries
       either all data from an eyelink .asc file or just the pupil area time
       series if using all three parameters
    x : None or TimeSeries
        x-axis trace if :e: is a TimeSeries
    y : None or TimeSeries
        y-axis trace if :e: is a TimeSeries

:Returns:
    : Events
        begininngs and ends of saccades

:SeeAlso:
    - Eyelink : class whose ['saccades'] field contains saccades detected by the
    eye tracker 

Notes
-----

This code is not finished yet, but seems to work well if edf2asc ran with the
option to output HREF angles (-sh) instead of default gaze position
    """
    if x==None and y==None:
        x=e['x'];y=e['y'];p=e['pupA']
    # use diff 
    dx = np.diff(x) #dx = x[1:] - x[:-1]
    dy = np.diff(y) #dy = y[1:] - y[:-1]
    boxcar = 31.*x.i.samplingrate/ 1000. #31 ms boxcar filter
    dx_s = (dsplib.smooth_timeseries(dx,.010))
    dy_s = (dsplib.smooth_timeseries(dy,.010))
    r = np.sqrt( dx_s**2 + dy_s**2) # magnitude of motion XXX: convert to degrees from pixels

    time = x.time()
    # discard velocities less than 3 deg/second
    theta = np.arctan2(dy_s,dx_s) *180.0/np.pi
    theta[theta<0] += 360.0 # theta < 0 same as dx_s>0 & dy_s<0 after previous line

    # 'rate-of-turn' indicator, check that movement changed less than 15 degrees
    #v_thresh is 3 deg/sec, but we're not in degrees yet, here
    v_thresh= 3000.0/x.i.samplingrate  # 3000 is 3 degrees, 1 seconds is 1000 Hz
    abs_dtheta= np.abs(np.diff(theta))
    #stopped = np.logical_or(
    stopped = r[1:] < v_thresh
              #  abs_dtheta > 10.0 )
    stopped = np.logical_xor(stopped[:-1],stopped[1:])

    saccades = time[stopped]
    # ok, now we have all saccade endings, extract the saccades themselves

    # XXX: this is totally ad hoc, just here to make plot_eeg example somewhat
    # meaningful create a bitmask of when dx or dy cross some threshold? then
    # calculate the combined distance traveled?  make sure it's not a blink
    
    return Events(saccades)

class AttribContainer(Container):
    """
    class where all attributes are stored in the dictionary.

    Used so that we can do things like:
        e.r.saccades.vpeak 
    instead of 
        e['r']['saccades'].vpeak

    XXX: incorporate this functionality into the base Container class. There's
    probably a more appropriate way of doing this
    """
    def __init__(self):
        Container.__init__(self)

    def __getattr__(self,key):
        return self[key]
    def __setattr__(self,key,value):
        self[key] = value

def __eye_repr__(eye):
    rep = "Eye data:"
    rep += "\n time: %.3f-%.3f (%.3f s)"% ( eye['x'][0],
            eye['x'][-1], eye['x'][-1]-eye['x'][0])
    for x in ['x','y','pupA']:
        rep += "\n %s data: shape=%s " %(x,str(eye[x].shape))
        rep += "%s0=%.1f" % (x, eye[x].v[0])
        rep += " %sf=%.1f" % (x, eye[x].v[-1])
    for x in ['fixepochs','sacepochs']:
        if len(eye[x]):
            rep +='\n %d %s' %(len(eye[x]),x)
            d = eye[x]['duration']
            rep +=", mean duration = %.3f +/- %.4f"%(d.mean(),d.std())
    rep+= "\n other data attached: discard epochs, fixations and saccades"
    rep+= " as parsed by Eyelink"
    return rep

class Eyelink(AttribContainer):
    """
class for Eyelink data (extends AttribContainer)

Has x,y,pupA time series, as well as saccade and discard epochs, frames
events and msgs strings extracted from an Eyelink .asc file. All members are
also attributes (you can use e.x and e.saccades as a shorthand for e['x'] and
e.['saccades']

:Members:
    e["x"] : TimeSeries
        The x-axis trace
    e["y"] : TimeSeries
        The y-axis trace
    e["pupA"] : TimeSeries
        Pupil area trace
            
    e["saccades"] : Epochs or None
        Epochs which contain saccades 
    e["blinks_l"] : Epochs or None
    e["blinks_r"] : Epochs or None
        Blinks from the left and right eyes
    e["discard"] : Epochs or None
        Epochs to be discarded because they contain a blink (see Notes below)

            
    e["frames"] : Events or None
        Onset of frames, e.g. e["frames"][29] returns when 29th frame was shown 
    e["msgs"] : tuple(str) or None
        Lines in Eyelink file which start with MSG (such as "Stimulus paused"),
        other than frame onsets (which are parsed out)
    e["raw"] : str
        Contains the text which was not processed into one of the above

:SeeAlso:
    - read_eyelink : function which returns an Eyelink object

Notes
-----

According to Eyelink manual 1.3 4.5.3.5 Blinks (p. 98) 

    "Blinks are always preceded and followed by partial occlusion of the and
    pupil, causing artificial changes in pupil position.  These are sensed by
    the EyeLink 1000 parser, and marked as saccades. The sequence of events
    produced is always:
    - start saccade (SSACC)
    - start blink   (SBLINK)
    - end blink     (EBLINK)
    - end saccade   (ESACC) 
    ... All data between SSACC and ESSAC events should be discarded. "

In the Eyelink 1.4 manual, p 108:

    Blinks are always embedded in saccades, caused by artificial motion as the
    eyelids progressively occlude the pupil of the eye. Such artifacts are best
    eliminated by labeling and SSACC...ESACC pair with one or more SBLINK
    events between them as a blink, not a saccade.

    It is also useful to eliminate any short (less than 120 millisecond
    duration) fixations that precede or follow a blink. These may be artificial
    or be corrupted by the blink.

XXX: We do not currently eliminate short (<120ms) fixations that precede or
follow a blink

    """
    def __init__(self,fname, binocular, have_right, have_left):
        Container.__init__(self)
        self.binocular = binocular
        self.have_right = have_right
        self.have_left = have_left
        self.r = AttribContainer()
        self.l = AttribContainer()
        self.r.saccades = AttribContainer()
        self.l.saccades = AttribContainer()
        self.r.fixepochs = None
        self.l.fixepochs = None
        self.r.sacepochs = None
        self.l.sacepochs = None
        self.r.discard = [] 
        self.l.discard = []
        #self["discard"] = self["saccades"] = None
        self["msgs"] = self["frames"] = self["raw"] = None
        self._fnamelong = fname
        self._fname = fname.split('/')[-1].split('.')[0]

    def __repr__(self):
        rep  =  "Eyelink Data:"
        if self.binocular: rep += " binocular"
        elif self.have_right: rep += " monocular (right eye)"
        else: rep += " monocular (left eye)"
        for x in ['frames','msgs']:
            if self[x] is not None:
                rep +="\n %d %s ["%(len(self[x]),x)
                rep += self[x][0].__repr__()
                rep += " ...]"
        raw = self['raw'].split('\n')

        if self.have_left:
            rep += "\nLeft " +  __eye_repr__(self.l)
        if self.have_right:
            rep += "\nRight "+ __eye_repr__(self.r)
        rep += "\n"+raw[0]
        rep += "\n"+raw[1]
        rep += "\n"+raw[5]
        return rep

    def replayer(self):
        return EyelinkReplayer(el=self)

    @property
    def eye_used(self):
        """ Return self.l or self.r, depending on which eye has data, raises an
        error if binocular"""

        if self.binocular:
            raise Exception, "Container is binocular, both eyes are available"
        if self.have_right:
            return self.r
        else:
            return self.l

    @property
    def experiment(self):
        """ the name of the python file that ran at the beginning of this
        recording.
        
        Just grabs it from the .msg[0] """
        return self.msgs[0].split('/')[-1]

    @property
    def fname(self):
        """short filename for this recording"""
        return self._fname

    @property
    def fnamelong(self):
        """Long (full) filename of this recording"""
        return self._fnamelong

def reprocess_eyelink_msgs(pattern, msgs, cols=(0,), dtype=None):
    """ 
    Takes the messages, creates a temporary buffer out of them, and reuses
    the machinery of findall_loadtxt to search for the pattern

    The pattern will be prepended with the `(?<=MSG.)\d+ ` regular expression,
    and appended with `.*` by default.

    If you don't supply a cols argument, you'll just get the message timestamp
    (column 0) as a float.

    Example:
    --------
        Grab just the timestamp of messages which start with pupArt:
        >>> reprocess_eyelink_msgs('pupArt', el.msgs)
     
    XXX: there's probably a better way of doing this than writing to a
    temporary buffer and reprocessing it, but this allows me to reuse machinery
    that's already there. and LTS. -pi
    """
    pat = "(?<=MSG.)\d+ "+pattern+".*"
    return findall_loadtxt(pat, "\n".join(msgs), cols, dtype)


def findall_loadtxt(pattern, raw, cols, dtype=None):
    matches = re.findall(pattern,raw, re.M)
    str = "\n".join(matches)
    tmp = StringIO(str)
    if tmp.len == 0:
        return np.array([]) 
    ret = np.loadtxt(tmp,dtype=dtype,usecols=cols) 
    tmp.close()
    return ret

def read_eyelink(filename='Eyelink/tim_imas.asc'):
    """ 
Read in Eyelink .asc file,

Returns obect containing eye position time series, along with pupil area,
saccades and blink data as parsed by Eyelink, and any messages (such as frame
information)


:Parameters:
    filename : str
        name of .asc file to be read

:Returns: 
    eye : Eyelink
        container with 'x','y','pupA' TimeSeries, 'saccade' and 'discard'
        Epochs, and 'frames' Events.

:SeeAlso:
  - Eyelink : class of the return object of this function 

Notes
-----

Assumes edf2asc ran with '-nflags -miss -9999.9'


Examples
--------

>>> eye = read_eyelink('data/BALDI009.asc')
>>> eye['x'], eye['y'], eye['saccades'][:2]
(TimeSeries([ 398.8,  398.8,  398.8, ...,  350.2,  355.5,  361.1]),
 TimeSeries([ 301.1,  301. ,  300.9, ...,  547.5,  512.4,  478.9]),
 Epochs([(18984432.0, 74.0, 0.0, 18984506.0),
       (18984800.0, 6.0, 0.0, 18984806.0)], 
      dtype=[('tstart', '<f8'), ('duration', '<f8'), ('t0_offset', '<f8'), ('tstop', '<f8')]))
>>>  print g['raw'][:200]
** CONVERTED FROM BALDI009.EDF using edfapi 3.0 Linux Jun 18 2008 on Wed Oct 22 16:47:03 2008
** DATE: Thu Mar  2 15:38:30 2006
** TYPE: EDF_FILE BINARY EVENT SAMPLE TAGGED
** VERSION: EYELINK II 1
**
"""
    f = file(filename)
    raw = f.read()

    # filter these redundant lines early, their E prefix counterparts contain all of the data
    # XXX: not true in HINAK075.ASC as parsed by before 2006 edf2asc - SSACS followed by blinks do not
    # have  a proper ESSAC after (but instead have a single sample long ESSAC)
    for x in ["SFIX","SSACC","SBLINK"]:
        raw = re.sub("\n"+x+".*","", raw)

    # from Eyelink manual 1.4 p 106:
    # Each type of event has its own line format. These use some of the data
    #items listed below. Each line begins with a keyword (always in uppercase) and
    #items are separated by one or more tabs or spaces.  
    #DATA NOTATIONS
    #<eye>                 which eye caused event ("L" or "R")
    #<time>                timestamp in milliseconds
    #<stime>               timestamp of first sample in milliseconds
    #<etime>               timestamp of last sample in milliseconds
    #<dur>                 duration in milliseconds
    #<axp>, <ayp>          average X and Y position
    #<sxp>, <syp>          start X and Y position data
    #<exp>, <eyp>          end X and Y position data
    #<aps>                 average pupil size (area or diameter)
    #<av>, <pv>            average, peak velocity (degrees/sec)
    #<ampl>                saccadic amplitude (degrees)
    #<xr>, <yr>            X and Y resolution (position units/degree)


    # frame data - MSG <time> <frame>
    # get first field for time, get second field as value dict
    frames = findall_loadtxt("(?<=MSG.)\d+\ \d+",raw,(0,1), dtype=int).astype(float)
    link_fix = findall_loadtxt("(?<=MSG.)\d+\ Fixation",raw,(0,), dtype=int).astype(float)
    
    gaze_dtype = np.dtype([('time', 'uint64'),
                           ('x','float64'),
                           ('y','float64')])

    gcdisp = findall_loadtxt("(?<=MSG.)\d+\ GCDISP.\d+.*",raw,(0,2,3),
            dtype=gaze_dtype)

    # get MSGs which are not like frames or GCDISP
    msgsstr = re.findall("^MSG.\d+\ [^\dG].*", raw, re.M) 

    gaze_cols = (0,1,2,3) #

    gaze_dtype = np.dtype([('time', 'uint64'),
                           ('x','float64'),
                           ('y','float64'),
                           ('pupA','float64')])

    #grab binocularity from the EVENTS GAZE line
    start = re.findall("^START.\d+\ [^\d].*", raw, re.M) 
    start = start[0].upper()

    binocular = False
    have_right = start.find('RIGHT') != -1
    have_left= start.find('LEFT') != -1
    if have_right and have_left:
        binocular = True

    if binocular:
        gaze_cols = (0,1,2,3,4,5,6) #
        gaze_dtype = np.dtype([('time', 'uint64'),
                               ('x','float64'),
                               ('y','float64'),
                               ('pupA','float64'),
                               ('x2','float64'),
                               ('y2','float64'),
                               ('pupA2','float64')])

    # EL1.4 p 103: 4.9.2 Sample Line Format (see DATA NOTATIONS for <codes>)
    # -----------------------
    #Sample lines contain time, position, and pupil size data. Optionally, velocity and
    #resolution data may be included. Several possible sample line formats are
    #possible. These are listed below.
    #Essentially, each sample line begins with a timestamp. Recordings done with a
    #2000 hz sampling rate will have two consecutive rows of the same time stamps.
    #The second row refers to the sample collected at 0.5 ms after the reported time
    #stamp. The time stamp field is followed by X and Y position pairs and pupil size
    #data for the tracked eye, and optionally by X and Y velocity pairs for the eye,
    #and resolution X and Y values. Missing data values are represented by a dot
    #("."), or the text specified by the "-miss" option to EDF2ASC.
    #         SAMPLE LINE FORMATS
    # Monocular:
    #        <time> <xp> <yp> <ps>
    # Monocular, with velocity
    #        <time> <xp> <yp> <ps> <xv> <yv>
    # Monocular, with resolution
    #     <time> <xp> <yp> <ps> <xr> <yr>
    # Monocular, with velocity and resolution
    #     <time> <xp> <yp> <ps> <xv> <yv> <xr> <yr>
    # Binocular
    #     <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr>
    # Binocular, with velocity
    #     <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr>
    # Binocular, with and resolution
    #     <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xr> <yr>
    # Binocular, with velocity and resolution
    #      <time> <xpl> <ypl> <psl> <xpr> <ypr> <psr> <xvl> <yvl> <xvr> <yvr> <xr> <yr>
    # -----------------------
    # XXX: for now, we'll only support monocular and binocular formats, will include others later
    
    #replace missing fields. with NaN XXX: should this NaN be user defined, instead of hardcoded?
    raw = re.sub('   .\t','  nan\t',raw)


    gaze = findall_loadtxt("^\d+[\t\ -]+.*",raw,gaze_cols,gaze_dtype) 
    #NOTE: If the eyetracking data contains calibrations, then saccade
    #and blinks times will be off. Time series assumes all data sampled
    #continuously, but no data is sent during a calibration. we'll take
    #the approach of figuring out when there's missing data and adding it
    #into the timeseries with zeros (so it stays one continuous thing)
    prev = 0
    tmp = {}
    for f in gaze.dtype.names:
        tmp[f] = np.array([])
    # tmp[0] will be our time - change its dtype to be unsigned int
    t = tmp["time"]
    t.dtype='uint64'
    missing_tstamp= np.array([])

    #Use the eyelink-reported samplerate
    dt = findall_loadtxt("RATE[\t ]*\d+.\d*",raw,(1,)) 
    dt = 1000/dt[0]

    #XXX: throw error if diff.gaze(['time']) is every either 0 or negative (samples repeated or out of order)
    for D in np.arange(len(gaze['time']))[np.diff(gaze['time']) != dt]:
        print "WARNING: discontinuity in eyetracker time series",
        print " at sample ",D,", time ",str(gaze['time'][D:D+2])
        missing_tstamp= np.concatenate((missing_tstamp,gaze['time'][D:D+2]))
        t = np.concatenate((t,gaze['time'][prev:D+1], 
                            np.arange(gaze['time'][D], 
                                      gaze['time'][D+1],
                                      dt,dtype='uint64')))
        # missing values stored as NaNs
        z = np.ones((gaze['time'][D+1]- gaze['time'][D])/dt) * np.nan
        # .names[1:] skips over the 'time' field
        for f in gaze.dtype.names[1:]:
            tmp[f] = np.concatenate((tmp[f],gaze[f][prev:D+1], z))
        prev = D+1

    # iterate over all fields 
    tmp['time'] = t
    for f in gaze.dtype.names:
        tmp[f] = np.concatenate((tmp[f],gaze[f][prev:]))

    gaze = np.zeros(len(tmp['time']),dtype=gaze_dtype)
    for f in gaze.dtype.names:
        gaze[f] = tmp[f]

    raw = re.sub("\n\d+.*","",raw) # get rid of lines which we've already
    raw= re.sub("\nMSG.*","", raw) # extracted

    # See Notes in the Eyelink class docstring
    # Basically, all we need to do is find those ESACC events which are
    # preceded by EBLINK events, and discard the entire saccade
    
    #TODO: a better way of doing this is to separate the streams (as we do
    #later anyway) and call a function to parse the two streams seperately...
    #...but life's too short for now

    # convert endblinks to fixed length (so we can use lookback)
    blinks_r = findall_loadtxt("(?<=EBLINK.)R\ \d+\t\d+",raw,(1,2)) 
    blinks_l = findall_loadtxt("(?<=EBLINK.)L\ \d+\t\d+",raw,(1,2)) 
    blinks_r.shape = -1,2
    blinks_l.shape = -1,2


    #XXX: previous approach does not worksince R and L can be staggered relative to one another
    # soln: split the left and right streams, and process them individually
    def grab_raw_monocular(raw,eye='R'):
        """ this function removes any lines that pertain to the oposite eye"""
        raw_ret = raw
        omit = " L"
        if eye=='L':
            omit = " R"
        for x in ["EBLINK","ESACC","EFIX"]:
            raw_ret = re.sub(x+omit+".*\n","",raw_ret)

        return raw_ret

    raw_l = grab_raw_monocular(raw,'L')
    raw_r = grab_raw_monocular(raw,'R')
    raw_l= re.sub("EBLINK.*","EBLINK",raw_l)
    raw_r= re.sub("EBLINK.*","EBLINK",raw_r)
    #1/0
    # lookback and return endsaccades which are preceded by an endblink
    discard_r = findall_loadtxt("(?<=EBLINK\nESACC...).*\d+\t\d+",raw_r,(0,1)) 
    discard_l = findall_loadtxt("(?<=EBLINK\nESACC...).*\d+\t\d+",raw_l,(0,1)) 
    # XXX: separate timestamp gaps with blinks (and maybe have a method that
    # reports the OR of all the crap
    discard_r = np.append(missing_tstamp,discard_r)
    discard_l = np.append(missing_tstamp,discard_l)
    discard_l.shape = -1,2
    discard_r.shape = -1,2

    #get rid of lines containing falsely reported ESACCS which were preceded by EBLINK
    # XXX: now I'm getting paranoid - are we excluding more ESACCs than we
    # should if there were multiple blinks ?
    raw_r = re.sub("(?<=EBLINK)\nESACC.*","",raw_r)
    raw_l = re.sub("(?<=EBLINK)\nESACC.*","",raw_l)
    # lookback and return endsaccades (see DATA NOTATIONS for <codes>)
    #ESACC   <eye> <stime> <etime> <dur> <sxp> <syp> <exp> <eyp> <ampl> <pv>

    saccades_r = findall_loadtxt("(?<=ESACC...).*",raw_r,(0,1,7,8,3,4,5,6)).astype(float)
    saccades_l = findall_loadtxt("(?<=ESACC...).*",raw_l,(0,1,7,8,3,4,5,6)).astype(float)
    saccades_r.shape = -1,8
    saccades_l.shape = -1,8

    # (see DATA NOTATIONS for <codes>)
    #EFIX <eye> <stime> <etime> <dur> <axp> <ayp> <aps>
    fix_r= findall_loadtxt("(?<=EFIX...).*",raw_r,(0,1,3,4,5)).astype(float)
    fix_l= findall_loadtxt("(?<=EFIX...).*",raw_l,(0,1,3,4,5)).astype(float)
    fix_r.shape = -1,5
    fix_l.shape = -1,5
    
    raw = re.sub("\n[ES]SACC.*","",raw)  # get rid of lines which we've already 
    raw = re.sub("\n[ES]BLINK.*","",raw) # extracted

    eye = Eyelink(filename, binocular, have_right, have_left)


    if binocular:
        # names[0] is 'time', everything else is x,y,pupA,
        for name in gaze.dtype.names[1:4]:
            eye.l[name] = Events(gaze['time']/1000.0,v=gaze[name]) # typecast as Events,
        # so far we've called right eye data x2,y2,pupA2 - rename these to x,y,pupA
        for name,field in zip(('x','y','pupA'),gaze.dtype.names[4:]):
            eye.r[name] = Events(gaze['time']/1000.0,v=gaze[field]) # typecast as Events,
    else:
        if have_right:
            for name in gaze.dtype.names[1:4]:
                eye.r[name] = Events(gaze['time']/1000.0,v=gaze[name]) # typecast as Events,
                eye.l[name] = None
        else:
            for name in gaze.dtype.names[1:4]:
                eye.l[name] = Events(gaze['time']/1000.0,v=gaze[name]) # typecast as Events,
                eye.r[name] = None


    # eye["x"] = TimeSeries(gaze['x'],
    #                       t0=gaze['time'][0]/1000.0,
    #                       samplingrate=500)#gaze['time'][1]-gaze['time'][0])
    # eye["y"] = TimeSeries(gaze['y'],
    #                       t0=gaze['time'][0]/1000.0,
    #                       samplingrate=500)#gaze['time'][1]-gaze['time'][0])
    # eye["pupA"] = TimeSeries(gaze['pupA'],
    #                          t0=gaze['time'][0]/1000.0,
    #                          samplingrate=500)#gaze['time'][1]-gaze['time'][0])

    #XXX: wrap this up into a loop for neatness / brevity
    if blinks_r.size:
        eye.r.blinks = Epochs(blinks_r[:,0]/1000.0,blinks_r[:,1]/1000.0)

    if blinks_l.size:
        eye.l.blinks = Epochs(blinks_l[:,0]/1000.0,blinks_l[:,1]/1000.0)

    if discard_l.size:
        eye.l.discard = Epochs(discard_l[:,0]/1000.0,discard_l[:,1]/1000.0)
    if discard_r.size:
        eye.r.discard = Epochs(discard_r[:,0]/1000.0,discard_r[:,1]/1000.0)

    if saccades_l.size:
        eye.l.saccades = Events(saccades_l[:,0]/1000.0, start=saccades_l[:,0]/1000.0,stop=saccades_l[:,1]/1000.0,
                amplitude=saccades_l[:,2], vpeak=saccades_l[:,3],
                #epochs=Epochs(saccades_l[:,0]/1000.0,saccades_l[:,1]/1000.0),
                xi=saccades_l[:,3], yi=saccades_l[:,4],
                xf=saccades_l[:,5], yf=saccades_l[:,6],
                )
        eye.l.sacepochs = Epochs(saccades_l[:,0]/1000.0,saccades_l[:,1]/1000.0)
    if saccades_r.size:
        eye.r.saccades = Events(saccades_r[:,0]/1000.0, start=saccades_r[:,0]/1000.0,stop=saccades_r[:,1]/1000.0,
                amplitude=saccades_r[:,2], vpeak=saccades_r[:,3],
                #epochs=Epochs(saccades_r[:,0]/1000.0,saccades_r[:,1]/1000.0),
                xi=saccades_r[:,3], yi=saccades_r[:,4],
                xf=saccades_r[:,5], yf=saccades_r[:,6],
                )
        eye.r.sacepochs = Epochs(saccades_r[:,0]/1000.0,saccades_r[:,1]/1000.0)
    if fix_l.size:
        eye.l.fixations = Events(fix_l[:,0]/1000.0, start=fix_l[:,0]/1000.0,stop=fix_l[:,1]/1000.0,
                xavg=fix_l[:,2], yavg=fix_l[:,3], pavg=fix_l[:,4]
                #epochs=Epochs(fix_l[:,0]/1000.0,fix_l[:,1]/1000.0))
                )
        eye.l.fixepochs=Epochs(fix_l[:,0]/1000.0,fix_l[:,1]/1000.0)
    if fix_r.size:
        eye.r.fixations = Events(fix_r[:,0]/1000.0, start=fix_r[:,0]/1000.0,stop=fix_r[:,1]/1000.0,
                xavg=fix_r[:,2], yavg=fix_r[:,3], pavg=fix_r[:,4],
                #epochs=Epochs(fix_r[:,0]/1000.0,fix_r[:,1]/1000.0)
                )
        eye.r.fixepochs=Epochs(fix_r[:,0]/1000.0,fix_r[:,1]/1000.0)
    
    if frames.size:
        eye["frames"] = Events(frames[:,0]/1000.0,v=frames[:,1])

    if link_fix.size:
        eye["link_fixations"] = Events(link_fix/1000.0)

    eye["msgs"] = msgsstr
    eye["raw"] = raw # contains all the lines which have not been processed
    if gcdisp.size:
        gcdisp =  gcdisp.view(np.recarray)
        eye["gcdisp"] = Events(gcdisp.time/1000.0, x=gcdisp.x, y=gcdisp.y)

    # XXX: kind of wish I could do some interval math on epochs (union,
    # intersection, a la pyinterval)

    #1/0
    return eye

def read_eyelink_cached(fname,d=None):
    """
    Read the asc file in `fname` into the dictionary `d` using
    read_eyelink(fname) and cache the results there. On
    subsequent calls, no reading is performed, and the cached results are
    returned. 
    """
    if d is None: return read_eyelink(fname)

    if d.has_key(fname) == False:
        d[fname] = read_eyelink(fname)
    else:
        print "Using cached version of ", fname
    return d[fname]
    
def ck_read_array(lines, dtype, separator='\t'):
    """ Read a list of line with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime

        modified from the scipy.io cookbook
    """
    # XXX: this is unnecessary and slow, should be removed
    cast = np.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in lines:
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return np.rec.array(data, dtype=dtype)

def plot_eyetrace(eye, xax,yax,option='full', ):
    xax.plot(eye.x,eye.x.v,label='X position')
    ts = eye.sacepochs['tstart']
    te = eye.sacepochs['tstop']
    x_ssl = ts.evt2idx(eye.x)
    x_esl = te.evt2idx(eye.x)
    y_ssl = ts.evt2idx(eye.y)
    y_esl = te.evt2idx(eye.y)
    xax.plot([ts[:],te[:]],[eye.x[x_ssl],eye.x[x_esl]],'rp-', linewidth=5, alpha =
            .5, label="saccades")
    yax.plot(eye.y,eye.y.v,label='Y position')
    yax.plot([ts[:],te[:]],[eye.y[y_ssl],eye.y[y_esl]], linewidth=5, alpha = .5, label="saccades")
    pass

class EyelinkReplayer(object):
    """ Class which implements the pylink API but plays back eyelink .asc files
    as though they are streaming in real time.

    This is useful for testing VisionEgg scripts which use the streaming data
    coming from the eyetracker, such as the gaze-contingent and the saccade
    prediction paradigms.

    Even though some of the methods should return "sample" type, we'll keep
    returning the reference to this object and implement all of the necessary
    methods as one flat object, for brevity.

    """
    def __init__(self, replayfile='', el=None, SAMPLE_TYPE=200):
        if replayfile=='' and el==None:
            return
        if el ==  None:
            el = read_eyelink(replayfile)
        self.el = el
        self.i=0
        self.bitflip=False


        self.x =  el.eye_used.x.v
        self.y =  el.eye_used.y.v
        self.t = el.eye_used.x #because x is an event array, indexing into it will just return single points in time

        self.SAMPLE_TYPE = SAMPLE_TYPE

    def startRecording(self, recSamples,recEvents,streamSamples,streamEvents):
        # I might have screwed up the order of parameters here
        return 0

    def getDataCount(self,samples=True):
        """Returns True every other time it is called (to give out one sample at
        a time
        """
        if self.bitflip:
            # get roughly 1000Hz replay
            time.sleep(.001) 
        else:
            # increment the internal index
            self.i += 1
        self.bitflip = not self.bitflip
        return self.bitflip

    def getNewestSample():
        # increment the internal index
        self.i += 1
        return self


    def getNextData(self):
        return self.SAMPLE_TYPE

    def getFloatData(self):
        return self
    
    def isRightSample(self):
        return self.el.have_right
        
    
    def isLeftSample(self):
        return self.el.have_left

    def getRightEye(self):
        return self

    def getLeftEye(self):
        return self

    def getGaze(self):
        if self.i < len(self.x):
            return self.x[self.i], self.y[self.i]
        else:
            raise RuntimeError("Index %d exceeds  length of eyelink data (%d)"
                    % (self.i, len(self.x)))

    def getTime(self):
        # convert back to milliseconds, which is how pylink represents time
        return int(self.t[self.i]*1000) 

    def sendMessage(self, txt):
        print "Replayer got Eyelink MSG: '%s'"%txt
