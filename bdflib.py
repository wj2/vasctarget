"""
Library of functions to load in BioSemi bdf binary format

format description available at:
http://www.biosemi.com/faq/file_format.htm

"""

import numpy as np 
import struct 
import sys
import pylab as P
import matplotlib.pyplot as plt

import neuropy 
import data
import eyelib
from data import Epochs,Events,TimeSeries,Container
import datalib as DL

__all__ = ['Biosemi', 'BiosemiPaul','read_bdf', 'from_hcd_to_biosemi', 'make_plots', 'make_individual_plots','plot_epochs', 's2int']


class Biosemi(Container):
    def __init__(self,h,c,d,use_SI=True,ref=None):
        Container.__init__(self)
        for name in h.dtype.names:
            self[name] = str(h[name]).strip()

        self["map"] = {}    #integer index 
        for i,name in enumerate(c.label):
            self.map[name] = i
            # XXX: all channels are still unreferenced!!
            # we can either subtract one channel (Cz equiv) or we can subtract
            # off the average of some subset of channels
            if name == "Status":
                self[i] = TimeSeries(d[i].flatten().astype(np.uint16), 
                                     samplingrate=int(c.num_samples[i]),
                                     label=name,
                                     transducer=c.transducer[i],
                                     physdimension = c.physdimension[i],
                                     phys_min = int(c.phys_min[i]),
                                     phys_max = int(c.phys_max[i]),
                                     dig_min = int(c.dig_min[i]),
                                     dig_max = int(c.dig_max[i]),
                                     prefilter = c.prefilter[i],
                                     reserved = c.reserved[i],
                                     gain = float(int(c.phys_max[i])-int(c.phys_min[i]))/( int(c.dig_max[i]) - int(c.dig_min[i]))
                                     )
            else: 
                self[i] = TimeSeries(d[i].flatten().astype(float), 
                                     samplingrate=int(c.num_samples[i]),
                                     label=name,
                                     transducer=c.transducer[i],
                                     physdimension = c.physdimension[i],
                                     phys_min = int(c.phys_min[i]),
                                     phys_max = int(c.phys_max[i]),
                                     dig_min = int(c.dig_min[i]),
                                     dig_max = int(c.dig_max[i]),
                                     prefilter = c.prefilter[i],
                                     reserved = c.reserved[i],
                                     gain = float(int(c.phys_max[i])-int(c.phys_min[i]))/( int(c.dig_max[i]) - int(c.dig_min[i]))
                                    )
                if use_SI:
                    self[i] = self.to_SI(i)
            self[name] = self[i]
        #reverse mapping
        self["rmap"] = dict(zip(self.map.values(),self.map.keys()))

        self['num_channels']=self.count_channels()

        if self.map.has_key("EyeX"):
            self.map["x"] = self.map["EyeX"]
            self.map["y"] = self.map["EyeY"]
            self.map["p"] = self.map["PupA"]
        self.reference(ref)
        #self.subtract_mean()
        #self.min_to_zero()

    def count_channels(self):
        cnt = 0
        for li in self.map.keys():
            if li[0].isalpha() and li[1:].isdigit():
                cnt+=1
        return cnt

    def subtract_mean(self):
        for x in range(self.num_channels):
            self[x] -= self[x].mean()

    def min_to_zero(self):
        for x in range(self.num_channels):
            self[x] -= self[x].min()

    def reference(self,ref=None):
        # crude referencing below
        #  average of all channels
        r = np.zeros_like(self[0])
        if ref == None:
            print "EEG Reference: Using unreferenced signal"
        elif ref == 'all' or ref == 'mean':
            print "EEG Reference: Referencing signal with respect to all %d channels" % self.num_channels
            for x in range(self.num_channels):
                r += self[x]
            r /= self.num_channels
            ref='mean'
        elif ref == 'cz':
            print "EEG Reference: Referencing signal with respect to Cz (B16) channel"
            # old Cz only ref
            r=self.b16.copy() #must copy, otherwise r==0 after b16 in self[x]-= r
        else:
            try:
                r = self[ref].copy()
                print "EEG Reference: Referencing signal with respect to (%s) channel"% ref
            except KeyError:
                print "EEG Reference: %s is an unknown reference configuration" % ref
                print "EEG Reference: Using unreferenced signal"
        for x in range(self.num_channels):
            self[x] -= r
            self[x].i.ref = ref

    def to_SI(self, i):
        return DL.raw_to_SI(self[i].astype(float), 
                            scalefactor = self[i].i.gain / 1e6  #makes nV into V
                            ,
                            offset=0.0,
                            overwrite=True)
    # alt scale factor
    # float((self[i].i.phys_max-self[i].i.phys_min ))/( self[i].i.dig_max-self[i].i.dig_min)
    def fix_raw_eyetracking(self):
        k = eyelib.analog_to_angles(max(diff(self.x)))


    # XXX: this is horribly ugly and causes me all sorts of grief, need to
    # abandon this for something cleaner
    def __getattr__(self,key):
        if self.has_key(key):
            return self[key]
        if type(key) == str:
            if self.map.has_key(key.upper()): #works for a27 -> A27
                key = self.map[key.upper()]
            if self.map.has_key(key):
                key = self.map[key]
        return self[key]

    def __setattr__(self,key,value):
        if self.has_key(key):
            self[key] = value
        if type(key) == str:
            if self.map.has_key(key.upper()): #works for a27 -> A27
                key = self.map[key.upper()]
            if self.map.has_key(key):
                key = self.map[key]
        self[key] = value

    def plot_mean(self,x,y,color='k',linewidth=1):
        return neuropy.plot_avg_timeseries(y, SEM=True);

    def plot_psd(self,x,y,color='k',linewidth=1):
        return P.psd(y.mean(axis=0),NFFT=256,Fs=y.i.samplingrate)
    
    def plot_specgram(self,x,y,color='k',linewidth=1):
        return P.specgram(y.mean(axis=0),NFFT=256,Fs=y.i.samplingrate)

    def plot_epochs(self,channels,epochs,fn=plot_mean):
        print " %s channels will be displayed" % channels
        pldim = np.sqrt(len(channels)).astype(int)
        P.title("Analog Eye Tracking Data:" + self.subject_id )
        P.figtext(.5,0,self.record_id, horizontalalignment='center')
        P.figtext(.5,1,"Analog Eye Tracking Data:" +
                self.subject_id + " ("+ str(len(epochs)) + " epochs)"
                ,horizontalalignment='center',
                verticalalignment='top' )
        for i,ch in enumerate(channels):
            P.title(self[ch].i.label)
            P.subplot(pldim+1,pldim+1, i+1)
            eeg = self[ch].flatten()
            peri_sac_eeg = eeg[epochs]
            #P.subplot(313)
            #r = [np.random.randint(eeg.shape[0]) for x in epochs]  #grab random time points 
            #rand_epochs = neuropy.epochs_from_onsets(eeg.time()[r],2.,offset=.5)
            #peri_rand_eeg = eeg[rand_epochs]
            l1 = fn(peri_sac_eeg.time(),peri_sac_eeg,'k')

            #l2 = P.plot(peri_sac_eeg.time(),np.mean(peri_rand_eeg,axis=0),'b',linewidth=2)

        #P.legend(("peri_sac_eeg (%s)"%epochs.shape[0],"peri_rand_eeg (%s)"%epochs.shape[0]))
        #P.figlegend((l1,l2),("peri_sac_eeg (%s)"%epochs.shape[0],"peri_rand_eeg (%s)"%epochs.shape[0]), 'lower right')

        P.draw()
        return peri_sac_eeg
        #P.show()

class BiosemiPaul(Biosemi):
    SWEEP = 0x4000
    MAXPOST = 0x8000-1
    BLANKSWEEP= 0x3ffe
    def __init__(self,h,c,d,use_SI=True,ref=None):
        """Initialization for the custom Biosemi class that implements aspects
        unique to Paul and Tim's work, and not Biosemi reading in general.
        This includes vsync times and frame information which was sent to the
        Biosemi status line via a DAQ on the stimulus computer
        """
        print "using custom BiosemiPaul class"
        Biosemi.__init__(self,h,c,d,use_SI=use_SI,ref=ref)
        # times when vsync bit was high
        self.vsync = ((self.Status & 0x8000 ) >> 15).astype(bool)
        # rising edge of vsync bit (to avoid duplicate samples of the same vsync)
        vsync_mask = np.append(self.vsync[0], np.logical_and(self.vsync[1:],np.diff(self.vsync)))
        self.synctime = self.vsync.time()[vsync_mask]
        #self.synctime = synctime[np.diff(synctime) > .001]
        # There is a fixed one frame latency between posting of the frame index
        # on the DT340 board and the onset of the photodiode transient  (i.e.
        # 5ms + distance down screen).  Will change the reading/parsing of
        # frame timestamps to reflect the *actual* time the frame is displayed.
        # Also adding a four BioSemi sample jitter 'fudgefactor'
        # (because sometimes the status line changes after the vsync occurs,
        # despite the fact that the code only writes out to that port right
        # before a vsync). 2010-04-22: photodiode tests suggest that actually
        # sometimes the the status line changes before the next vsync, so I
        # changed the sign from .005-fudgefactor, to fudgefactor - This means
        # that the 1 frame delay is gone! that's weird!!!
        fudgefactor = 2./ self.a1.i.samplingrate
        #fudgefactor = 0
        frame = Events(self.synctime[:],
                v=self.Status[self.synctime-(fudgefactor)] & 0x7FFF)
        self.status = Events(self.Status.time(),
                v=self.Status & 0x7FFF)
        
        # status.v == SWEEP is sweep zero
        duringsweep = np.logical_and(frame.v >= self.SWEEP, frame.v < self.MAXPOST)
        duringblank= np.logical_or( frame.v == self.BLANKSWEEP, 
                frame.v == self.BLANKSWEEP + self.SWEEP)

        self.frame = frame[duringsweep]
        self.frame.v &= 0x3FFF

        #experiments should always be between all bits floating high (MAXPOST)
        # a diff on a binary array will be non-zero on all transitions
        experiment_marks=np.diff(frame.v==self.MAXPOST) > 0
        exp_times = np.array(frame[experiment_marks])
        exp_times.shape=(-1,2)
        self.experiment = neuropy.Epochs(exp_times.T[0], exp_times.T[1])


    def plot_raw_frameinfo(self):
        regular_status = self.status.v <  0x7FFF 
        top_plot=plt.subplot(211)
        plt.plot(self.status,self.status.v,'b-')
        plt.title("Plotting every sample (transitions >1 in red)")
        trans_not_one =  np.diff(self.status.v) > 1
        plt.plot(self.status[trans_not_one],self.status.v[trans_not_one],'ro')
        self.trans_not_one = trans_not_one
        plt.subplot(212,sharex=top_plot)
        plt.plot(self.status[regular_status],self.status.v[regular_status],'b')

    def plot_frameinfo(self):
        regular_frame = self.frame.v <  0x7FFF 
        top_plot=plt.subplot(211)
        plt.title("Plotting the frames (transitions >1 in red)")
        plt.plot(self.frame,self.frame.v,'b-')
        trans_not_one =  np.diff(self.frame.v) > 1
        plt.plot(self.frame[trans_not_one],self.frame.v[trans_not_one],'ro')
        self.trans_not_one = trans_not_one
        plt.subplot(212,sharex=top_plot)
        plt.scatter(self.frame[regular_frame],self.frame.v[regular_frame])

    def plot_vsync_hist(self):
        plt.subplot(111)
        histogram = plt.hist(np.diff(self.synctime))
    
    def plot_vsync_photodiod(self):
        ax1 = plt.subplot(211)
        #plt.plot(self.synctime, np.ones_like(self.synctime))
        plt.plot(self.vsync.time(), self.vsync*200)
        plt.plot(self.frame, self.frame.v)
        trans_one=  np.diff(self.frame.v) == 1
        regular_frame = self.frame.v <  0x7FFF 
        plt.scatter(self.frame[1:][trans_one],self.frame.v[1:][trans_one])
        plt.subplot(212,sharex=ax1)
        plt.plot(self.Ana1.time(), self.Ana1)
        
    def frame_lag_test(self,photo_ch,stim_start,stim_end,thresh):
        '''
        Generates histogram of time differences between repeated single frame
        flashes and the actual onset of the flash as measured by a photodiode
        photo_ch is the biosemi channel with the photodiode timeseries,
        
        'stim_start' & 'stim_end' delineate the start of the flash stimulus,
        'thresh' should be set to capture all of the photodiode transients.
        
        Use 'plot_vsync_photodiod' to determine appropriate function parameters.
        '''
         
        # detect photodiode flashes
        se = Epochs((stim_start, stim_end))
        x=neuropy.cut_time_series(self[photo_ch],se)
        a=(x[1:]-x[:-1])>0
        b=(x[2:]-x[1:-1])<0
        c=a[:-1] & b # mask with 'True' one sample to the left of local maxima
        d=x[1:-1]>thresh # threshold photodiode transients
        xx=P.find(c & d != 0) # maxima of photodiode transients

        pdt = xx/float(self[photo_ch].i.samplingrate)
        frt = self.status[P.find(np.diff(self.status.v)==1)+1]
        lag = pdt-frt+stim_start-0.001 #subtract 1ms to get flash onset rather than maxima
        plt.figure(); plt.hist(lag, label='from status line change to photodiode onset')
        frt = self.frame[P.find(np.diff(self.frame.v)==1)+1]
        lag = pdt-frt+stim_start-0.001 #subtract 1ms to get flash onset rather than maxima
        plt.hist(lag, label='from frame change to photodiode onset')
        plt.title("Histogram of delays"); plt.legend(loc='upper center')
        plt.grid('on'); plt.show()
        return lag
        

def from_hcd_to_biosemi(h,c,d):
    return Biosemi(h,c,d)


def make_plots(filename):
    b = read_bdf(filename)
    sac = eyelib.detect_saccades_trivial(b.p,b.x,b.y)
    epochs = neuropy.epochs_from_onsets(sac,2.0,offset=.5)
    P.clf()
    b.plot_epochs(range(32),epochs,fn=b.plot_mean)
    P.savefig('for_tim/'+filename.split('/')[-1].split('.')[0]+"_meanSEM_A1-32.png",dpi=400)
    P.clf()
    b.plot_epochs(range(32,64),epochs,fn=b.plot_mean)
    P.savefig('for_tim/'+filename.split('/')[-1].split('.')[0]+"_meanSEM_B1-32.png",dpi=400)

def make_individual_plots(filename):
    b = read_bdf(filename)
    sac = eyelib.detect_saccades_trivial(b.p,b.x,b.y)
    epochs = neuropy.epochs_from_onsets(sac,2.0,offset=.5)
    print " plots will be based on ", len(epochs), " epochs"
    for i in range(64):
        P.clf()
        peri_sac_eeg = b[i][epochs]
        P.specgram(peri_sac_eeg.mean(axis=0),NFFT=256,Fs=b[i].i.samplingrate)
        P.colorbar()
        P.savefig('for_tim/'+filename.split('/')[-1].split('.')[0]+"_spec_"+str(b.rmap[i])+".png",dpi=200)
        P.clf()
        neuropy.plot_avg_timeseries(peri_sac_eeg) #gaze['x'] is now TimeSeries
        P.savefig('for_tim/'+filename.split('/')[-1].split('.')[0]+"_mean_"+str(b.rmap[i])+".png",dpi=200)

    
def plot_epochs(h,c,d,channels,epochs):
    print " %s channels will be displayed" % channels
    pldim = np.sqrt(len(channels)).astype(int)
    #P.figure()
    P.title("Analog Eye Tracking Data:" + h.subject_id )
    P.figtext(.5,0,h.record_id, horizontalalignment='center')
    P.figtext(.5,1,"Analog Eye Tracking Data:" +
            h.subject_id,horizontalalignment='center',
            verticalalignment='top' )
    for i,ch in enumerate(channels):
        P.subplot(pldim+1,pldim+1, i+1)
        P.title(c.label[ch])
        eeg = d[ch].flatten()
        peri_sac_eeg = eeg[epochs]
        #P.subplot(312)
        #P.plot(peri_sac_eeg.time(),N.mean(peri_sac_eeg,axis=0),'k',linewidth=2)
        #P.plot(peri_sac_eeg.time(),N.min(peri_sac_eeg,axis=0),'k',linewidth=1)
        #P.plot(peri_sac_eeg.time(),N.max(peri_sac_eeg,axis=0),'k',linewidth=1)
        #P.errorbar(peri_sac_eeg.time(),N.mean(peri_sac_eeg,axis=0),yerr=peri_sac_eeg.std(axis=0),color=(.1,.1,.1,.1),linewidth=1)

        #P.subplot(313)
        r = [np.random.randint(eeg.shape[0]) for x in epochs]  #grab random time points 
        rand_epochs = neuropy.epochs_from_onsets(eeg.time()[r],2.,offset=.5)
        peri_rand_eeg = eeg[rand_epochs]
        l1 = P.plot(peri_sac_eeg.time(),np.mean(peri_sac_eeg,axis=0),'k',linewidth=1)

        #l2 = P.plot(peri_sac_eeg.time(),np.mean(peri_rand_eeg,axis=0),'b',linewidth=2)

    #P.legend(("peri_sac_eeg (%s)"%epochs.shape[0],"peri_rand_eeg (%s)"%epochs.shape[0]))
    #P.figlegend((l1,l2),("peri_sac_eeg (%s)"%epochs.shape[0],"peri_rand_eeg (%s)"%epochs.shape[0]), 'lower right')

    P.draw()
    #P.show()


def s2int(t):
    """ convert tuple of stringable things to integers"""
    r = []
    for x in t:
        r.append(int(str(x)))
    return r

def read_bdf(filename='timcal.bdf', legacy=None, use_SI=True, ref=None,
        Biosemi=BiosemiPaul):
    """
    loads data classes from BioSemi bdf binary file
    XXX: not actually reading in data classes yet, just raw data
    """
    try:
        f = file(filename)
    except IOError:
        sys.stderr.write("Could not open file '"+filename+"'");
        return None;


    # XXX: it's not clear to me if i should be 
    # using dtype stuff or struct unpack stuff
    id_code 			= np.dtype('S8')
    subject_id 			= np.dtype('S80')
    record_id 			= np.dtype('S80') 
    startdate			= np.dtype('S8')
    starttime			= np.dtype('S8')
    header_sz 			= np.dtype('S8') 
    data_format 		= np.dtype('S44')
    num_records 		= np.dtype('S8')
    duration_record		= np.dtype('S8')
    num_channels 		= np.dtype('S4')

    file_hdr = np.dtype(
            {
                'names': [
                    'id_code'
                    ,'subject_id'
                    ,'record_id'
                    ,'startdate'
                    ,'starttime'
                    ,'header_sz'
                    ,'data_format'
                    ,'num_records'
                    ,'duration_record'
                    ,'num_channels'
                    ],
                'formats': [ 
                    id_code 
                    ,subject_id 
                    ,record_id  
                    ,startdate
                    ,starttime
                    ,header_sz 
                    ,data_format 
                    ,num_records 
                    ,duration_record
                    ,num_channels 
                    ]
                }
            )
    # total record length is num_records * duration_record 
        
    f.seek(0)
    #test = fromstring(f.read(256), dtype=file_hdr)
    test = np.array(f.read(file_hdr.itemsize), dtype=file_hdr)

    N = int('17')
    N = int(str(test['num_channels']))

    label 		    	= np.dtype(('S16',N))
    transducer 			= np.dtype(('S80',N))
    physdimension 		= np.dtype(('S8',N))
    phys_min 			= np.dtype(('S8',N))
    phys_max 			= np.dtype(('S8',N))
    dig_min 			= np.dtype(('S8',N))
    dig_max 			= np.dtype(('S8',N))
    prefilter 			= np.dtype(('S80',N))
    num_samples         = np.dtype(('S8',N))
    reserved 			= np.dtype(('S32',N))

    ch_hdr = np.dtype(
            {
                'names': [ 
                    'label',
                    'transducer',
                    'physdimension',
                    'phys_min',
                    'phys_max',
                    'dig_min',
                    'dig_max',
                    'prefilter',
                    'num_samples',
                    'reserved'
                    ],
                'formats': [
                    label ,
                    transducer ,
                    physdimension ,
                    phys_min ,
                    phys_max ,
                    dig_min ,
                    dig_max ,
                    prefilter ,
                    num_samples,
                    reserved ,
                    ]
                }
            )
    chtest =  np.array(f.read(ch_hdr.itemsize), dtype=ch_hdr)

    hdr = test.view(np.recarray)
    chl = chtest.view(np.recarray)


    chl.num_samples

    rawdata = f.read() #header done, read the rest to EOF 
    f.close()
    dtmp = np.fromstring(rawdata,dtype='B')
    d = np.zeros(len(rawdata) + len(rawdata)/3,'B')
    if np.array(1,'i') != np.array(1,'<i'): #if not little endian machine
        dtmp = dtmp.byteswap()
    d[0::4] = dtmp[0::3]
    d[1::4] = dtmp[1::3]
    d[2::4] = dtmp[2::3]
    d[3::4][ dtmp[2::3] > 127  ] = 255
    # XXX: not sure if this works as is on big endian machines, might need another byteswap here
    d.dtype = 'i'
    # XXX: this assumes that all channels are sampled at the same rate, which is not necessarily true
    d.shape = s2int( (hdr.num_records, hdr.num_channels, chl.num_samples[0], ))
    d = d.transpose(1,0,2)  # move the channels to be first dimension
    # XXX: we have to .copy() d in order to do the next step, since array is no longer contiguous. Alternatively, we'll have to call .flatten() onevery channel (see examples/plot_eeg.py) 
    #d.shape = d.shape[0],d.shape[1]*d.shape[2] # records are meaningless, all we have are samples


    d = data.TimeSeries(d, samplingrate=int(chl.num_samples[0]))

    #XXX: cache the resulting parsed file, as this takes some time
    if legacy:
        return hdr, chl, d

    return Biosemi(hdr, chl, d, use_SI, ref)
