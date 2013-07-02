"""
This module contains customized classes and methods
derived from data.py and datalib.py.
"""
import numpy as N
import data as D
import dsplib as DSP
import datalib as DL
import spikelib as SL
import plotlib as PL
import utils as U
import os
import struct
import re

from neuropy import datapath
from neuropy import cachepath

__all__ = ['CatRecording','MonkeyRecording','denoise_lfp','compute_strf']

# translation table for variable names
translation = {'REFRESHRATE':'refreshrate'}

                        
class CatRecording(D.Recording):
    """
    Container class for Tim's cat recordings.

    Inputs:
      cat_id   - a string specifying the cat,  e.g. '15'
      track_id - a string specifying the recording track,  e.g. '7c'
      rec_id   - a string specifying the recording file, e.g. '83'
      cond_id  - an integer specifying the stimulus condition
    """
    def __init__(self,cat_id,track_id,rec_id,cond_id=0):
	info = D.Info()
        info.cond_id = cond_id  # stimulus condition
        info.cat_id = cat_id
        info.track_id = track_id
        info.rec_id = rec_id
        info.cond_id = cond_id

        # set path variables
        info.rec_path = os.path.join(datapath,'cat_%s','track_%s','rec_%s')%(
            cat_id,track_id,rec_id)
        info.rec_name = 'cat_%s-track_%s-rec_%s-%d'%(
            cat_id,track_id,rec_id,cond_id) #unique name
        info.infofile = os.path.join(info.rec_path,'rec_info.inf')
        info.cache_path = os.path.join(cachepath,info.rec_name)
        info.spk_path = os.path.join(info.rec_path,'spk_data')
        info.lfp_path = os.path.join(info.rec_path,'lfp_data')
        info.stim_path = os.path.join(info.rec_path,'stim_data')
        info.fig_path = os.path.join(info.rec_path,'figures')
        info.res_path = os.path.join(info.rec_path,'results')

        # parse rec_info file
        self.i = D.Info(filename=info.infofile,translation=translation)
        # add path variables
        self.i.update(info)
        if hasattr(self.i,'rec_epochs'):
            (self.i.t_start,self.i.t_stop) = self.i.rec_epochs[cond_id]

        # create empty containers for spikes, lfps, etc.
        self.spk=D.Container() # dict contains spike timings
        self.lfp=D.Container() # dict contains LFPs
        self.comments={} # dict contains comments during data collection
        
        # intialize cache
        try:            
            U.make_dir(cachepath)
        except:
            assert 0, "could not create cache dir %s"%cachepath

    def load_lfp(self,id_lfp=None):
        '''
        loads local field potential data

        Input:
          id_lfp - an id that identifies the LFP recording
        '''
        return load_lfp(self,id_lfp)

    def load_spk(self,id_spk=None):
        '''
        loads time-stamps for spikes from a .spk file and convert to seconds.
        
        Input:
          id_spk - int that identifies the LFP recording     
        '''
        return load_spk(self,id_spk)

    def load_stim(self):
        '''
        loads stimulus timing data          
        '''
        return load_stim(self)

    def load_movie(self, **args): #filename=None, invert=False, dtype=N.uint8):
        '''
        loads a stimulus movie

        Input:
          filename - the file containing the raw movie data, if blank
                     trys the filename specified in: self.stim.i.fname 
        '''
        return load_movie(self, **args) #filename=filename, invert=invert, dtype=dtype)


class MonkeyRecording(D.Recording):
    """
    Container class for Tim's cat recordings.

    Inputs:
      monkey_id   - a string specifying the monkey,  e.g. '15'
      track_id - a string specifying the recording track,  e.g. '7c'
      rec_id   - a string specifying the recording file, e.g. '83'
      cond_id  - an integer specifying the stimulus condition
    """
    def __init__(self,monkey_id,track_id,rec_id,cond_id=0):        
        info = D.Info()
        info.cond_id = cond_id  # stimulus condition
        info.monkey_id = monkey_id
        info.track_id = track_id
        info.rec_id = rec_id
        info.cond_id = cond_id

        # set path variables
        info.rec_path = os.path.join(datapath,'monkey_%s','track_%s','rec_%s')%(
            monkey_id,track_id,rec_id)
        info.rec_name = 'cat_%s-track_%s-rec_%s-%d'%(
            monkey_id,track_id,rec_id,cond_id) #unique name
        info.infofile = os.path.join(info.rec_path,'rec_info.inf')
        info.cache_path = os.path.join(cachepath,info.rec_name)
        info.spk_path = os.path.join(info.rec_path,'spk_data')
        info.lfp_path = os.path.join(info.rec_path,'lfp_data')
        info.stim_path = os.path.join(info.rec_path,'stim_data')
        info.fig_path = os.path.join(info.rec_path,'figures')
        info.res_path = os.path.join(info.rec_path,'results')

        # parse rec_info file
        self.i = D.Info(filename=info.infofile,translation=translation)
        # add path variables
        self.i.update(info)
        if hasattr(self.i,'rec_epochs'):
            (self.i.t_start,self.i.t_stop) = self.i.rec_epochs[cond_id]

        # create empty containers for spikes, lfps, etc.
        self.spk=D.Container() # dict contains spike timings
        self.lfp=D.Container() # dict contains LFPs
        self.comments={} # dict contains comments during data collection
        
        # intialize cache
        try:            
            U.make_dir(cachepath)
        except:
            assert 0, "could not create cache dir %s"%cachepath
            
    def load_lfp(self,id_lfp=None):
        '''
        loads LFP data

        Input:
          id_lfp - an id that identifies the LFP recording
        '''
        return load_lfp(self,id_lfp)
    
    def load_spk(self,id_spk=None):
        '''
        Load time-stamps for spikes from a .spk file and convert to seconds.
        
        Input:
          id_spk - int that identifies the LFP recording
          
        '''
        return load_spk(self,id_spk)

def load_lfp(rec,id_lfp=None):
    try:
        filenames = rec.lfp.filenames
        info = rec.lfp.i
    except: # parse directory
        infofile = os.path.join(rec.i.lfp_path,'lfp_info.inf')
        info = D.Info(filename=infofile,translation=translation)
        # default parameters
        params = dict(filename_prefix=None,filename_suffix='lfp')
        params.update(info)
        # get filenames
        filenames = U.parse_dir(rec.i.lfp_path,
                                prefix=params['filename_prefix'],
                                suffix=params['filename_suffix'])
        rec.lfp.filenames = filenames
        rec.lfp.i = info
        
    if id_lfp != None:
	try:    # data already loaded?
            lfp = rec.lfp[id_lfp]
        except: # load data
            fn_lfp = os.path.join(rec.i.lfp_path,filenames[id_lfp])            

            # load raw data
            lfp = DL.load_time_series(fn_lfp,samplingrate=info.samplingrate,
                                      dtype=info.datatype,
                                      nsamples=info['nsamples'])

            # convert data into SI units
            if hasattr(info,'units'):
                zero = 2**(info.resolution-1)
                lfp_key = re.sub('\.\w+$','',filenames[id_lfp])
                if isinstance(info.units_multiplier,dict):
                    units_multiplier = info.units_multiplier[lfp_key]
                else:
                    units_multiplier = info.units_multiplier
                DL.raw_to_SI(lfp,scalefactor=units_multiplier,offset=zero)
            else:
                info.units = 'raw'
                
            # cut lfp
            if (hasattr(rec.i,'t_start') & hasattr(rec.i,'t_stop')):
                lfp = lfp[D.Epochs((rec.i.t_start,rec.i.t_stop))]
            lfp.i.update(info)
            lfp.i.id_lfp=id_lfp
            if info.units != 'raw':
                lfp.i.units_multiplier = units_multiplier
            rec.lfp[id_lfp] = lfp
            return lfp
    else: # load all of them
        id_lfp = filenames.keys()
        for i in id_lfp:
            rec.load_lfp(i)

def load_spk(rec,id_spk):
    try:
        filenames = rec.spk.filenames
        info = rec.spk.i
    except: # parse directory
        infofile = os.path.join(rec.i.spk_path,'spk_info.inf')
        info = D.Info(filename=infofile,translation=translation)
        # default parameters
        params = dict(filename_prefix=None,filename_suffix='spk')
        params.update(info)
        # get filenames
        filenames = U.parse_dir(rec.i.spk_path,
                                prefix=params['filename_prefix'],
                                suffix=params['filename_suffix'])
        rec.spk.i = info
        rec.spk.i.update(filenames)
        
    if id_spk != None:
        try:    # data already loaded?
            spk = rec.spk[id_spk]
        except: # load data
            fn_spk = os.path.join(rec.i.spk_path,filenames[id_spk])

            # load raw data
            spk = DL.load_event_timings(fn_spk,dtype=info.datatype)
            
            # convert data into SI units
            if hasattr(info,'units'):
                DL.raw_to_SI(spk,scalefactor=info.units_multiplier,offset=0)
            else:
                info.units = 'raw'

            # cut spk
            if (hasattr(rec.i,'t_start') & hasattr(rec.i,'t_stop')):
                spk = spk[D.Epochs((rec.i.t_start,rec.i.t_stop))]
            spk.i.update(info)
            spk.i.id = id_spk
            rec.spk[id_spk] = spk
        return spk
    else: # load all of them
        id_spk = filenames.keys()
        for i in id_spk:
            rec.load_spk(i)

def load_stim(rec):
    try:    # data already loaded?
        stim = rec.stim
    except :# load data
        infofile = os.path.join(rec.i.stim_path,'stim_info.inf')
        info = D.Info(filename=infofile,translation=translation)

        # load and parse meta data
        filename_prefix = info.get('filename_prefix','.+\.s') # load all stimuli if no prefix given
        infofiles = U.parse_dir(rec.i.stim_path,prefix=filename_prefix,suffix='inf')
        infofile = os.path.join(rec.i.stim_path,infofiles[rec.i.cond_id])
        moreinfo = D.Info(filename=infofile,translation=translation)
        info.update(moreinfo)

        #
        # this is incomplete; the parser has to deal with vision egg files
        # and extract i.t_trial, i.t_start, i.start_times
        #
        # start_times = frm.select(rec.info["frm_ind"],rec.info["frm_skip"])
        # if hasattr(rec.info,"t_trial"):
        #     if len(start_times) > 1:
        #         rec.info["t_trial"] = round(N.diff(start_times)[0])
        # rec.info["t_start"] = start_times[0]
        # rec.info["start_times"] = start_times

        # load raw stimulus data
        fn_stim = os.path.join(rec.i.stim_path,info.frm_filename)
        stim = D.Stimulus(DL.load_events(fn_stim,dtype_t=info.timestamp_datatype,
                                         dtype_v=info.value_datatype))
    
        # convert data into SI units
        if hasattr(info,'units'):
            DL.raw_to_SI(stim,scalefactor=info.units_multiplier,offset=0)
        else:
            info.units = 'raw'

        # cut din
        if (hasattr(rec.i,'t_start') & hasattr(rec.i,'t_stop')):
            stim = stim[D.Epochs((rec.i.t_start,rec.i.t_stop))]
        stim.i.update(info)
        stim.i.id = rec.i.cond_id

        # in case of movies with identical repeats: change frame index
        if hasattr(stim.i,'nrepeats'):
            if stim.i.nrepeats != stim.i.nruns:
                max_frm = stim.i.nsweeps/stim.i.nrepeats
                stim.v = N.mod(stim.v,max_frm)
        else:
            stim.i.nrepeats = stim.i.nruns
            
        # compute framerate
	if stim.i.has_key('sweeptimeMsec'): # temporary - all exported times should be in s
	    stim.i.framerate = 1000./stim.i.sweeptimeMsec
	else:
	    stim.i.framerate = 1./stim.i.sweeptimeSec 

        rec.stim = stim
    return stim

def load_movie(rec, filename=None, dtype=N.uint8, invert=False, normalize=False, **args):
    """
    TO DO: get additional movie information parsed from stim_info.inf
    """
    try:    # already in memory?
        movie = rec.movie
    except: # load data
        if filename == None:
            filename = os.path.join(rec.i.stim_path, rec.stim.i.fname)
        f = file(filename, 'rb')
        hdrstr = f.read(5)
        if hdrstr == 'movie': # parse movie header, if present
            hdrlen = 11
            nx, = struct.unpack('H', f.read(2))
            ny, = struct.unpack('H', f.read(2))
            nframes, = struct.unpack('H', f.read(2))
        else: # no header, use default values
            hdrlen = 0
            nx = ny = 64
            nframes = 6000
        f.close()
        
        # load raw movie data
        movie = DL.load_movie(filename, dtype, nx, ny, nframes, hdrlen)

        # process movie
        if normalize:
            if movie.dtype != N.uint8: raise 'normalize only works for 8 bit unsigned movies'
            movie = (movie.astype(float) -128) / 128 # pixel values from -1 to 1
        if invert:
            if movie.dtype == N.uint8: raise 'inversion only works for signed value movies'
            movie = -movie
        
        # update movie info
        movie.i.nx = nx
        movie.i.ny = ny
        movie.i.nframes = nframes
        movie.i.filename = filename

    rec.movie = movie
    
    return movie

def denoise_time_series(rec, id_lfp=None):
    # load stimulus in order to get framerate etc.
    rec.load_stim()

    freqs_refresh = [rec.stim.i.refreshrate]
    print freqs_refresh
    freqs_frame = [rec.stim.i.framerate,2*rec.stim.i.framerate,3*rec.stim.i.framerate]
    print freqs_frame
    noise = D.Info(freqs_line=[60,120,180,240,300,360],sf_line=.075,
                   freqs_refresh=freqs_refresh,sf_refresh=.01,
                   freqs_frame=freqs_frame,sf_frame=.025)

    if id_lfp==None:
        ids = rec.lfp.keys() # denoise all lfp channels
    else:
        ids = [id_lfp]

    for id_lfp in ids:
        rec.lfp[id_lfp] = DSP.denoise_time_series(rec.lfp[id_lfp],**noise)

    if id_lfp != None:
        return rec.lfp[id_lfp]

def compute_strf(rec,id_spk=None,plot_strf=False,**args):
    """
    compute strf from recording rec (e.g. m-sequence)
    """
    
    # load stimulus information
    rec.load_stim()

    # load movie (e.g. m-sequence checker-board) stimulus
    rec.load_movie(args.get('invert',True)) # invert if req'd for correct STRF polarity

    # load neuron spiketimes
    if id_spk == None: 
        NotImplementedError # we need id_spk. later: implement strf for all cells
    spk = rec.load_spk(id_spk)

    # zoom in on ROI, specified in pixel coordinates
    roi = args.get('roi',None) # [xmin, xmax, ymin, ymax]

    # compute STA and plot STRF
    tslices = args.get('tslices',N.linspace(0, 28, 29)) # tslices specified in refresh interval units
    if plot_strf:
        strf = PL.plot_strf(spk, rec.stim, rec.movie, tslices=tslices, roi=roi, normalize=False)
    else:
        strf = SL.strf(spk, rec.stim, rec.movie, tslices=tslices, roi=roi, normalize=False)

    return strf
