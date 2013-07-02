"""
This module contains customized classes and methods
derived from data.py and datalib.py.
"""
import numpy as N
import data as D
import datalib as DL
import utils as U
import os
import re

from neuropy import datapath
try: from neuropy import datapath_usreylab as datapath
except: pass

from neuropy import cachepath

__all__ = ['OpticTractRecording','PairedCellRecording']

class OpticTractRecording(D.Recording):
    """
    Container class for optic tract single cell recordings.

    Inputs:
      celldate - a string specifying the recording date, e.g. '20050803'
      cellid   - a string specifying the cell, e.g. 'F01'
    """
    def __init__(self,celldate=None,cellid=None,ch_ret=None):
        self.i = D.Info()
        self.i.celldate = celldate
        self.i.cellid = cellid
        self.i.ch_ret = ch_ret
        self.i.t_frame = 0.007139184
        
        # set path variables
        self.i.data_path = os.path.join(datapath)
        self.i.rec_path = os.path.join(self.i.data_path,'OpticTractData',celldate)
        self.i.rec_filenames = U.parse_dir(self.i.rec_path,
                                           suffix='mat',key_type=str)
        self.i.movie_path = os.path.join(self.i.data_path,'stimulus')
        # find recording id
        rec_ids = filter(lambda k: cellid in k, self.i.rec_filenames.keys())
        if len(rec_ids) != 1:
            print '%d recordings found'%len(rec_ids)
            self.i.rec_id = None
        else:
            self.i.rec_id = rec_ids[0]
        
        # create empty containers for spikes, lfps, etc.
        # self.spks=D.Container() # dict contains spike timings
        # self.frm=D.Container() # dict contains auditory stimuli

        if self.i.rec_id != None:
            self.load()
        
    def load(self):
        '''
        loads spike and stimulus data

        '''
        import scipy.io as io
        fn = os.path.join(self.i.rec_path,self.i.rec_filenames[self.i.rec_id])
        mat = io.loadmat(fn)
        self.mat = mat

        info = D.Info()
        for key in ['__header__','__version__']:
            info[key] = mat[key]
        self.i.update(info)
        
        self.i.event_ids = filter(lambda k: self.i.rec_id in k, mat.keys())
        for event_id in self.i.event_ids:
            evt = D.Events(mat[event_id].times)
            evt.i.units = 's'
            # retina
            if self.i.ch_ret != None:
                if re.match('.*Ch%d$'%self.i.ch_ret,event_id):
                    self.ret = evt
            else:
                if 'Cortex' in mat[event_id].title:
                    self.ret = evt
            # stimulus
            if 'Trigger' in mat[event_id].title:
                self.frm = evt
                # fix trigger in case of m_sequence
                if len(self.frm) == 327669: # 10 m-sequences
                    print 'adding trigger t0'
                    self.frm = D.Events(N.concatenate(([0],evt)),i=evt.i)
                if len(self.frm) == 10:
                    print 'inserting triggers'
                    t_frame = self.i.t_frame
                    self.frm = D.Events(N.concatenate([N.arange(32767)*t_frame+ti for ti in evt]),i=evt.i)
                self.i.start_times = self.frm[::32767]
                self.i.t_trial = N.diff(self.i.start_times)[0]
                self.i.epochs = DL.epochs_from_onsets(self.i.start_times,self.i.t_trial)

    def load_movie(self, filename='Big_M.mat', varname='mseq'):
        '''
        loads a stimulus movie
        '''
        return load_movie(self, filename=filename, varname=varname)


class PairedCellRecording(D.Recording):
    """
    Container class for Retina-LGN paired cell recordings.

    Inputs:
      celldate - a string specifying the recording date, e.g. '20050803'
      cellid   - a string specifying the cell, e.g. 'F01'
    """
    def __init__(self,celldate=None,cellid=None,ch_ret=None,ch_lgn=None):
        self.i = D.Info()
        self.i.celldate = celldate
        self.i.cellid = cellid
        self.i.ch_ret = ch_ret
        self.i.ch_lgn = ch_lgn
        self.i.t_frame = 0.007139184
        # set path variables
        self.i.rec_path = os.path.join(datapath,'RetinaLGNPairData',celldate)
        self.i.rec_filenames = U.parse_dir(self.i.rec_path,
                                           suffix='mat',key_type=str)        
        # find recording id
        rec_ids = filter(lambda k: cellid in k, self.i.rec_filenames.keys())
        if len(rec_ids) != 1:
            print '%d recordings found'%len(rec_ids)
            self.i.rec_id = None
        else:
            self.i.rec_id = rec_ids[0]
        
        # create empty containers for spikes, lfps, etc.
        # self.spks=D.Container() # dict contains spike timings
        # self.frm=D.Container() # dict contains auditory stimuli

        if self.i.rec_id != None:
            self.load()
        
    def load(self):
        '''
        loads spike and stimulus data

        '''
        import scipy.io as io
        fn = os.path.join(self.i.rec_path,self.i.rec_filenames[self.i.rec_id])
        mat = io.loadmat(fn)
        # self.mat = mat

        info = D.Info()
        for key in ['__header__','__version__']:
            info[key] = mat[key]
        self.i.update(info)
        
        self.i.event_ids = filter(lambda k: self.i.rec_id in k, mat.keys())
        for event_id in self.i.event_ids:
            evt = D.Events(mat[event_id].times)
            evt.i.units = 's'
            # lgn
            if self.i.ch_lgn != None:
                if re.match('.*Ch%d$'%self.i.ch_lgn,event_id):
                    self.lgn = evt
            else:
                if 'LGN' in mat[event_id].title:
                    self.lgn = evt
            # retina
            if self.i.ch_ret != None:
                if re.match('.*Ch%d$'%self.i.ch_ret,event_id):
                    self.ret = evt
            else:
                if 'Cortex' in mat[event_id].title:
                    self.ret = evt
            # stimulus
            if 'Trigger' in mat[event_id].title:
                self.frm = evt


def load_movie(rec, filename=None, varname=None):
    # load m-sequence movie
    try:    # data already loaded?
        movie = rec.movie
    except: # load data
        import scipy.io as io
        fn_movie = os.path.join(rec.i.movie_path,filename)
        mat = io.loadmat(fn_movie)
        # self.mat = mat

        info = D.Info()
        for key in ['__header__','__version__']:
            info[key] = mat[key]
        
        movie = D.Image(mat[varname])
        nframes = movie.shape[0]
        if len(movie.shape) == 2:
            nx = int(N.sqrt(movie.shape[1]))
            movie = N.reshape(movie,(nframes,nx,nx))
        nframes,nx,ny = movie.shape
        
        # update movie info
        movie.i.update(info)
        movie.i.nx = nx
        movie.i.ny = ny
        movie.i.nframes = nframes
        movie.i.filename = fn_movie

    rec.movie = movie
    
    return movie


if __name__ == "__main__":
    celldate = '20021016'
    cellid = 'a04b'
    rec = OpticTractRecording(celldate=celldate,cellid=cellid)
    print 'single cell recording:'
    print '%d retina spikes, %d frames'%(len(rec.ret),len(rec.frm))

    celldate = '20050803'
    cellid = 'F01'
    rec = PairedCellRecording(celldate=celldate,cellid=cellid)
    print '\npaired recording:'
    print '%d retina spikes, %d lgn spikes, %d frames'%(len(rec.ret),len(rec.lgn),len(rec.frm))
