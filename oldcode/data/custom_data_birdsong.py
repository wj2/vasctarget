"""
This module contains customized classes and methods
derived from data.py and datalib.py.
"""
import numpy as N
import data as D
import datalib as DL
import utils as U
import os

from neuropy import datapath
from neuropy import cachepath

__all__ = ['BirdRecording']

class BirdRecording(D.Recording):
    """
    Container class for birdsong recordings.

    Inputs:
      cellid   - a string specifying the cell, e.g. 'blabla0713_3_B_data'
    """
    def __init__(self,cellid=None):        
        self.i = D.Info()
        # set path variables
        self.i.rec_path = os.path.join(datapath,'birdsong','fieldL_data')
        self.i.rec_filenames = U.parse_dir(self.i.rec_path,
                                           suffix='mat',key_type=str)
        # create empty containers for spikes, lfps, etc.
        self.spks=D.Container() # dict contains spike timings
        self.resp=D.Container() # dict contains averaged response
        self.stim=D.Container() # dict contains auditory stimuli

        if cellid != None:
            self.load(cellid)
        
    def load(self,cellid=None):
        '''
        loads stimulus timing data

        Input:
          id_spk - an id that identifies the spike recording
          
        '''
        import scipy.io as io
        fn_stim = os.path.join(self.i.rec_path,self.i.rec_filenames[cellid])
        mat = io.loadmat(fn_stim)

        info = D.Info()
        for key in ['__header__','__version__','cellid']:
            info[key] = mat[key]
        self.i.update(info)        
        self.vstim = D.TimeSeries(mat['vstim'],32*10**3,info=info)
        
        i = 0
        for trial in mat['est_data']:
            self.stim[i] =  D.TimeSeries(trial.stim,32*10**3,info=info)
            self.spks[i] =  D.TimeSeries(trial.resp_raw.T,10**3,info=info)
            self.resp[i] =  D.TimeSeries(trial.resp,10**3,info=info)
            i+=1

