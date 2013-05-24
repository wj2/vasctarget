"""
This module contains customized classes and methods
derived from data.py and datalib.py.
"""
import numpy as N
import data as D
import utils as U
import datalib as DL
import plotlib as PL
import dsplib as DSP
import os,sys
from scipy import io

try:
    import tables as T
except:
    print "could not import pytables"

try:
    from neuropy import datapath_hirschlab as datapath
except:
    from neuropy import datapath
from neuropy import cachepath

__all__ = ['WholeCellRecording','combine_recordings']

class WholeCellRecording(D.Recording):
    """
    Container class for whole-cell recordings.

    Inputs:
      cellid   - a string specifying the cell, e.g. 'blabla0713_3_B_data'
    """
    def __init__(self,cell_id,rec_id,**kargs):
        defaultargs = dict(load_tic=True,load_tmg=True,load_sig=False,split=True,
                           neg=False,version=0,cc=True)
        defaultargs.update(kargs)

        self.i = D.Info(defaultargs)
        self.i.rec_id = rec_id
        self.i.cell_id = cell_id
        self.i.rec_path = os.path.join(datapath,cell_id)
        try:
            self.i.update(load_m(self))
        except:
            print 'coul not load m-file'
        if self.i.load_tic:
            try:
                load_tic(self)
            except:
                print "warning: could not stimulus time stamps"
        if self.i.load_tmg:
            try:
                self.load_spk()
            except:
                print "warning: could not load spk timings"
            if self.i.cc:
                try:
                    self.load_psp()
                except:
                    print "warning: could not load epsp timings"
            else:
                try:
                    self.load_psc()
                except:
                    print "warning: could not load epsc timings"
        if self.i.load_sig:
            try:
                self.load_sig()
            except:
                print "warning: could not load signal"
        if self.i.split:
            try:
                self.split_trials()
            except:
                print "warning: could not split trials"

    def cellfile(self,ext):
        ''' returns filename based on recording name '''
        fn = os.path.join(self.i.rec_path,self.i.rec_id+'.'+ext)
        return fn

    def get_trials_epochs(self):
        self.i.start_times = N.zeros(self.i.repeats,'d')
        # for i in range(self.i.repeats):
        #     self.i.start_times[i] = i* self.i.t_trial +(
        #         self.i.t1-self.i.lost_time)
        self.i.t0 = self.i.lost_time-self.i.t1 # t0 is the time at which the recording started
        for i in range(self.i.repeats):
            self.i.start_times[i] = i* self.i.t_trial
        if self.i.year == 2003:
            self.i.start_times[1:] -= self.i.t_start
        # self.i.epochs = DL.epochs_from_onsets(self.i.start_times,self.i.t_trial)
        self.i.epochs = DL.epochs_from_onsets(self.i.start_times,self.i.t_trial-self.i.t0) #shorten trials by offset lost at start

    def split_trials(self):
        self.i.split = 1
        self.get_trials_epochs()
        self.spks = self.spk[self.i.epochs]
        self.spks.i.t_max = self.i.t_trial
        if self.i.cc:
            self.psps = self.psp[self.i.epochs]
            self.psps.i.t_max = self.i.t_trial
        else:
            self.pscs = self.psc[self.i.epochs]
            self.pscs.i.t_max = self.i.t_trial
        if hasattr(self,'sig'):
            self.sigs = self.sig[self.i.epochs]
            # self.sig = sig_split_trials(self,self.sig[::10])

    def load_evt(self,name='evt',version=None):
        if version != None:
            fn_txt = self.cellfile('%s.%d.txt'%(name,version))
        else:
            fn_txt = self.cellfile('%s.txt'%name)
        evt = load_txt(fn_txt,load_t1=True,sorted=True)
        print ('%d '+name+' events') % len(evt)
        return evt

    def save_evt(self,timings,name='evt',version=None):
        if version != None:
            fn_txt = self.cellfile('%s.%d.txt'%(name,version))
        else:
            fn_txt = self.cellfile('%s.txt'%name)
        # t0 = self.i.lost_time-timings.i.t1
        # save_txt(fn_txt,timings-t0,t1=timings.i.t1)
        save_txt(fn_txt,timings,t1=timings.i.t1)

    def load_spk(self,**kargs):
        args = dict(version=self.i.version)
        args.update(kargs)
        spk = self.load_evt(name='spk',**args)
        # t0 = self.i.lost_time-spk.i.t1
        # self.spk = spk + t0
        self.spk = spk
        self.i.t1 = spk.i.t1
        return spk

    def load_psp(self,**kargs):
        args = dict(version=self.i.version)
        args.update(kargs)
        psp = self.load_evt(name='psp',**args)
        self.psp = psp
        self.i.t1 = psp.i.t1
        return psp

    def load_psc(self,**kargs):
        args = dict(version=self.i.version)
        args.update(kargs)
        psc = self.load_evt(name='psc',**args)
        self.psc = psc
        self.i.t1 = psc.i.t1
        return psc

    def load_sig(self):
        fn_rcd = self.cellfile('rcd')
        sig = load_rcd(fn_rcd)
        # sig.i.t0 = self.i.lost_time-sig.i.t1
        sig.i.t0 = 0.
        self.sig = sig
        self.i.t1 = sig.i.t1
        return sig

    def load(self):
        ''' load recording data in hdf5 format '''
        # try:
        if 1:
            fn_h5 = self.cellfile('%i.h5'%self.i.version)
            name = self.i.rec_path+'-'+self.i.rec_id
            fil = T.openFile(fn_h5,mode="r")
            params = fil.getNode(fil.root,"params")
            data = fil.getNode(fil.root,"data")
            for array in fil.listNodes(params,classname="Array"):
                # print "...loading %s"%array._v_name
                var = array._v_name
                self.i[var]=array.read()
            for array in fil.listNodes(data,classname="Array"):
                # print "...loading %s"%array._v_name
                var = array._v_name
                setattr(self.i,var,array.read())
            fil.close()
        # except:
        #     print "could not load cell data"

    def save(cell,all=0):
        ''' save recording data in hdf5 format '''
        raise NotImplementedError
        if cell.cellname == 'combined':
            all=1
        fn_h5 = cellfile(cell,np.str(cell.version)+'.h5')
        name = cell.cellpath+'-'+cell.cellname
        fil = T.openFile(fn_h5,mode="w",title=name)
        params = fil.createGroup(fil.root,"params","cell params")
        for var in cell.params.keys():
            fil.createArray(params,var,cell.params[var])
        data = fil.createGroup(fil.root,"data","cell data")
        save_vars = ['f','sf','skappa','pkappa','smu','pmu','spk_num','psp_num',
                     'inf','inft','repeats',
                     'spkmean','pspmean','f_cut','info','d_info','p_info']
        if all:
            save_vars.extend(['spktrain','psptrain','asigs','signals'])
        for var in save_vars:
            if cell.__dict__.has_key(var):
                if cell.__dict__[var] != None:
                    if np.isscalar(cell.__dict__[var]):
                        fil.createArray(data,var,np.float(cell.__dict__[var]))
                    else:
                        fil.createArray(data,var,cell.__dict__[var])
        fil.close()

    def get_asig(self,f,sf,neg=None,win_type='z',samplingrate=10**3):
        if neg == None: neg = self.i.neg
        self.i.f = f
        self.i.sf = sf
        # kilian: check if t0 in the following line is needed
        nbins = int(samplingrate*max(self.psp.max(),self.spk.max(),self.i.epochs['tstop'].max()))+1
        self.asig = DSP.analytic_signal(f,sf,self.psp.asraster(samplingrate,t0=self.i.start_times[0],nbins=nbins),
                                        neg=neg,win_type=win_type)
        if self.i.split:
            self.asigs = self.asig[self.i.epochs]

    def get_phases(self,psp=True,spk=True):
        if psp:
            self.psp.phase = N.angle(self.asig[self.psp.asidx(self.asig.i.samplingrate)])
            self.psp.power = N.abs(self.asig[self.psp.asidx(self.asig.i.samplingrate)])**2
        if spk:
            self.spk.phase = N.angle(self.asig[self.spk.asidx(self.asig.i.samplingrate)])
            self.spk.power = N.abs(self.asig[self.spk.asidx(self.asig.i.samplingrate)])**2
        if self.i.split:
            if psp:
                self.psps.phase = N.angle(self.asigs[self.psps.asidx(self.asigs.i.samplingrate)])
                self.psps.power = N.abs(self.asigs[self.psps.asidx(self.asigs.i.samplingrate)])**2
            if spk:
                self.spks.phase = N.angle(self.asigs[self.spks.asidx(self.asigs.i.samplingrate)])
                self.spks.power = N.abs(self.asigs[self.spks.asidx(self.asigs.i.samplingrate)])**2


def load_txt(filename,load_t1=None,sorted=False,data_type=float):
    ''' load event timings from txt file '''
    # return timings in [s]
    fil = open(filename,'r')
    if load_t1:
        t1 = float(fil.readline())
        if fil.readline() != '\n': assert 0, 'format error'
    else:
        t1 = 0
    # timings = D.Events(map(data_type,fil.readlines()))
    timings = D.Events(map(data_type,fil.readlines()),dtype=data_type)
    if sorted:
        timings.sort()
    fil.close()
    timings.i.t1 = t1
    return timings

def save_txt(filename,timings,t1=0):
    ''' save event timings in txt file '''
    fil = open(filename,'w')
    if t1:
        fil.write(N.str(t1)+'\n\n')
    for t in timings:
        fil.write(N.str(t)+'\n')
    fil.close()

def load_rcd(filename):
    '''
    save intracellular signal
    hirschlab binary format
    '''
    from scipy import io
    byteswap = (sys.byteorder=='big')
    fil = open(filename,'r')
    sig_size = io.fread(fil,1,'i','d',byteswap)
    sig = io.fread(fil,sig_size,'h','d',byteswap)/ 32768*5;
    sig = D.TimeSeries(sig,samplingrate=10**4)
    sig.i.t1 = float(io.fread(fil,1,'f','d',byteswap))
    fil.close()
    if len(sig)!=sig_size: assert 0, 'wrong file size'
    return sig

def save_rcd(filename,sig):
    '''
    save intracellular signal
    hirschlab binary format
    '''
    from scipy import io
    fil = open(filename,'w')
    sig_size = len(sig)
    byteswap = (sys.byteorder=='big')
    io.fwrite(fil,1,N.array([sig_size]),'i',byteswap)
    io.fwrite(fil,sig_size,sig*32768/5,'h',byteswap)
    io.fwrite(fil,1,N.array([sig.i.t1]),'f',byteswap)
    fil.close()

def load_m(rec):
    import re
    params = D.Info()
    fn = os.path.join(rec.i.rec_path,rec.i.rec_id+'.m')
    matchstr = re.compile(r"""\s*(\w+)\s*=\s*([\d\.-]+)\s*;""",re.VERBOSE)
    fil = open(fn,'r')
    for key,value in matchstr.findall(fil.read()):
        if int(float(value))==float(value):
            params[key] = int(value)
        else:
            params[key] = float(value)
    params.is_movie = ('startfilen' in params.keys())
    params.is_discs = ('out_diameter' in params.keys())
    params.year = int('20'+rec.i.cell_id[:2])
    if params.is_movie:
        params.in_terms = params.endfilen-params.startfilen+1
        params.out_terms = params.preterm
        params.f_rate = get_framerate(params.refreshrate,params.frmperterm)
        params.lost_time = float(1./params.f_rate)
        params.trial_terms = params.in_terms+params.out_terms
        params.start_terms = params.out_terms
        params.t_trial = float(params.trial_terms)/params.f_rate
        params.t_start = float(params.start_terms)/params.f_rate
    if params.is_discs:
        params.in_terms = params.sigtime
        params.out_terms = params.pretime
        params.lost_time = float(params.out_terms)
        params.t_trial = float(params.in_terms+2*params.out_terms)
        params.t_start = float(params.out_terms)
    return params

def load_tic(rec):
    fn =rec.cellfile('.tic.txt')
    if not os.access(fn,os.F_OK):
        print 'no stimulus time stamps available...'
        return
    tics = load_txt(fn)
    print 'old t_trial: %8.6f, old t_start: %8.6f' %(rec.i.t_trial,rec.i.t_start)
    if rec.i.is_movie:
        t_term = (tics[-1]-tics[0])/(len(tics)-1)
        rec.i.t_trial = rec.i.trial_terms*t_term
        rec.i.t_start = rec.i.start_terms*t_term
        rec.i.t_term = t_term
        rec.i.f_rate = 1/t_term
        rec.i.lost_time = t_term
    print 'new t_trial: %8.6f, new t_start: %8.6f' %(
        rec.i.t_trial,rec.i.t_start)

def get_framerate(refreshrate,frmperterm):
    if ((refreshrate == 133) & (frmperterm == 7)):
        return 19.00113
    if ((refreshrate == 144) & (frmperterm == 3)):
        return 48.05727
    elif ((refreshrate == 150) & (frmperterm == 5)):
        return 29.95780
    elif ((refreshrate == 150) & (frmperterm == 3)):
        return 49.92972
    elif ((refreshrate == 160) & (frmperterm == 5)):
        return 32.017395
    elif ((refreshrate == 160) & (frmperterm == 4)):
        return 40.021758
    else:
        assert 0, 'bad termrate'


def sig_split_trials(cell,data,samplingrate=0):
    raise NotImplementedError
    if not samplingrate:
        samplingrate = cell.i.samplingrate
        nt = int(cell.t_trial*cell.i.samplingrate)+2
    sig = zeros((cell.repeats,nt),dtype=data.dtype)
    # first start-time is negative
    lost_data= int(-cell.start_times[0]*cell.i.samplingrate+.5)
    if cell.repeats>1:
        sig [0,lost_data:nt] = data[:nt-lost_data]
        for i in range(1,cell.repeats):
            start_time = int(cell.start_times[i]*cell.i.samplingrate+.5)
            nmax = min(nt,len(data)-start_time)
            sig[i,:nmax] = data[start_time:start_time+nmax]
    else:
        assert 0, 'only 1 trial'
    return s.Signal(sig,samplingrate=cell.samplingrate)

def combine_recordings(rec1,rec2):
    rec2.spks.trials += rec1.i.repeats
    rec2.psps.trials += rec1.i.repeats
    rec1.spks = rec1.spks.concatenate(rec2.spks)
    rec1.spks.i.ntrials += rec2.i.repeats
    rec1.psps = rec1.psps.concatenate(rec2.psps)
    rec1.psps.i.ntrials += rec2.i.repeats
    if hasattr(rec1,'asigs'):
        rec1.asigs = rec1.asigs.concatenate(rec2.asigs)
    if hasattr(rec1,'sigs'):
        rec1.sigs = rec1.sigs.concatenate(rec2.sigs)
    rec1.i.repeats += rec2.i.repeats
    rec1.i.rec_id = 'combined'
    return rec1



