"""
Library for functions related to spike train analysis
"""

import numpy as N
import pylab as P
import dsplib as DSP
import data as D
import datalib as DL
import plotlib as PL

__all__ = ['isi', 'acor', 'xcor', 'psth', 'jpsth', 'find_events', 'sta',
           'fano', 'jitter', 'tuning', 'ssrevcor', 'dejitter_spks','oscore','cscore']

def get_best_shift(data1,data2,nmax):
    """
    determines by how many bins data2 has to be shifted in order to produce the best overlap.
    if data2 has trials, an array is returned.
    """
    if data2.has_trials():
        return N.array([get_best_shift(data1,data,maxshift) for data in data2],'i')
    i0 = (len(data1)+len(data2)-1)/2
    nmax = min(nmax,i0)
    c = N.correlate(data1,data2,'full')[i0-nmax:i0+nmax+1]
    return c.argmax()-nmax

def dejitter_spks(spks,maxshift=.02,samplingrate=1000,maxiter=10,spkwidth=.001):
    """
    dejitter spikes a la Optican and Richmond [ref]

    returns latencies and shifted spikes.
    maxshift is the maximal time to dejitter in s.
    spkwidth is the width of the smoothing used to compute the mean rate.
    maxiter is the maximal number of iterations to find best shift.
    """
    ntrials = spks.has_trials()
    shiftspks = spks.copy()
    shifts = N.zeros(ntrials,'i')
    nmax = int(maxshift*samplingrate)
    deltashifts = N.ones(ntrials,'i')
    while (deltashifts!=0).any():
        maxiter -= 1
        assert maxiter, 'too many iterations'
        shiftraster = shiftspks.asraster(samplingrate,t0=-maxshift)
        shiftrate = shiftspks.asrate(samplingrate,width=spkwidth*samplingrate,t0=-maxshift).mean(axis=0)
        shifttime = shiftrate.time()
        oldshifts = shifts.copy()
        shifts = get_best_shift(shiftrate,raster,nmax=nmax)
        deltashifts = shifts-oldshifts
        shiftspks.shift(deltashifts.astype('d')/samplingrate)
    latencies = -shifts.astype('d')/samplingrate
    return latencies,shiftspks

def isi_intervals(events):
    """
    Returns interspike-intervals (ISI).
    """
    # handle 1D or multi-trial events
    if events.has_trials():
        trlist = []
        for tr in events:
            trlist.extend(N.diff(tr))
        intervals = N.array(trlist,'d')
    else:
        intervals = N.diff(events)
    return intervals

def isi(events, xmax=0.1, nbins=100):
    """
    Returns an interspike-interval (ISI) distribution.

    Optional arguments: 'xmax' (s), 'nbins'
    Defaults: xmax = 100ms, nbins = 100.
    """
    intervals = isi_intervals(events)
    t = P.linspace(0,xmax,nbins+1)
    values,edges = N.histogram(intervals,t,new=False)
    #TODO: adapt function to use new histogram conventions

    return values[:-1], edges[:-1] # dump spurious bin > xmax

def acor(events,**kargs):
    """
    Auto-correlogram of events

    Optional arguments: 'xmax' (s), 'nbins'
    Defaults: xmax = 100ms, nbins = 200.
    """
    return xcor(events,events,**kargs)    

def xcor(events1,events2,xmax=.1,nbins=200,shift=False):
    """
    Cross-correlogram of events1 with events2
    ie. events from events1 occur at t=0 on the xcorrelogram.

    Optional arguments: 'xmax' (s), 'nbins', 'shiftpred',
    Defaults: xmax = 100ms, nbins = 200, shift = False.
    """
    time = P.linspace(-xmax,xmax,nbins+1)
    counts = N.zeros(nbins+1,'i')
    nt = events1.has_trials()
    if nt:
        events2list = [ev for ev in events2]
        for i,tr in enumerate(events1):
            for t in tr:
                 values,edges = N.histogram(events2list[i][N.where(abs(events2list[i]-t)<xmax)]-t,time,new=False)
                 #TODO: adapt function to use new histogram conventions
                 counts += values
                 if shift:
                    i2 = N.mod(i+nt/2,nt) # shift predictor by half of trials
                    v_shift,edges = N.histogram(events2list[i2][N.where(abs(events2list[i2]-t)<xmax)]-t,time,new=False)
                    counts -= v_shift
    else:
        if shift:
            raise ('cannot compute shift predictor without trials') 
        else:
            for t in events1:
                values,edges = N.histogram(events2[N.where(abs(events2-t)<xmax)]-t,time,new=False)
                counts += values
    return counts, edges # dump spurious bin > xmax (not necessary with numpy.histogram 1.3+ )

def psth(events, duration = None, binwidth = None, nbins = None, meanrate = True, smooth = True):
    """
    Generates a post-stimulus time histogram (PSTH).

    Passed spike times cut into trials, returns the population of each bin 'y'
    and the left-hand edges of each time bin 'bins'. If specified, 'nbins' takes
    precedence over 'binwidth'. The 'meanrate' flag determines if the mean spike
    rate is returned, otherwise it's just the spike count in each bin.
    If 'smooth' is True, the histogram values are gaussian smoothed.

    Optional arguments:  'binwidth' (s), 'nbins', 'meanrate', 'smooth'.
    Defaults: binwidth = 1ms, nbins = None, meanrate = True, smooth = True.
    """
    
    # histogram arguments
    if duration == None:
        duration = N.ceil(events.max())
    if binwidth == None:
        binwidth = 0.001 # 1ms
    if nbins == None:
        nbins = int(N.ceil(duration/binwidth))
    else: duration = nbins*binwidth
    
    # build histogram
    bins = N.arange(0,duration+binwidth,binwidth)
    #bins = bins[:nbins] # strip extraneous bin in cases where: duration mod binwidth ~ 0
    nt = events.has_trials() # get number of trials

    y = N.zeros(nbins) # need to initialize histogram as float or smoothg throws an error
    y+= N.histogram(events,bins)[0]
    print N.shape(y)

    if meanrate and nt: y/=float(nt)*binwidth
    if smooth:
        gstd = 0.01/binwidth
        winsz = 6*gstd
        y=DSP.smoothg(y,winsz) # to be changed when smoothg takes seconds, not samples

    return y, bins

def jpsth(n1_times, n2_times, nbins): # return the jpsth matrix
    """
    Generate a Joint Peri/Post-Stimulus Time Histogram

    inputs:
    n1_times and n2_times should contain spike times cut into trials
    nbins is the number of histogram bins to use

    outputs:
    returns a dictionary containing entries:
    'jpsth_raw': the raw JPSTH
    'jpsth_poisson': the JPSTH predicted from the outer product of 1-dimensional histograms
    'n1_times': bin labels for neuron 1
    'n2_times': bin labels for neuron 2
    'n1_hist': 1-dimensional histogram (PSTH) for neuron 1
    'n2_hist': 1-dimensional histogram (PSTH) for neuron 2
    """
    
    # n1_times and n2_times - 1 dimensional arrays holding times for all the spikes of interest.
    # nbins - number of bins along each axis
    
#    from jpsth import *

    # fill in the times for the binning boundaries
    min_time = min( n1_times.min(), n2_times.min() )
    max_time = max( n1_times.max(), n2_times.max() )

    d1_hist_bounds = N.arange( nbins + 1 ).astype(float)
    d1_hist_bounds = d1_hist_bounds * (max_time - min_time) / (d1_hist_bounds.max() - d1_hist_bounds.min())
    d1_hist_bounds = d1_hist_bounds + min_time
#    print min_time, max_time, d1_hist_bounds

    d1_hist_indexes = N.arange( nbins )
    n1_hist_bounds = N.empty( [nbins,nbins,2], float )
    n2_hist_bounds = N.empty( [nbins,nbins,2], float )
    jpsth_matrix = N.zeros( [nbins, nbins] )
    n1n2_outer = N.zeros( [nbins, nbins] )

    # populate the bounds in the full size bounds arrays
    for ii in d1_hist_indexes:
        n1_hist_bounds[:,ii,0] = d1_hist_bounds[ii]
        n1_hist_bounds[:,ii,1] = d1_hist_bounds[ii+1]
        n2_hist_bounds[ii,:,0] = d1_hist_bounds[ii]
        n2_hist_bounds[ii,:,1] = d1_hist_bounds[ii+1]

    n1_hist = N.zeros( [nbins] )
    n2_hist = N.zeros( [nbins] )

    # now for each trial . . . .
    ntrials = len(n1_times)
    for trial_i in xrange( ntrials ):
        # ******* this bit should probably be redone with the builting 2d histogramming function
        # find 1d indexes for each of the neurons
        n1_ins = N.empty( [len(n1_times[trial_i])] )
        for i in xrange(len(n1_times[trial_i])):
            n1_ins[i] =  N.max( d1_hist_indexes[d1_hist_bounds[0:-1]<=n1_times[trial_i][i]] )

        n2_ins = N.empty( [len(n2_times[trial_i])] )
        for i in xrange(len(n2_times[trial_i])):
#           print i, n2_times[trial_i][i], d1_hist_indexes[d1_hist_bounds[0:-1]<=n2_times[trial_i][i]]
            n2_ins[i] =  N.max( d1_hist_indexes[d1_hist_bounds[0:-1]<=n2_times[trial_i][i]] )

        # step through all the input times
        for n1_ind in n1_ins:
            for n2_ind in n2_ins: psth_matrix[n1_ind,n2_ind] = jpsth_matrix[n1_ind,n2_ind] + 1

        # calculate 1 dimensional histograms
        for n1_ind in n1_ins:
            n1_hist[n1_ind] = n1_hist[n1_ind] + 1
        for n2_ind in n2_ins:
            n2_hist[n2_ind] = n2_hist[n2_ind] + 1

    # and now the outer-product poisson case . . .
    for n1_ind in xrange( nbins ):
        for n2_ind in xrange( nbins ):
            n1n2_outer[n1_ind, n2_ind] = n1_hist[n1_ind] * n2_hist[n2_ind]
    n1n2_outer = n1n2_outer / ntrials

    # the output will be a dictionary holding the JPSTH matrix and the start and stop time bins    
    results = {'jpsth_raw': jpsth_matrix, 'n1_times': n1_hist_bounds, 'n2_times': n2_hist_bounds, 'n1_hist':n1_hist, 'n2_hist':n2_hist, 'jpsth_poisson':n1n2_outer }

    return results


def find_events(psth,edges,threshold=None):
    """
    Finds events in a PSTH. 

    Given a 'psth', returns 'peaks' and 'epochs' corresponding to 
    maximal firing events and the boundaries of these events. Event
    boundaries are defined as significant (p < 0.05) minima between
    adjacent peaks as described in several Warland and Reinegal
    papers [refs needed]. Minima are defined relative to the
    geometric mean of adjacent maxima, that is: sqrt(Mi*M(i+1)) <= 3
    
    TODO: remove hard-coded hacks that assume 1000sps PSTH!
    """
    # find PSTH peaks
    maxidx = DSP.arglocmax(psth)
    pkval  = psth[maxidx]
    pktime = edges[maxidx]

    # find PSTH epochs
    mint=[0]
    minv=[0]
    peaksv=[]
    peakst=[]
    for j in xrange(len(maxidx)-1):
        thresh = N.sqrt(pkval[j]*pkval[j+1])/3
        mi = psth[maxidx[j]:maxidx[j+1]].argmin()+maxidx[j]
        m  = psth[mi]
        if m <= thresh:
            minv.append(m)
            mint.append(edges[mi])
            peaksv.append(pkval[j])
            peakst.append(pktime[j])
            
    # append last event of trial ~ hack!
    peaksv.append(pkval[-1])
    peakst.append(pktime[-1])
    idx = int(edges[-1]*1000)
    while psth[idx] == 0:
        idx-=1
    mint.append(edges[idx])
    
    # convert consecutive minima into event epochs
    epochs = []
    for i in xrange(len(mint)-1):
        # next four lines adjust epoch start times to onset of PSTH peak
        idx = int(mint[i]*1000)
        while psth[idx] == 0:
            idx+=1
        mint[i]=idx/1000.
        epochs += [(mint[i], mint[i+1])]

    peaks  = D.Events(peakst)
    peaks.v = N.array(peaksv)
    return peaks, D.Epochs(epochs, 'Peak')


def sta(spks, stim, movie, tslices = 10, dim = None, roi = None, spkrate = True):
    """
    Spike-triggered average based spatiotemporal receptive field (STRF).
 
    Uses spike-triggered average of responses to a movie stimulus (typically
    an m-sequence, white noise, or sparse bar movie) to generate the linear
    receptive field estimate of a neuron at a fixed number of 'tslices', or
    at arbitrary times if a vector is passed to 'tslices'.

    TO DO: - would be better if this function was independent of metadata in stim object
           
    Optional arguments: tslices, list of desired timepoints (frame/refresh intervals),
    roi=[xmin, xmax, ymin, ymax], crop to region of interest (in pixels), normalize (spikes/s).
    """
    # strf arguments
    # if dim == None: dlen = 1
    # else: d = stim.i[dim]; dlen = len(d)
    
    # blank frames of movie if they don't match specified stimulus dimension value
    if dim != None: movie[N.unique(stim[N.where(stim['brightness']==0)].v)] = 0
    
    if N.isscalar(tslices): nbins = tslices
    else: nbins = len(tslices)

    # crop strf if roi specified
    if roi != None : movie = DL.crop(movie,roi)
    sta = N.zeros((nbins, movie.shape[1], movie.shape[2]), 'd')
    #if dlen == 1:
    
    # computes sta in two (equivalent) ways depending on whether
    # the number of frames or specific time slices are specified
    if N.isscalar(tslices):
        frbins = stim[N.nonzero(N.diff(stim.v))[0]+1] # movie frame bins
        frbins = N.concatenate(([stim[0]], frbins, [stim[-1]+1/stim.i.refreshrate])) # first and last frame bins
        binspks, blah = N.histogram(spks, frbins, new=True) # bin spikes to frames
        binspks[0:nbins] = 0; # discard pre-stimulus spikes
        binspks = binspks[:-1] # -1 is a temporary hack to deal with spurious frame in mseq rec.stim.din
        for i in N.nonzero(binspks)[0]:
             sta += movie[i-nbins+1:i+1]*binspks[i]
    else: # use this method is for tslices that are discontiguous (about 50% slower)             
        binspks = spks.evt2idx(stim) # movie frame bins
        binspks = binspks[binspks < len(stim.v)-1] # -1 is a temporary hack to deal with spurious frame in mseq rec.stim.din
        binspks = binspks[binspks > max(tslices)] # discard pre-stimulus spikes
        for j, t in enumerate(tslices):
            k = nbins-1-j
            for i in binspks:
                sta[k] += movie[stim.v[i-t]]    
    
    # else: # multi-dimensional stimulus, e.g. sparse noise with both white and black spots
    #     binspks = spks.evt2idx(stim)# index spike times to movie frames
    #     binspks = binspks[binspks < len(stim.v)-1] #-1 is a temporary hack to deal with spurious frame in mseq rec.stim.din
    #     binspks = binspks[binspks > max(tslices)] # discard pre-stimulus spikes
    #     for idx, stimdim in enumerate(d):   # with ndenumerate, it would be possible to...
    #         for j, t in enumerate(tslices): # remove tslice loop and add extra dimensions
    #              condmask = stim[dim][binspks-t] == stimdim
    #              for i in binspks[condmask]:
    #                  sta[idx,j] += movie[stim.v[i-t]]
    #sta = sta.squeeze() # collapse first dimension when dim not specified
    if spkrate: sta = sta/float(stim[-1]-stim[0]) # spikes/s, more meaningful to have integrated area = 1?
    
    return D.Image(sta)


def fano(events):
    """
    Computes the fano factor for a set of events with multiple trials.
    """
    if events.has_trials():
        x = events.events_per_trial()
        return N.var(x)/N.mean(x)
    else:
        assert 0, 'fano factor undefined for single trials'
    

def jitter(events, ms = True):
    """
    Computes the standard deviation of first spike times
    in a spike event with multiple trials.

    If ms is true, jitter values are return in milliseconds.
    """
    if events.has_trials():
        if ms : multiplier = 1000
        else : multiplier = 1
        return N.std([sp[0]*multiplier for sp in events.iter_dim()])
    else:
        assert 0, 'spike time jitter is undefined for single trials'


def tuning(stim, events, dim, toff=0, durn=None, cond_mask=N.array(True), metrics=True):
    """
    Generates a neuronal tuning curve.

    Generates a neuronal tuning curve for parametric stimuli based
    on firing rates. Expects a NP.Stimulus object,'stim', that has
    been processed with BuildSweepTable, spike times in 'events', and
    a list of the stimulus dimensions, 'dim', to evaluate combinatorially.

    An error is raised if no stimuli match user-specified mask conditions.

    Optional arguments: 'toff', in s, shifts counting window (wrt stimulus onsets).
                        'durn', in s, specifies counting window (from stim onsets).
                        'mask', boolean array restricts evaluation to these conditions.
                        'metrics', return trial-averaged event count & standard deviation.
    TODO: would be much faster to boolean slice with stim.v rather than i.dim, and,
          in addition, blank stimuli cannot be analysed here.
    """
    dlen = tuple([len(stim.i[d]) for d in dim])
    if N.ndim(cond_mask) > 1 : cond_mask = N.logical_and.reduce(cond_mask)
    assert cond_mask.any(),'No stimulus with those parameters exists'   
    frms_per_sweep = round(stim.i.sweeptimeSec*stim.i.refreshrate) #assumes identical sweep lengths
    ncond = N.multiply.reduce(dlen) 
    avgcount = N.zeros(ncond, float)
    std_dev = N.zeros(ncond, float)
    cspks, i = None, 0
    for idx in N.ndindex(dlen): # better to use ndenumerate? native numpy, faster?
        mask = [stim[d]==stim.i[d][idx[j]] for j,d in enumerate(dim)] # mapping to v then making a single
        mask = N.logical_and.reduce(cond_mask & mask)                 # boolean slice of stim.v would be much faster!
        if mask.any() == False : continue # stimulus condition does not exist
        onsets  = stim[mask][::frms_per_sweep] + toff
        if durn == None:
            offsets = stim[mask][frms_per_sweep-1::frms_per_sweep]+1/stim.i.refreshrate
            epochs = D.Epochs(onsets,offsets,name='trials')
        else : epochs = DL.epochs_from_onsets(onsets, durn)
        spks = DL.cut_events(events, epochs)  
        if metrics : 
            spkcount = spks.events_per_trial()
            avgcount[i] = N.mean(spkcount)
            std_dev[i]  = N.std(spkcount)
        spks.cond = N.ones(len(spks),int)*i
        if cspks == None : cspks = spks
        else : cspks = cspks.concatenate(spks)
        i+=1
    if metrics : return avgcount, std_dev, cspks
    else: return cspks


def ssrevcor(spks, stim, tslices = range(11), normalize = False):
    """
    Sub-space reverse correlation.
    
    Use reverse correlation to estimate neuronal tuning curves of arbitrary dimensionality.

    Optional arguments: 'tslices', list of desired timepoints (in refreshes),
                        'normalize', returns spike count in Hz (spikes/s).
    """
    nslices = len(tslices)
    ncond = len(set(stim.v))
    spkcount = N.zeros((ncond,len(tslices)), float)
    stim.v[stim.v==65535] = ncond-1 # map blank stim to last v condition

    # remove any spikes that occured before or after the stimulus
    spks = spks[(spks > (stim[0]+nslices/stim.i.refreshrate)) & (spks < stim[-1])]

    # index spike times to stimulus frames
    idxspks = spks.evt2idx(stim)
    
    # compute ndim,t histogram
    for i in idxspks:
        for t in tslices:
            dim = stim.v[i-t]
            spkcount[dim,t]+=1 
    if normalize: spkcount/=stim.i.nrepeats*stim.i.sweeptimeSec # spikes/s 
    stim.v[stim.v==ncond-1] = 65535 # restore blank stim
    return spkcount


def get_best_shift(data1,data2,nmax):
    """
    determines by how many bins data2 has to be shifted in order to produce the best overlap.
    if data2 has trials, an array is returned.
    """
    if data2.has_trials():
        return N.array([get_best_shift(data1,data,nmax) for data in data2],'i')
    i0 = min(len(data1),len(data2))-1#(len(data1)+len(data2)-1)/2
    nmax = min(nmax,i0)
    c = N.correlate(data1,data2,'full')[i0-nmax:i0+nmax+1]
    best_shift = c.argmax()-nmax
    if c.argmin() == best_shift+nmax:
        return 0
    else:
        return best_shift


def dejitter_spks(spks,maxshift=.02,samplingrate=1000,maxiter=1,spkwidth=.001):
    """
    dejitter spikes a la Optican and Richmond [ref]

    returns shifts that best align trials (given a maximal shift of maxshift).
    maxshift is the maximal time to dejitter in s (+/- maxshift).
    spkwidth is the width of the smoothing used to compute the mean rate.
    maxiter is the maximal number of iterations to find best shift.
    """
    ntrials = spks.has_trials()
    shiftspks = spks.copy()
    shifts = N.zeros(ntrials,'i')
    nmax = int(N.abs(maxshift*samplingrate))
    
    raster = spks.asraster(samplingrate,t0=-2*maxshift)
    nbins = int((spks.max()+4*maxshift)*samplingrate+1)
    
    deltashifts = N.ones(ntrials,'i')
    while (deltashifts!=0).any():
        maxiter -= 1
        #assert maxiter, 'too many iterations'
        shiftrate = shiftspks.asrate(samplingrate,width=spkwidth*samplingrate,t0=-2*maxshift,nbins=nbins).mean(axis=0)
        shifttime = shiftrate.time()
        oldshifts = shifts.copy()
        shifts = get_best_shift(shiftrate,raster,nmax=nmax)
        deltashifts = shifts-oldshifts
        shiftspks.shift(deltashifts.astype('d')/samplingrate)
        print maxiter, (shifts*1000.0)/samplingrate
        if maxiter < 1 : break # one pass only
    return shifts.astype('d')/samplingrate


def oscore(events,fmin,fmax,samplingrate):
    '''
    oscillation score

    see Muresan RC, Jurjut OF, Moca VV, Singer W, Nikolic D.
    The Oscillation Score: An Efficient Method for Estimating Oscillation Strength in Neuronal Activity.
    J Neurophysiol. 2007 Dec 26;
    '''
    # fmin : lower bound for frequency band of interest
    # fmax : upper bound for frequency band of interest
    fc = samplingrate # correlogram frequency in Hz (1/binwidth[s])
    # w : width of one flank of the ACH 
    # W : width of one flank of the buffers = 2*w
    # sgm_fast  : standard deviation of fast Gaussian kernel (in bins)
    # S_fast : fast smoothed ACH 
    # sgm_slow  : standard deviation of slow Gaussian kernel (in bins)
    # S_slow : slow smoothed ACH

    # (0) iterate through trials
    if events.has_trials():
        fos=[]
        Os =[]
        for dat in events.iter_dim('trials'):
            fosi,Osi = oscore(dat,fmin,fmax,samplingrate)
            fos.append(fosi)
            Os.append(Osi)
        return fos, Os

    # (1) compute auto-correlation histogram (ACH)
    w = int(2**(N.floor(max(N.log2(3.*fc/fmin),N.log2(fc/4.)))+1))
    W = 2*w
    values,edges = acor(events,xmax=float(w)/fc,nbins=W)
    ACH = D.TimeSeries(values.astype('d')[:-1],samplingrate=samplingrate,t0=edges[1])
    i0 = w # index of time 0
    ACH_half = D.TimeSeries(ACH.copy()[i0:],samplingrate=samplingrate,t0=0)
    
    # (2) smoothing
    sgm_fast = min(2,134./(1.5*fmax))*fc/1000. # obviously in bins
    S_fast = DSP.smoothg(ACH,std=sgm_fast/fc)

    # (3) removing the central peak
    sgm_slow = 2*134./(1.5*fmin)*fc/1000. # obviously in bins
    S_slow = DSP.smoothg(ACH,std=sgm_slow/fc)
    slope = N.diff(S_slow)*W/S_slow[i0]
    for imin in N.arange(i0,0,-1):
        if slope[imin-1] <= N.tan(N.pi*10/180):
            break
    SP_fast = S_fast.copy()
    SP_fast[imin:-imin] = SP_fast[imin]

    # (4) applying FFT
    SFFT = N.abs(N.fft.fft(SP_fast))
    freqs = N.fft.fftfreq(len(SP_fast),d=1./samplingrate)

    # (5) estimation of oscillation score
    f_ind = ((freqs>fmin) & (freqs<fmax))
    i_max = N.abs(SFFT[f_ind]).argmax()
    fos = freqs[f_ind][i_max]
    Mpeak = N.abs(SFFT[f_ind][i_max])
    Mavg = N.abs(SFFT).mean()
    Os = Mpeak/Mavg
    return fos,Os

def cscore(events):
    '''
    confidence score
    
    see Muresan RC, Jurjut OF, Moca VV, Singer W, Nikolic D.
    The Oscillation Score: An Efficient Method for Estimating Oscillation Strength in Neuronal Activity.
    J Neurophysiol. 2007 Dec 26;
    '''
    dat = N.array(events)
    cv = dat.std()/dat.mean()
    return 1/(1+cv)

def xc(neuron1,neuron2,stims,varargin):
    """
    [n,t] = xc(neuron1,neuron2,stims,varargin);
    Cross correlates spikes from neuron1 with spike from neuron2,
    ie. spikes from neuron1 occur at t=0 on the xcorrelogram. Expects
    neuron data structures (see loadneuron()).

    If you're building up the xc over trials and/or you want to do
    shuffle correction, each stim's sweepranges will be used. Make
    sure the sweeps in sweepranges are meaningful. If not (such as for
    a movie stim, where the sweep indices are just frame indices),
    regroup them into something meaningful using regroupsweepis().

    Optional args: range, res, units, newfigure, sweepis, shufflecorrect.
    """
    pass


def xcs(neurons,stims,varargin):
    """
    xcs(neurons,stims,varargin);

    Cross correlates spikes from all neurons within a neuron structure array
    "neurons" (generated using loadneurons()) and plots 'em.

    Click a plot for a more detailed popup view: on the popup, use
    arrow keys to pan and zoom, +/- control bin size < and > page left
    and right, u changes units, r reverses the xcorrelogram 0
    recenters on t==0.

    If you're building up the xc over trials and/or you want to do
    shuffle correction, each stim's sweepranges will be used. Make
    sure the sweeps in sweepranges are meaningful. If not (such as for
    a movie stim, where the sweep indices are probably just  frame
    indices), regroup them into something meaningful using
    regroupsweepis().

    Optional args: range, res, units, minupperylimit, maxlowerylimit,
    sweepis, shufflecorrect.

    """
    pass
    
def stc(neuron,stim,rcmovie,varargin):
    """
    [rf,tms] = stc(neuron,stim,rcmovie,varargin);

    Uses spike-triggered covariance of reponses to a movie stim
    (generally an msequence movie in a cell array) to generate the
    receptive field of a neuron at the specified revcor times. Returns
    matrix rf, or plots the RFs as a series of images. Optionally
    plays them as a movie. Click an image and hit the spacebar to
    toggle between fully scaled colour (to accentuate image structure)
    and properly normalized colour.

    Optional args: t (list of desired revcor timepoints in ms), nrfs,
    scale (screen pix per movie pix), playrf, updateperiod, cmap,
    returnmatrix.
    """
    pass


def rate(neuron,varargins):
    """
    [r,t] = rate(neuron,varargins);

    Calculates and plots firing rate (in Hz) using a sliding window,
    or bins, or variable width bins. Passed neuron structure 'neuron'.
    Returns the firing rate 'r' at each timepoint 't'.

    Optional arguments: range, tres, winfn, sig, nspb, units, newfigure, ylock.

    Defaults:
    t0 = neuron.spikes(1).
    tend = neuron.spikes(end).
    range = [t0 tend]; %[t0 t0+3e6]; % interval in us over which to calculate rate.
    tres = 50e3; % time resolution to use, in us.
    winfn = 'vbin'; %; % gauss, rect, vbin, bin.
    sig = 200e3; % in us, width of 'rect', sigma of 'gauss'.
    nspb = 4; % number of spikes per bin, when using 'vbin' as the winfn.
    units = 's'; % display units for time axis.
    newfigure = 1; % 1==create a new figure, otherwise==just update existing one (no effect if expecting an output arg).
    ylock = 0; % sets the y limits if non-zero (pass a 2 column vector), otherwise, leaves matlab to autoscale the y axis.
    extrapolatedrate = 0; % when necessary, the extrapolated rate set when doing interpolation.
    """
    pass


def ratePdf(neuron,stims,varargins):
    """
    [p,edges] = ratePdf(neuron,stim(s),varargins);

    Returns a firing rate PDF of a neuron. 'p' are the probabilities, and
    'edges' are the boundaries of the corresponding spike rate bins.
    """
    pass


def rateJPdf(neuron1,neuron2,stims,varargin):
    """
    [N,entropy] = rateJPdf(neuron1,neuron2,stim(s),varargin);
    
    Returns the firing rate joint probability distribution function of
    2 neurons in matrix N, along with the entropy.

    """
    pass


def rateJPdfs(neurons,stims,varargin):
    """
    [N,entropy] = rateJPdfs(neurons,stim(s),varargin);

    Plots all combinations of ratejpdfs in neurons structure.
    """
    pass


def rateFPdf(neuron1,neuron2,stims,varargin):
    """
    [N,entropy] = rateFPdf(neuron1,neuron2,stim(s),varargin);

    Returns the theoretical firing rate factorial probability density
    function of 2 neurons in matrix N by taking the outer product of
    their respective firing rate PDFs.
    """
    pass


def rateFPdfs(neurons,stims,varargin):
    """
    [N,entropy] = rateFpdfs(neurons,stim(s),varargin);

    Plots all combinations of ratefpdfs in neurons structure.
    """
    pass


def isiRate(spikes,varargin):
    """
    [r,t] = isiRate(spikes,varargin);

    Calculates and plots firing rate according to inter spike interval.
    Passed spike times in a row vector 'spikes'.
    Returns the firing rate 'r' at each timepoint 't'.

    Optional arguments: timepoints vector 't'.
    """
    pass


def hist2(s1,s2,varargin):
    """
    [N,edges1,edges2] = hist2(s1,s2,varargin);

    Generates a 2D histogram. Passed equilength simultaneous signals
    in vectors s1 and s2. Each edges vector holds the bin edges of the
    respective signal. Returns 2D histogram in matrix N. Passing
    'edges' args overrides passing of 'nbins' args.

    Optional args: nbins1, nbins2, edges1, edges2.
    """
    pass

