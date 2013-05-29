"""
Library for functions related to plotting
"""

import numpy as np
import data,model,utils
import spikelib as spk
from neuropy import pylab as P

__all__ = ['multiplot','plot_raster','plot_timeseries','plot_isi','plot_xcor','plot_acor',
           'plot_phasedist','plot_psth','plot_sta','plot_avg_timeseries','plot_saccade_hist',
           'normed_hist']

def _get_next_cycle_color(figure=None):
    if figure==None:
        figure = P.gcf()
    ax = figure.gca()
    if color==None:
        color=ax._get_lines._get_next_cycle_color()
    return color

def multiplot(numRows,numCols,**moreargs):
    """
    Creates multiple axes (like subplot) and cycles through them

    When called the first time, it attaches itself to the current figure.
    """

    # default arguments
    args = dict(figure=P.gcf())
    args.update(moreargs)

    fig = args['figure']
    if not hasattr(fig,'_multidx'):
        fig._multidx=0
        # attach multiplot function to figure
        fig.next = lambda: multiplot(numRows,numCols,**moreargs)
    else:
        fig._multidx+=1
        # create subplot
        ax = P.subplot(numRows,numCols,fig._multidx)
        # set current axis
        fig.sca(ax)

def plot_raster(events,overplot=False,trange=None,**args):
    """
    Generates a raster plot of event times.

    y-axis offset is incremented automatically.
    to prevent offset incrementation, define overplot=True.

    Optional arguments: 'spacing', 'spikeheight','duration',
                        'showlines', 'units'.

    TO DO: - use figure-centric rather than data-centric offsets for row spacing
           - add "sub-group" capability for sub-grouping of trials
    """

    # default plot arguments
    plotargs = dict(ms=5,mew=1,linestyle='None',marker='|',figure=P.gcf())
    plotargs.update(args)

    # get current axis
    ax = plotargs['figure'].gca()

    nt = events.has_trials() # get number of trials
    t0 = events.i.get('t0',0)

    if not hasattr(ax,'_rowidx'):
        ax._rowidx = -1
    else:
        ax._rowidx -= nt+1

    if trange == None:
        y = ax._rowidx - events.trials
        P.plot(events+t0,y,**plotargs)
    else:
        epoch = data.Epochs(trange,align=False)
        y = ax._rowidx - (events+t0)[epoch].trials
        P.plot((events+t0)[epoch],y,**plotargs)
        P.xlim(trange)
    P.xlabel(events.i.get('xlabel','time (s)'))
    P.ylim([ax._rowidx-nt,0])
    P.yticks([])
    if overplot: ax._rowidx += nt+1


def plot_timeseries(timeseries,multiline=True,events=None,gain=None,**args):
    """
    Generates a multi-line plot of time series
    """

    # default plot arguments
    plotargs = dict(figure=P.gcf())
    plotargs.update(args)

    # get current axis
    ax = plotargs['figure'].gca()

    if not hasattr(ax,'_rowidx'): ax._rowidx=0
    nt = max(1,timeseries.has_trials()) # get number of trials
    if nt == 1: timeseries = np.reshape(timeseries,(1,len(timeseries))) # single trial
    if gain == None: gain = .25/timeseries[0].std()

    for tr in range(nt):
        ax._rowidx -= 1
        P.plot(timeseries.time(),gain*timeseries[tr]+ax._rowidx*multiline,**plotargs)
        if events!=None:
            color = ax.lines[-1].get_color()
            rasterargs = plotargs.copy()
            rasterargs.update(dict(ms=5,mew=2,linestyle='None',marker='|',color=color))
            # rasterargs.update(dict(ms=10,mew=2,linestyle='None',marker='|',color='k',alpha=1))
            evts = events.select_dim(tr)
            P.plot(evts,ax._rowidx*np.ones_like(evts),**rasterargs)
            # P.plot(evts,gain*timeseries[tr][evts]+ax._rowidx,**rasterargs)
    # P.xlabel(timeseries.xlabel())
    # P.ylabel(timeseries.ylabel())
    #P.ylim([0,nt+1])
    P.yticks([])


def plot_avg_timeseries(timeseries, show_error = True, SEM = False, **args):
    """
    Plots averaged time series +/- standard deviation or SEM
    """

    # default plot arguments
    plotargs = dict(figure=P.gcf(), color = 'k')
    plotargs.update(args)
    plotargs['color'] = args.get('facecolor',plotargs['color'])

    timeseries_avg = timeseries.mean(axis=0)
    t = timeseries.time()
    if show_error :
        timeseries_std = timeseries.std(axis=0)
        if SEM : 
            if timeseries.has_trials():
                timeseries_std /= np.sqrt(timeseries.has_trials())
            else:
                timeseries_std /= np.sqrt(timeseries.shape[0])
        upper_er = timeseries_avg + timeseries_std
        lower_er = timeseries_avg - timeseries_std
        x = np.concatenate((t, t[::-1]))
        y = np.concatenate((upper_er, lower_er[::-1]))
        P.fill(x, y, alpha = 0.2)
        #P.fill(x, y, facecolor = 'b', edgecolor = 'b', alpha = 0.2)
    P.plot(t, timeseries_avg, **plotargs)

    P.xlabel('time (s)')
    if timeseries.has_trials():
        P.ylabel('average (n=%s)'%timeseries.has_trials())
    else:
        P.ylabel('average (n=%s)'%timeseries.shape[0])

def plot_isi(events, xmax=0.1, nbins=100, fit_gamma=False, **args):
    """
    Plots an interspike-interval (ISI) distribution.

    Optional arguments: 'xmax' (s), 'nbins'
    Defaults: xmax = 100ms, nbins = 100.
    """
    # default plot arguments
    plotargs = dict(figure=P.gcf())
    plotargs.update(args)

    # get current axis
    fig = plotargs['figure']
    if hasattr(fig,'next'): # select next subplot
        fig.next()
    ax = fig.gca()

    # compute and plot histogram
    values,edges = spk.isi(events, xmax=xmax, nbins=nbins)
    P.bar(edges, values, width=xmax/nbins, **args)

    if fit_gamma:
        k,xmean = model.fit_gamma(spk.isi_intervals(events))
        db = edges[1]-edges[0]
        if (k<=1): edges = edges[edges>0]
        P.plot(edges,db*len(events)*model.gamma(edges,k,xmean),'k',linewidth=2)

    P.xlim([0,xmax])
    if hasattr(events.i,'units'):
        P.xlabel('interval (%s)'%events.i.units)
    else:
        P.xlabel('interval')
    return values,edges

def plot_acor(events,**kargs):
    counts,edges = plot_xcor(events,events,**kargs)
    centerpeak = counts.argmax()
    ymax = utils.round_up(counts[centerpeak+1:].max())
    P.ylim(0,ymax)
    return counts,edges

def plot_xcor(events1,events2,xmax=.1,nbins=200,shift=False,**kargs):
    counts,edges = spk.xcor(events1,events2,xmax=xmax,nbins=nbins,shift=shift)

    plotargs = dict(color='k')
    plotargs.update(kargs)
    dt = np.diff(edges)[0]
    P.bar(edges-dt/2,counts,dt,**plotargs)
    ymax = utils.round_up(counts.max())
    P.axis([-xmax,xmax,0,ymax])
    # P.title(title)
    if hasattr(events1.i,'units'):
        P.xlabel('time delay (%s)'%events1.i.units)
    else:
        P.xlabel('time delay')
    return counts,edges

def plot_phasedist(phases,bins=37,center_pi=False,**args):

    # default plot arguments
    plotargs = dict(figure=P.gcf(),color='b',alpha=.75)
    plotargs.update(args)

    # get current axis
    fig = plotargs['figure']
    if hasattr(fig,'next'): #select next subplot
        fig.next()
    ax = fig.gca()

    if center_pi:
        phases = np.mod(phases,2*np.pi)
        (nphases, bins, patches) = P.hist(phases,P.linspace(0,2*np.pi,bins+1),normed=1)
    else:
        phases = utils.smod(phases,2*np.pi)
        (nphases, bins, patches) = P.hist(phases,P.linspace(-np.pi,np.pi,bins+1),normed=1)
    P.setp(patches, 'facecolor', plotargs['color'], 'alpha', plotargs['alpha'])
    ymax = utils.round_up(nphases.max())
    (kappa,mu,p,n) = model.phasedist(phases)
    if center_pi:
        P.plot (P.linspace(0,2*np.pi,100),model.mises(P.linspace(0,2*np.pi,100),
                                                 kappa,mu),'k',linewidth=2)
        P.axis([0,2*np.pi,0,ymax])
    else:
        P.plot (P.linspace(-np.pi,np.pi,100),model.mises(P.linspace(-np.pi,np.pi,100),
                                                   kappa,mu),'k',linewidth=2)
        P.axis([-np.pi,np.pi,0,ymax])
    P.xlabel('phase (rad)')
    labels = P.getp (P.gca(), 'xticklabels')
    P.yticks ([])
    P.title ('[kappa=%4.2f, mu=%3.1f, N=%d]' % (kappa,180.*mu/np.pi,n))
    print 'kappa = %4.2f, mu = %1.2f, N = %d, p = %10.10f /100' % (
        kappa,mu,n,100.*p)
    return (kappa,mu,p,n)

def plot_psth(events, duration = None, binwidth = None, nbins = None, meanrate = True, smooth = True, **args):
    """
    Plots a post-stimulus time histogram (PSTH).

    Calls psth. See spikelib module for details.

    Optional arguments:  'binwidth' (s), 'nbins', 'meanrate', 'smooth'.
    Defaults: binwidth = 1ms, nbins = None, meanrate = True, smooth = True.
    """

    # compute histogram
    values, edges = spk.psth(events, duration, binwidth, nbins, meanrate, smooth)
    edges = edges[:-1]
    # update user-defined plot arguments
    plotargs = dict(figure=P.gcf())
    plotargs.update(args)

    if meanrate:
        ylabel ='spike rate (Hz)'
    else:
        ylabel ='spike count'

    # plot it!
    P.plot(edges,values,**plotargs)
    P.xlabel('time (%s)'%events.i.units)
    P.ylabel(ylabel)

    return values, edges

def plot_sta(spks, stim, movie, tslices = 10, dim = None, roi = None, spkrate = True, **plotargs):
    """
    Plot spike-triggered average derived spatiotemporal receptive field (STRF).

    Calls sta. See spikelib module for details.

    Optional arguments: tslices (scalar number of frames, or list of specific time-slices),
    roi (bounding box region-of-interest), spkrate (normalize sta in units of spikes/s).
    """
    # compute sta
    sta =  spk.sta(spks, stim, movie, tslices = tslices, dim = dim, roi = roi, spkrate = spkrate)
    
    # default plot arguments
    fig = dict(figure=P.gcf())
    plotargs.update(fig)

    # get current axis
    ax = plotargs['figure'].gca()

    # plot kernel
    if np.isscalar(tslices):
        tslices = np.multiply(range(tslices),stim.i.sweeptimeMsec)
    else:
        tslices = np.multiply(tslices,1000/stim.i.refreshrate)    
    nbins = len(tslices)
    if ax.numCols != nbins: P.clf()
    (vmin, vmax) = (sta.min(), sta.max())
    P.ioff()
    for i, j in enumerate(tslices):
        P.subplot(ax.numRows, nbins, i+1)#, sharex=ax, sharey=ax)# link axes
        if nbins == 1: P.imshow(sta, vmin=vmin, vmax=vmax, **plotargs)
        else: P.imshow(sta[nbins-i-1], vmin=vmin, vmax=vmax, **plotargs)
        P.xticks([]); P.yticks([])
        P.xlabel('t=%.0f ms'%j)
    P.show()
    return sta

def plot_saccade_scatter(sac):
    dx = sac.xf - sac.xi
    dy = sac.yf - sac.yi
    radii = sac.amplitude
    thetas = np.arctan2(dy, dx)
    ax = P.gca(polar=True)
    ax.scatter(thetas,radii,alpha=.5,s=sac.vpeak/10.0)

def plot_saccade_hist(sac,bins=36):
    dx = sac.xf - sac.xi
    dy = sac.yf - sac.yi
    radii = sac.amplitude
    thetas = np.arctan2(dy, dx)
    theta_bins = np.linspace(-np.pi,np.pi,bins+1)
    theta_hist = np.histogram(thetas,bins=theta_bins)
    #XXX: I tried using multiplot here, but failed to set polar=True and gave
    # up. -pi
    P.subplot(221,polar=True)
    P.jet()
    P.title("individual saccades")
    # XXX: there's an interesting artifact if we use sac.amplitude as the
    # sac_length in the scatter plot below
    sac_length = np.sqrt((dx**2+dy**2))
    P.scatter(thetas,sac_length,alpha=.5,c=sac.amplitude,s=sac.vpeak/10.0)
    P.colorbar()
    P.subplot(222,polar=True)
    P.title("Directional histogram")
    bars = P.bar(theta_bins[:-1],theta_hist[0],width=(2*np.pi)/bins, alpha=.5)
    P.subplot(223)
    P.scatter(dx,dy,alpha=.5,c=sac.amplitude,s=sac.vpeak/10.0)

    #P.hist(sac_length,bins=100)
    #P.xlabel("saccade lengths (pixels)")
    P.subplot(224)
    P.hist(sac.vpeak,bins=100)
    P.xlabel("saccade peak velocities")

def rasters(neurons,stims,varargins):
    """
    raster(neurons,stims,varargins)

    Given a cell array of neuron structures, and a stim structure or a
    cell array of stim structures, generates a multi-unit raster plot.

    Use arrow keys to pan and zoom. Use < and > keys to page left and
    right. If an 'axishandle' is specified, doesn't create a new figure
    window.

    Optional arguments: 'range', 'spacing', 'spikeheight', 'duration',
                        'showlines', 'units'.
    """
    pass


def rasterspophist(neurons,stims,varargins):
    """
    rasterspophist(neurons,stims,varargins)

    Given a cell array of neuron structures, and a stim structure or a
    cell array of stim structures, generates a multi-unit raster plot
    and a population histogram below it.

    Use arrow keys to pan and zoom. Use < and > keys to page left and right
    If an 'axishandle' is specified, doesn't create a new figure window

    Optional arguments: 'range', 'spacing', 'spikeheight', 'duration',
                        'showlines', 'units'.
    """
    pass


def rastersweeps(neuron,stim,varargin):
    """
    rastersweeps(neuron,stim,varargin)

    Generates and labels multiple raster plots, one per stimulus sweep
    Expects a neuron structure (see loadneuron) and a stim structure (see
    loadstim). Plots all sweeps specified in stim.sweeptable, unless
    optional argument 'sweeps' is passed, which specifies in a 1d array
    the sweeps you want to create raster plots for.
    """
    pass


def rawplot(fname):
    """
    rawplot(fname)

    Plot time-series data from an UInt16 binary file. Arrow keys zoom and pan.
    """
    pass

def plot_gaze_pos_hist():
    pass

def plot_fixation_vs_link_vs_photodiode(bdf,eye):
    ax1 = P.subplot(211)
    nf = bdf.frame.v < 0x7FFF
    bdft0 = bdf.frame[nf][0] +1.0 # minus extra second, because there was a 1 second prestimulus
    P.plot(bdf.frame[nf]-bdft0,bdf.frame.v[nf]+.08,'bo',label='frames')
    P.plot(bdf.Ana1.time()-bdft0,bdf.Ana1,label="Photodiode")
    P.legend()
    P.subplot(212,sharex=ax1)
    eyet0 = eye.frames[0]
    P.plot(eye.frames-eyet0, eye.frames.v,'bo',label='frames')
    P.plot(eye.l.fixations.start-eyet0,np.ones_like(eye.l.fixations.start),'go',label="Fixations")
    P.plot(eye.link_fixations-eyet0,np.ones_like(eye.link_fixations)+1,'ro',label="Link Fixs")
    P.legend()

def normed_hist(vals, bins=10, cumulative=False, orientation='', ax=None, **kwargs):
    if ax == None:
        ax = P.gca()
    H,edges = np.histogram(vals, bins=bins )
    h = H.astype(float)  / H.sum().astype(float)
    if 'color' not in kwargs:
        kwargs['color'] = ax._get_lines._get_next_cycle_color()
    if cumulative:
        h = h.cumsum()

    if orientation.lower() in ['h', 'hor', 'horizontal']:
        bar = ax.barh
        kwargs.update({"height":edges[1]-edges[0]})
    else:
        bar = ax.bar
        kwargs.update({"width":edges[1]-edges[0]})

    bar(edges[:-1],h,**kwargs)
