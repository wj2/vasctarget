"""
This is a module for the simulation of spikt trains.
"""

import numpy as np
import dsplib as dsp
import spikelib as spk
import plotlib as plt
import model
import utils
import data

from matplotlib.pyplot import *
reload(utils)
reload(data)
reload(model)


def simulate_poisson(rate,repeats=None,samplingrate=None):
    """Simulate Poisson point process"""
    return simulate_gamma(rate,k=1,repeats=repeats,samplingrate=samplingrate)

def simulate_gamma(rate,k,repeats=None,samplingrate=None):
    """Simulate Gamma point process"""
    if samplingrate is None:
        samplingrate = rate.i.samplingrate
    oversampling = int(samplingrate/rate.i.samplingrate)
    assert float(samplingrate)/rate.i.samplingrate == oversampling
    timings = []

    if rate.has_trials():
        assert repeats is None
        for trial in rate:
            tmg = simulate_gamma(trial,k,repeats,samplingrate)
            timings.append(tmg)
    else:
        nt = oversampling*len(rate)
        density = np.repeat(rate/samplingrate,k*oversampling)
        time = np.arange(k*nt,dtype=float)/(k*samplingrate)
        kstart = np.random.randint(0,k,repeats or 1)
        for tr in range(repeats or 1):
            tmg = time[density > np.random.rand(k*nt)][kstart[tr]::k]
            if repeats is None:
                timings.extend(tmg.tolist())
            else:
                timings.append(tmg.tolist())
    return data.Events(timings)


def simulate_phase(nt,samplingrate,f,sf=0):
    time = np.arange(nt,dtype=float)/(samplingrate)
    if sf != 0:
        noise = data.TimeSeries(np.random.randn(nt),samplingrate)
        phase = data.TimeSeries(np.angle(dsp.analytic_signal(f,sf,noise)),samplingrate)
    else:
        phase = data.TimeSeries(utils.smod(2*np.pi*f*time),samplingrate)
    return phase

def simulate_modulated_gamma(meanrate,k,f,kappa,mu=0,sf=0,repeats=None,samplingrate=None,random_phase=True):
    """Simulate modulated Gamma point process"""
    if samplingrate is None:
        samplingrate = meanrate.i.samplingrate
    oversampling = int(samplingrate/meanrate.i.samplingrate)
    assert float(samplingrate)/meanrate.i.samplingrate == oversampling
    timings = []
    phases = []
    phi = 0

    nt = oversampling*len(meanrate)
    phase = simulate_phase(nt,samplingrate,f,sf)

    for trial in range(repeats or 1):
        if not np.mod(trial-9,10):
            print trial+1, ' ',
            utils.flush() # flush stdout
        if kappa != None:
            if random_phase:
                phi0=2*np.pi*np.random.rand()
            else:
                phi0=0
        tphase = data.TimeSeries(utils.smod(phase+phi0),phase.i.samplingrate)
        rate =  data.TimeSeries(np.repeat(meanrate,oversampling)*2*np.pi*model.mises(tphase,kappa,mu),samplingrate)
        tmg = simulate_gamma(rate,k=2,repeats=None,samplingrate=samplingrate)
        if repeats is not None:
            timings.append(tmg)
            phases.append(tphase[tmg])
        else:
            timings.extend(tmg)
            phases.extend(tphase[tmg])
    return data.Events(timings,phase=phases)




if __name__ == '__main__':
    # simple tests
#     rate = data.TimeSeries(10*np.random.rand(10**3),samplingrate=10**3)
#     sim_spikes = simulate_gamma(rate,k=2,repeats=10**3,samplingrate=10**3)
#     sim_rate = sim_spikes.asrate(10**3,nbins=len(rate)).mean(0)
#     np.testing.assert_almost_equal(rate,sim_rate,decimal=-1)


#     rate = data.TimeSeries(100*np.ones(10**3),samplingrate=10**3)
#     sim_spikes = simulate_gamma(rate,k=3,repeats=10**3,samplingrate=10**4)
#     isi = np.concatenate([np.diff(tmg) for tmg in sim_spikes.iter_dim()])
#     k,xmean = model.fit_gamma(isi)
#     np.testing.assert_almost_equal(k,3,decimal=1)

    kappa0 = .3
    f0 = 50
    rate = data.TimeSeries(40*np.random.rand(10**5),samplingrate=10**3)
    sim_spikes = simulate_modulated_gamma(rate,k=2,kappa=kappa0,f=f0,sf=2,samplingrate=10**3)

    close('all')
    plt.plot_phasedist(sim_spikes.phase)
    show()

    kappa,mu,p,n = model.phasedist(sim_spikes.phase)
    # sim_rate = sim_spikes.asrate(10**3,nbins=len(rate)).mean(0)
    np.testing.assert_almost_equal(kappa,kappa0,decimal=0)
    f, Os = spk.oscore(sim_spikes,40,80,10**3)
    print 'frequency: ',np.mean(f),'   oscillation score: ',np.mean(Os)
    # np.testing.assert_almost_equal(mean(f),f0,decimal=-1)
    # np.testing.assert_almost_equal(mean(Os),4*kappa0**2,decimal=-1)

    sf=2
    sdt = 1./(2.*np.pi*sf)
    simrate = sim_spikes.asrate(1000,width=.001)
    slowrate = dsp.smoothg(simrate,std=sdt,normed='energy',ns=5)

    asig = dsp.analytic_signal(50,2,simrate)

    nonzero = slowrate>0
    coherence = abs(asig)**2
    coherence[nonzero] /= slowrate[nonzero]**2

