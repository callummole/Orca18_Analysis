import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.interpolate
import scipy.optimize
import scipy.stats
import scipy.special

from reconstruct_trajectory import get_automation_trajectory
expdata = pd.read_parquet(Path(__file__).parent/'../Data/orca_raw_longformat.parq')

trialcols = "ppid radius yawrate_offset cogload block count".split()

def get_track(radius, is_mirrored):
    myStraight  = simTrackMaker.lineStraight(startpos = [0,0], length= 16)
    x_dir = -1 if is_mirrored else 1
    myBend = simTrackMaker.lineBend(startpos = myStraight.RoadEnd, rads=radius, x_dir=x_dir, road_width=3.0) 
    Course_midline = np.vstack((myStraight.midline, myBend.midline))
    return Course_midline, myBend.CurveOrigin

def get_untouched_trajectories():
    cpath = 'untouchedhackcache.pickle'
    try:
        return pickle.load(open(cpath, 'rb'))
    except FileNotFoundError:
        pass
    
    data = {}
    for tg, td in expdata.groupby(trialcols):
        traj = get_automation_trajectory(td)
        #plt.plot(traj.World_x, traj.World_z)
        #plt.plot(td.World_x, td.World_z)
        td['ts'] = np.cumsum(td.dt)
        try:
            onset_t = td.ts.values[td.timestamp_trial.values.searchsorted(td.OnsetTime.iloc[0])]
        except IndexError:
            onset_t = np.inf

        try:
            tot = td.ts[td.AutoFlag != 1].iloc[0]
        except IndexError:
            tot = np.inf
        data[tg] = {
                'traj': traj,
                'takeover_time': tot,
                'onset_time': onset_t,
                **dict(zip(trialcols, tg))
                }
        #auto = traj.ts.values < tot
        #plt.plot(traj.ts[auto], lane_bias[auto])
    pickle.dump(data, open(cpath, 'wb'))
    return data
    #plt.plot(midline[:,0], midline[:,1])

bias_at_takeover = []

untouched = get_untouched_trajectories()
for g, td in untouched.items():
    traj = td['traj']
    
    traj['lane_bias_change'] = np.gradient(traj.lane_bias.values)/traj.dt.values
    closing_edge = np.sign(traj.lane_bias_change.values)*1.5
    margin = closing_edge + traj.lane_bias
    traj['ttc'] = margin/traj.lane_bias_change


times, trials, onsets, offsets = zip(*[
    (v['takeover_time'], v['traj'], v['onset_time'], v['yawrate_offset'])
    for g, v in untouched.items()
    if g[3] != 'None'
    ])
times = np.array(times)

hackinf = 5.0
def thtimer(ts, val):
    valid = np.isfinite(val) # HACK!
    ts = ts[valid]
    val = val[valid]
    max_seen = np.maximum.accumulate(val)
    changes = np.flatnonzero(np.ediff1d(max_seen) > 0) + 1
    changes = [0] + list(changes)
    tint = scipy.interpolate.interp1d(max_seen[changes], ts[changes],
        fill_value=(ts[0], np.inf), bounds_error=False)
    return tint

timers = [thtimer(t.ts.values, 1.0/t.ttc.values)
    for t in trials]

densitiers = []
for t in trials:
    ts = t.ts.values
    val = 1.0/t.ttc.values
    valid = np.isfinite(val) # HACK!
    ts = ts[valid]
    val = val[valid]
    max_seen = np.maximum.accumulate(val)
    changes = np.flatnonzero(np.ediff1d(max_seen) > 0) + 1
    changes = [0] + list(changes)
    def dens(thm, ths, rtm, rts, fudge=0.001, ts=ts, val=val, changes=changes):
        thd = scipy.stats.lognorm(s=ths, scale=thm)
        rtd = scipy.stats.lognorm(s=rts, scale=rtm)
        fudged = scipy.stats.expon(scale=1.0/fudge)
        decs = thd.cdf(val[changes])
        decs = np.ediff1d(decs, to_begin=decs[0])
        #decs *= fudged.sf(ts[changes])
        noresp = 1 - np.sum(decs)
        #noresp -= fudged.cdf(ts[-1])
        noresp += np.sum(rtd.sf(ts[-1] - ts[changes])*decs)
        def d(t):
            sts = t - ts[changes].reshape(-1, 1)
            ps = (rtd.cdf(sts) - rtd.cdf(sts - 1/60))*decs.reshape(-1, 1)
            #ps += fudged.pdf(t) # TODO: Not right at all
            finites = np.isfinite(t)
            ps = np.sum(ps, axis=0)*finites + noresp*(~finites)
            res = np.reshape(ps, np.shape(t))
            return res
        return d
    densitiers.append(dens)
    #d = dens(0.5, 1.0, 0.5, 0.5)(t.ts.values)
    #plt.plot(t.ts.values, d)
    #plt.twinx()
    #plt.plot(t.ts.values, t.lane_bias.values, color='red')
    #plt.show()

#times[~np.isfinite(times)] = hackinf # Hack!
def predict_takeovers(th, rt):
    preds = []
    for timer in timers:
        preds.append(timer(th) + rt)
    return np.array(preds)

def param_likelihood(x):
    print(x)
    return -np.sum(np.log([
        densitier(*x)(t) + 1e-3 for densitier, t in zip(densitiers, times)
        ]))


#loss = lambda x: times - predict_takeovers(*np.exp(x))
#wtf = scipy.optimize.least_squares(loss, np.log((0.4, 0.5)), loss='soft_l1')

#params = [0.04204588, 0.04009958, 5.32956664, 0.54384475]
#params = [0.29232104, 1.0171227 , 0.36539256, 1.8746681]
params = [0.1417864243934838, 0.30644850892735503, 0.7559849093141972, 0.6570338953631301]
params = [0.19827915328181286, 0.13989151611317846, 0.46692017625776416, 0.3286570149898443]
mangle = lambda x: list(np.log(x)) #+ [scipy.special.logit(x[-1])]
demangle = lambda x: list(np.exp(x)) #+ [scipy.special.expit(x[-1])]

"""
wtf = scipy.optimize.minimize(
        lambda x: param_likelihood(demangle(x)),
        mangle((0.29232104, 1.0171227 , 0.36539256, 1.8746681)),
        method='powell')
wtf.x = demangle(wtf.x)
print(wtf)
params = wtf.x
"""

thm, ths, rtm, rts, *_ = params
th_mode = np.exp(np.log(thm) - ths**2)
rt_mode = np.exp(np.log(rtm) - rts**2)
thd = scipy.stats.lognorm(s=ths, scale=thm)
rtd = scipy.stats.lognorm(s=rts, scale=rtm)

"""
rng = np.linspace(0.0, thd.ppf(0.9), 1000)
plt.plot(rng, thd.pdf(rng))
plt.axvline(th_mode, color='red')
plt.figure()
rng = np.linspace(0.0, rtd.ppf(0.9), 1000)
plt.plot(rng, rtd.pdf(rng))
plt.axvline(rt_mode, color='red')
plt.show()
"""

for dens, traj, tot, ot in zip(densitiers, trials, times, onsets):
    d = dens(*params)(traj.ts.values)
    plt.plot(traj.ts.values, dens(*params)(traj.ts.values), color='blue', label="Crossing probability")
    plt.axvline(tot, color='black', label='Participant takeover')
    plt.axvline(ot, color='black', linestyle="dashed", label='Offset onset')
    plt.ylabel("Crossing probability")
    plt.title(f"Non-cross probability {dens(*params)(np.inf):.2f}")
    plt.legend()
    plt.twinx()
    plt.plot(traj.ts, 1.0/traj.ttc, color='red')
    plt.xlabel("Time (seconds)")
    plt.ylabel("1/TTC (1/seconds)")
    plt.show()

times -= onsets

times[~np.isfinite(times)] = hackinf
def loss(x):
    preds = predict_takeovers(*np.exp(x)) - onsets
    preds[~np.isfinite(preds)] = hackinf
    return preds - times
    
wtf = scipy.optimize.least_squares(loss, np.log((0.4, 0.5)), loss='soft_l1')
wtf.x = np.exp(wtf.x)
preds = predict_takeovers(th_mode, rt_mode) - onsets

preds[~np.isfinite(preds)] = hackinf
plt.scatter(times, preds, c=offsets, marker='.', alpha=0.6)
#plt.plot([6, 14], [6, 14])
plt.ylabel("Predicted time since offset (seconds)")
plt.xlabel("Measured time since offset (seconds)")
plt.show()

"""
for g, td in get_untouched_trajectories().items():
    if g[3] != 'None':
        continue
    traj = td['traj']
    dt = 1/60.0
    lane_bias = scipy.interpolate.interp1d(traj.ts.values, traj.lane_bias.values)(nt)

    traj['lane_bias_change'] = np.gradient(traj.lane_bias)/dt
    closing_edge = np.sign(lane_bias_change)*1.5
    margin = closing_edge + lane_bias
    ttc = margin/lane_bias_change

    plt.plot(nts, 1.0/ttc)
    plt.twinx()
    plt.plot(nts, traj.lane_bias, color='red')
    plt.axvline(td['takeover_time'], color='black')
    plt.show()
    tot = td['takeover_time']
    if not np.isfinite(tot):
        continue
    bias_at_takeover.append(
            traj.lane_bias_change.values[traj.ts.values.searchsorted(tot)]
            )

plt.hist(np.abs(bias_at_takeover))
"""
plt.show()

