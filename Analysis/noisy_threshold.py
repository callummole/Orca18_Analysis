import numpy as np
from scipy.interpolate import interp1d
from scipy import stats
import scipy.optimize
import itertools

def rt_likelihood(ts, vals, threshold, noise, lag=0.4):
    # TODO: Lag as a distribution!
    # TODO: Fumbling
    # TODO: Compute as a continuous process (for constant or linear segments)
    dt = np.median(np.diff(ts)) # Hack
    # TODO: Check if the dt belongs within the sqrt and whether it should
    # be multiplied or divided!!?!
    std = noise*np.sqrt(dt)
    alive = np.exp(np.cumsum(stats.norm.logcdf(threshold, loc=vals, scale=std)))
    #assert np.all(np.isfinite(alive))
    prob = np.ediff1d(1 - alive, to_begin=0)
    interp = interp1d(ts + lag, prob, fill_value=(np.nan, alive[-1]), bounds_error=False)
    interp.cumulative = 1 - alive
    return interp

from reconstruct_trajectory import get_untouched_trajectories
import matplotlib.pyplot as plt
def plot_trajs():
    data = get_untouched_trajectories()
    for g, td in data.items():
        traj = td['traj']
        otraj = td['otraj']

        plt.plot(traj['ts'], traj['yawrate_seconds'])
        plt.plot(otraj['ts'], otraj['yawrate_seconds'])
        plt.show()

# Derived using quadratic approximations of circular
# trajectories and path edges. Trusting numpy to turn
# imaginary solutions to nans and thus get rid of them.
# Could be a lot nicer.
def lane_ttlc(a0, b, v, d0, w):
    sol1 = -a0/b - np.sqrt(a0**2*v - 2*b*d0 - 2*b*w)/(b*np.sqrt(v))
    sol2 = -a0/b + np.sqrt(a0**2*v - 2*b*d0 - 2*b*w)/(b*np.sqrt(v))
    sol1 = np.sqrt(sol1)**2
    sol2 = np.sqrt(sol2)**2
    sol = np.fmin(sol1, sol2)
    sol[np.abs(d0) > np.abs(w)] = 0.0
    return sol

# Hacks the positive real solution if one exists
def path_ttlc(a0, b, v, d0, w):
    return np.fmin(
        lane_ttlc(a0, b, v, d0, w),
        lane_ttlc(a0, b, v, d0, -w)
        )

def fit_trials(trials):
    speed = 8
    w = (3 - 1.5)/2

    # TODO: This really shouldn't be done here
    for trial in trials:
        radius = trial['otraj'].radius.iloc[0]
        direction = trial['otraj'].bend.iloc[0]
        traj = trial['traj']
        a0 = np.radians(traj['world_yaw'].values - traj['tangent_yaw'].values)
        a0 = (a0 + np.pi) % (2*np.pi) - np.pi
        b = (np.radians(traj['yawrate_seconds'].values) - direction*speed/radius)
        
        d0 = -traj['lane_bias'].values
        
        b *= direction
        a0 *= direction
        traj['ttlc'] = path_ttlc(a0, b, speed, d0, w)
        traj['thresholdee'] =(1/np.maximum(1e-5, traj['ttlc'].values))**(1/3)
        #traj['thresholdee'] = np.abs(traj['lane_bias'].values)
        
        continue
        ax = plt.subplot(2,1,1)
        plt.plot(traj['ts'], np.degrees(d0), alpha=0.5)
        #plt.plot(traj['ts'], traj['lane_bias'])
        #plt.axhline(-w, color='black')
        #plt.axhline(w, color='black')
        #plt.ylim(-w*2, 2*w)

        #plt.plot(traj['ts'], traj['ttc'])
        plt.subplot(2,1,2, sharex=ax)
        plt.plot(traj['ts'], 1/traj['ttlc'], alpha=0.5)
        #plt.plot(traj['ts'], 1.0/traj['ttc'])
        #if direction == 1:
        plt.show()

    trajs, tots = list(zip(*((v['traj'], v['takeover_time']) for v in trials)))
    liks = [
            lambda *args, traj=traj, tot=tot: rt_likelihood(traj['ts'].values, traj['thresholdee'].values, *args)(tot)
            for traj, tot in zip(trajs, tots)
            ]

    def loss(args):
        ls = np.array([l(*np.exp(args)) for l in liks])
        return -np.sum(np.log(ls + 1e-9))
    
    init = np.array([1.0, 1.0])
    res = scipy.optimize.minimize(loss, np.log(init))
    res.x = np.exp(res.x)
    print(res)
    print((2*(len(init) + res.fun))/len(liks))
    return res.x

# Transforms between TTLC and yaw rate offset on circular
# trajectories, assuming centered initial position. Based on
# small angle approximation of cosine, but should be in practice
# very accurate on reasonable yaw rates.
#
# b is the yaw rate offset in radians (negative is understeering)
# t is the TTLC (negative is understeering)
# w is half the road width
# r is the path radius
# v is the velocity (m/s)
def ttlc_from_offset(b, w, r, v):
    return np.sign(b)*np.sqrt(w*(2*r + np.sign(b)*w)/(np.abs(b)*r*v))

def offset_from_ttlc(t, w, r, v):
    return np.sign(t)*w*(2*r + np.sign(t)*w)/(r*t**2*v)

def fit_per_participant():
    data = get_untouched_trajectories()
    for s, sds in itertools.groupby(data.items(), lambda kv: kv[0][0]):
        if s in [4, 19]: continue # HACK!!
        gs, sds = list(zip(*sds))
        sds = [t for t in sds if (
            t['traj']['ts'].values[0] <= t['takeover_time'] # Remove very early responses. TODO: Should be done before!
            and t['cogload'] == "None"
        )]

        params = fit_trials(sds)
        print(params)

        sab_preds = {}
        fits = []
        tots = []
        sabs = []
        for sd in sds:
            if sd['otraj'].design.iloc[0] != 'balanced': continue
            traj = sd['traj']
            lik = rt_likelihood(traj.ts.values, traj['thresholdee'], *params)
            sab = sd['otraj']['sab'].iloc[0]
            if sd['otraj'].design.iloc[0] == 'balanced' and sab not in sab_preds:
                sab_preds[sab] = lik
            fits.append(lik.x[np.argmax(lik.y)] - sd['onset_time'])
            tots.append(sd['takeover_time'] - sd['onset_time'])
            sabs.append(sab)

            plt.plot(lik.x, lik.y)
            plt.axvline(sd['takeover_time'])
            plt.twinx()
            plt.plot(traj.ts.values, traj['thresholdee'], color='red', alpha=0.5)
            plt.ylim(0, 2)
            plt.show()
        radius = sds[0]['otraj']['radius'].iloc[0]
        speed = 8
        margin = 3/2
        onset = 6
        
        sab_modes = []
        for sab, lik in sab_preds.items():
            sab_modes.append((
                    lik.x[lik.cumulative.searchsorted(0.05)],
                    lik.x[lik.cumulative.searchsorted(0.25)],
                    #lik.x[np.argmax(lik.y)],
                    lik.x[lik.cumulative.searchsorted(0.5)],
                    lik.x[lik.cumulative.searchsorted(0.75)],
                    lik.x[lik.cumulative.searchsorted(0.95)],
                    ))
        sab_modes = np.array(sab_modes) - onset
        

        ttlcs = ttlc_from_offset(np.radians(list(sab_preds.keys())), margin, radius, speed)
        color = plt.plot(np.abs(ttlcs), sab_modes[:,(len(sab_modes))//2])[0].get_color()
        plt.fill_between(np.abs(ttlcs), sab_modes[:,0], sab_modes[:,-1], color=color, alpha=0.3)
        plt.fill_between(np.abs(ttlcs), sab_modes[:,1], sab_modes[:,-2], color=color, alpha=0.3)
        
        ttlcs = ttlc_from_offset(np.radians(sabs), margin, radius, speed)
        plt.plot(np.abs(ttlcs), tots, '.')

        
        plt.plot([0, 10], [0, 10], color='black')
        #plt.xlabel("Real takeover time (seconds since onset)")
        #plt.ylabel("Predicted takeover time (seconds since onset)")
        plt.show()


def demo():
    import matplotlib.pyplot as plt
    dt = 1/60.0
    noise = 1.0
    threshold = 1.0
    dur = 10

    val = lambda t: 0.1*t

    ts = np.arange(0, dur, dt)
    #plt.plot(ts, val(ts))
    for i, noise in enumerate([0.5, 1.0, 2.0]):
        dxs = np.linspace(0.1, 10, 100)
        lows = []
        highs = []
        modes = []
        for dx in dxs:
            val = dx*ts
            likelihood = rt_likelihood(ts, val, threshold, noise)
            mode = likelihood.x[np.argmax(likelihood.y)]
            low = likelihood.x[likelihood.cumulative.searchsorted(0.25)]
            high = likelihood.x[likelihood.cumulative.searchsorted(0.75)]
            #plt.plot(ts, likelihood(ts))
            modes.append(mode)
            lows.append(low)
            highs.append(high)
        
        plt.fill_between(1/dxs, lows, highs, color=f"C{i}", alpha=0.2)
        plt.plot(1/dxs, modes, color=f"C{i}", label=f"Noise {noise}")

    plt.legend()
    plt.show()

if __name__ == '__main__':
    fit_per_participant()
    #plot_trajs()
    #demo()
