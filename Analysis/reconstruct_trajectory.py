import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path
import pickle

SPEED = 8.0
DT = 1/60.0
MAX_YAW_RATE = 35.0

trajs = Path(__file__).parent.glob('../Data/Midline*.csv')
trajs = {t.name: pd.read_csv(t) for t in trajs}

import simTrackMaker
from shapely.geometry import LineString, Point

def get_track(radius, is_mirrored):
    myStraight  = simTrackMaker.lineStraight(startpos = [0,0], length= 16)
    x_dir = -1 if is_mirrored else 1
    myBend = simTrackMaker.lineBend(startpos = myStraight.RoadEnd, rads=radius, x_dir=x_dir, road_width=3.0) 
    Course_midline = np.vstack((myStraight.midline, myBend.midline))
    return Course_midline, myBend.CurveOrigin

f = lambda t: np.exp(-1/t)*(t > 0)
_smooth_step = lambda t: f(t)/(f(t) + f(1 - t))
eps = np.finfo(float).eps
def smooth_step(t):
    out = _smooth_step(t)
    out[t <= 0+eps] = 0.0
    out[t >= 1-eps] = 1.0
    return out

transition_duration = .5

trialcols = "ppid radius sab cogload trialn".split()

def get_automation_trajectory(td):
    traj = trajs[td.autofile.iloc[0]]
    
    n = min(len(td), len(traj))
    traj = traj.iloc[:n]
    td = td.iloc[:n]
    
    #is_mirrored = np.sum(np.abs(td.world_x - td.world_x_mirrored)) > 1e-6
    direction = td.bend.iloc[0]
    is_mirrored = direction < 0
    orig_traj = traj
    traj = traj.copy()
    
    onset_i = td.timestamp_trial.values.searchsorted(td.onsettime.iloc[0])
    yaw_rate = traj.YawRate_seconds.values
    yaw_rate += smooth_step((td.timestamp_trial.values - td.onsettime.iloc[0])/transition_duration)*td.sab.iloc[0]
    yaw_rate[:] = np.clip(yaw_rate, -MAX_YAW_RATE, MAX_YAW_RATE)
    if is_mirrored:
        yaw_rate *= -1

    xs = []
    zs = []
    bs = []
    
    b = td.world_yaw.iloc[0]
    x = td.world_x.iloc[0]
    z = td.world_z.iloc[0]
    dt = td.dt.values

    for i in range(len(yaw_rate)):
        if i < len(yaw_rate) - 1 and i > 0:
            b += yaw_rate[i]*dt[i - 1]
        xs.append(x)
        zs.append(z)
        bs.append(b)
        if i == len(yaw_rate) - 1: continue
        x += np.sin(np.radians(b))*SPEED*dt[i+1]
        z += np.cos(np.radians(b))*SPEED*dt[i+1]

    traj = pd.DataFrame(dict(
        ts=np.cumsum(td.dt.values),
        dt=td.dt,
        world_yaw=bs,
        world_x=xs,
        world_z=zs,
        yawrate_seconds=yaw_rate,
        is_mirrored=is_mirrored
        ))
    
    midline, origin = get_track(td.radius.iloc[0], traj.is_mirrored.iloc[0])
    midshape = LineString(midline)
    
    pos = traj[['world_x', 'world_z']].values
    dists = [midshape.project(Point(*p)) for p in pos]
    cps = np.array([midshape.interpolate(d).xy for d in dists]).reshape(-1, 2)
    
    off_origin = cps - origin
    norma = lambda a: (a + np.pi) % (2*np.pi) - np.pi
    # TODO: Could compute also the linear part correctly!
    traj['tangent_yaw'] = np.degrees(
            (np.arctan2(off_origin[:,0], off_origin[:,1]) + np.pi*(float(is_mirrored) + 1/2))
            )
    
    """
    import matplotlib.pyplot as plt
    plt.plot(traj['ts'], np.degrees(norma(np.radians(traj['world_yaw'] - traj['tangent_yaw']))), color='green')
    plt.show()
    """

    # This is just bizarre!
    errors = pos - cps
    dirs = np.sign(np.einsum('ij,ij->i', off_origin, errors))
    lane_bias = np.sign(dirs)*np.linalg.norm(errors, axis=1)
    traj['lane_bias'] = lane_bias

    return traj

def _get_untouched_trajectories():
    cpath = 'untouchedhackcache.pickle'
    try:
        return pickle.load(open(cpath, 'rb'))
    except FileNotFoundError:
        pass
    expdata = pd.read_parquet(Path(__file__).parent/'../Data/orca_rerun.parq')
    data = {}
    for tg, td in expdata.groupby(trialcols):
        traj = get_automation_trajectory(td)
        #plt.plot(traj.World_x, traj.World_z)
        #plt.plot(td.World_x, td.World_z)
        td['ts'] = np.cumsum(td.dt)
        try:
            onset_t = td.ts.values[td.timestamp_trial.values.searchsorted(td.onsettime.iloc[0])]
        except IndexError:
            onset_t = np.inf

        try:
            tot = td.ts[td.autoflag != 1].iloc[0]
        except IndexError:
            tot = np.inf
        data[tg] = {
                'traj': traj,
                'otraj': td.copy(),
                'takeover_time': tot,
                'onset_time': onset_t,
                **dict(zip(trialcols, tg))
                }
        #auto = traj.ts.values < tot
        #plt.plot(traj.ts[auto], lane_bias[auto])
    pickle.dump(data, open(cpath, 'wb'))
    return data

def get_untouched_trajectories():
    untouched = _get_untouched_trajectories()
    for g, td in untouched.items():
        traj = td['traj']
        traj['lane_bias_change'] = np.gradient(traj.lane_bias.values)/traj.dt.values
        closing_edge = np.sign(traj.lane_bias_change.values)*1.5
        margin = closing_edge - traj.lane_bias.values
        ttc = margin/(traj.lane_bias_change.values)
        ttc[ttc < 1e-5] = 1e-5 # Hack!
        traj['ttc'] = ttc
        # Hack: There's something funny going on in the beginning
        # causing peaks in the ttc
        traj = traj.iloc[traj.ts.values.searchsorted(2.5):]
        td['traj'] = traj
    return untouched
