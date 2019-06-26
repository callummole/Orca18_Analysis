import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path

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

def get_automation_trajectory(td):
    traj = trajs[td.AutoFile.iloc[0]]
    
    n = min(len(td), len(traj))
    traj = traj.iloc[:n]
    td = td.iloc[:n]
    
    is_mirrored = np.sum(np.abs(td.World_x - td.world_x_mirrored)) > 1e-6
    
    orig_traj = traj
    traj = traj.copy()
    
    onset_i = td.timestamp_trial.values.searchsorted(td.OnsetTime.iloc[0])
    yaw_rate = traj.YawRate_seconds.values
    yaw_rate[onset_i:] += td.yawrate_offset.iloc[0]
    yaw_rate[:] = np.clip(yaw_rate, -MAX_YAW_RATE, MAX_YAW_RATE)
    if is_mirrored:
        yaw_rate *= -1

    xs = []
    zs = []
    bs = []
    
    b = td.WorldYaw.iloc[0]
    x = td.World_x.iloc[0]
    z = td.World_z.iloc[0]
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
        WorldYaw=bs,
        World_x=xs,
        World_z=zs,
        YawRate_seconds=yaw_rate,
        is_mirrored=is_mirrored
        ))
    
    midline, origin = get_track(td.radius.iloc[0], traj.is_mirrored.iloc[0])
    midshape = LineString(midline)
    
    pos = traj[['World_x', 'World_z']].values
    dists = [midshape.project(Point(*p)) for p in pos]
    cps = np.array([midshape.interpolate(d).xy for d in dists]).reshape(-1, 2)
    errors = pos - cps
    dirs = np.sign(np.einsum('ij,ij->i', cps - origin, errors))
    lane_bias = np.sign(dirs)*np.linalg.norm(errors, axis=1)
    traj['lane_bias'] = lane_bias

    return traj
