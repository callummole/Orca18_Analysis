import numpy as np
import pandas as pd 
import TrackSimulation as tsim
				

#####recreates the experiment, with an off-tangent failure########


track = pd.read_csv("..//Data//track_with_edges_orca_80.csv")

def remove_straight(cols):

	road = np.array([track[cols[0]].values, track[cols[1]].values]).T
	road[:,1] = road[:,1] - 16 #remove straight.
	mask = road[:,1]>= 0
	road = road[mask,:]	
	return road

midline =  remove_straight(['midlinex','midlinez'])


bend_yr = np.rad2deg(8.0 / 80.0)
yr_array = np.repeat(bend_yr, 900)

course = [(0,0), midline, (80,0)]
 
sabs = [-5.72957795, -1.19868047, -0.52191351, -0.3039716]
for sab in sabs:
    Car, t = tsim.runSimulation(course, yr_array, 80, sab, 0)		
    print('time to line crossing:', t)

def off_tangent(w=1.5, r = 80, v = 8):
    return np.sqrt(w* (2.0*r +w) / (v**2))

print('formula:', off_tangent())

