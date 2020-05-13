import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.patches import Polygon
import drivinglab_projection as dp
import pandas as pd 
import TrackSimulation as tsim
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

track = pd.read_csv("..//Data//track_with_edges_orca_80.csv")

# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"

def remove_straight(cols):

	road = np.array([track[cols[0]].values, track[cols[1]].values]).T
	road[:,1] = road[:,1] - 16 #remove straight.
	mask = road[:,1]>= 0
	road = road[mask,:]	
	return road

inside_edge = remove_straight(['insidex','insidez'])
outside_edge = remove_straight(['outsidex','outsidez'])
midline =  remove_straight(['midlinex','midlinez'])


angles_limits_bottom = dp.screen_to_angles([0,0])[0]

#print(angles_limits_bottom)
angles_limits_top = dp.screen_to_angles([1,1])[0]

pixels_limits_top = [1920,1080]

viewpoint = (0, 0) #end of straight, start of curve.
yaw = 0      

#compute track from viewpoint.
inside_edge_angles,depth = dp.world_to_angles_through_screen(inside_edge, viewpoint, yaw)
outside_edge_angles,depth = dp.world_to_angles_through_screen(outside_edge, viewpoint, yaw)    

#remove any above the horizon
inside_edge_angles = inside_edge_angles[depth > 0,:]
outside_edge_angles = outside_edge_angles[depth > 0,:]

#calculate steering bias so it is smooth and not prey to discretisatiokn
def smooth_steering_bias():

	bend_yr = 8.0/80
	print(bend_yr)
	print(np.rad2deg(bend_yr))
	#yr_bias = np.deg2rad(-5) #-5 degrees per second.	
	yr_bias = -bend_yr

	speed = 8.0
	dt = 1.0/60
	current_yr_bias = 0
	time = 0
	run_time = 2
	pos_hist = []
	current_sb = 0
	timestep = []

	while (time < run_time):

		#do loop
		time += dt
		current_yr_bias += yr_bias*dt #each time step the yr diff increases

		sb_change = speed * dt * np.sin(current_yr_bias)
		current_sb += sb_change
		pos_hist.append(current_sb)
		timestep.append(time)


	fig, ax = plt.subplots()
	#ax.axhline(y = -1.5, color = (.5, .5, .5), linestyle = '--')
	#ax.axhline(y = 1.5, color = (.5, .5, .5), linestyle = '--')
	#ax.axhline(y = 0, color = (.9, .9, .9), linestyle = '--')
	
	pos_hist = np.array(pos_hist)
	ax.plot(timestep, pos_hist)
	plt.show()

	


def plot_road(ax):

	verts = np.vstack([outside_edge_angles, np.flipud(inside_edge_angles)])
	

	midline_angles,depth = dp.world_to_angles_through_screen(midline, viewpoint, yaw)    	
	midline_angles = midline_angles[depth > 0,:]
	ax.plot(midline_angles[:,0], midline_angles[:,1], 'k-', alpha = .6)    

	poly = Polygon(verts, facecolor='0.9', edgecolor='0.9')
	ax.add_patch(poly)

	ax.set_ylim(angles_limits_bottom[1],angles_limits_top[1])
	ax.set_xlim(angles_limits_bottom[0],angles_limits_top[0])

	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xlabel("Horizontal Angle ($^\circ$)", fontsize = 10)
	ax.set_ylabel("Vertical Angle ($^\circ$)", fontsize = 10)


	return ax

def plot_trajectories(Cars, ax, col = 'xkcd:blue'):

	for car in Cars:

		#compute pos from viewpoint.
		positions = np.array(car.pos_history)
		angles,depth = dp.world_to_angles_through_screen(positions, viewpoint, yaw)
		
		#remove any above the horizon
		angles = angles[depth > 0,:]

		ax.plot(angles[:,0], angles[:,1], color = col, alpha = .6)

	return ax

def simulate_offsets(course, yr_array, onset_time, sabs, smooth = True):
	
	Cars = []
	for yr in sabs:
		Car, _ = tsim.runSimulation(course, yr_array, 80, yr, onset_time, smooth)		
		Cars.append(Car)
	return Cars
				

def simulate_sobol_offsets():

	filename = "SimResults_samplesobol_onsettimes.csv"
	#columns are: yr_offset, file_i, onsettime, predicted_time_til_crossing
	#sobol_condition_list = np.genfromtxt(filename, delimiter=',')

	sobol_condition_list = pd.read_csv(filename, 
					sep=',', 
					names=["sab", "autofile_i", "onsettime", "simulated_ttlc"])


#smooth_steering_bias()

bend_yr = np.rad2deg(8.0 / 80.0)
bend_yr_array = np.repeat(bend_yr, 900)
#yawrate_readout[0] = 0


filename = "Midline_80_1.csv"
playbackdata = pd.read_csv("..//Data//"+filename) 	
yawrate_readout = playbackdata.get("YawRate_seconds").values
#yawrate_readout = yawrate_readout[mask]
world_z = playbackdata.get("World_z").values
yawrate_readout = yawrate_readout[world_z>= 16]
#midline = np.array([track['midlinex'].values, track['midlinez'].values]).T

course = [(0,0), midline, (80,0)]

Cars_ideal_bal = simulate_offsets(course, bend_yr_array, .5, sabs = [-5.72957795, -1.19868047, -0.52191351, -0.3039716], smooth = False)

filename = "..//Data//SimResults_samplesobol_onsettimes.csv"
#columns are: yr_offset, file_i, onsettime, predicted_time_til_crossing	
sobol_condition_list = pd.read_csv(filename, sep=',', names=["sab", "autofile_i", "onsettime", "simulated_ttlc"])
sabs = sobol_condition_list.get('sab').values
Cars_ideal_sobol = simulate_offsets(course, bend_yr_array, .5, sabs = sabs, smooth = False)


cmap = cm.get_cmap('tab10')
rgbas = cmap([6,1,2,4])


###plot idealised####
fig_cm = np.array([13.2,10])
fig_inc = fig_cm /2.54 

fig, axes = plt.subplots(2,2, figsize = fig_inc, constrained_layout = True)

#plot lane position

axes[0,0].axhline(y = -1.5, color = (.7, .7, .7), linestyle = '--')
axes[0,0].axhline(y = 1.5, color = (.7, .7, .7), linestyle = '--')


for car in Cars_ideal_sobol:
	sb = car.error_history
	axes[0,0].plot(np.arange(0, len(sb)) / 60 - .5, sb, color = 'xkcd:light grey')

	axes[0,0].plot((len(sb)-1)/60 - .5, sb[-1], color = 'xkcd:light grey', marker = ".", markersize = 5)

	yr = car.yawrate_history
	axes[0,1].plot(np.arange(0, len(yr)) / 60 - .5, np.rad2deg(yr), color = 'xkcd:light grey')

	axes[0,1].plot((len(yr)-1)/60 - .5, np.rad2deg(yr[-1]), color = 'xkcd:light grey', marker = ".", markersize = 5)


for car, col in zip(Cars_ideal_bal, rgbas):
	sb = car.error_history
	axes[0,0].plot(np.arange(0, len(sb)) / 60 - .5, sb, color = col)

	axes[0,0].plot((len(sb)-1)/60 - .5, sb[-1], color = col, marker = ".", markersize = 5)

	yr = car.yawrate_history
	axes[0,1].plot(np.arange(0, len(yr)) / 60 - .5, np.rad2deg(yr), color = col)

	axes[0,1].plot((len(yr)-1)/60 - .5, np.rad2deg(yr[-1]), color = col, marker = ".", markersize = 5)


axes[0,0].set_xlabel('Time after Failure Onset (s)', fontsize = 10)
axes[0,0].set_ylabel('Lane Position (m)', fontsize = 10)
axes[0,0].invert_yaxis()
axes[0,0].set_yticks([-1.5, 0, 1.5]) 
axes[0,0].set_ylim([1.7,-1.7])  

axes[0,1].axhline(y =np.rad2deg(8.0 / 80), color = (.2, .2, .2), linestyle = ':')
axes[0,1].set_xlabel('Time after Failure Onset (s)', fontsize = 10)
axes[0,1].set_ylabel('Yaw rate ($^\circ$/s)', fontsize = 10)

####Including autofile and onset time noise.
sobol_list = ["Midline_80_2.csv","Midline_80_3.csv","Midline_80_4.csv","Midline_80_5.csv"]
sobol_YRs = []
for file in sobol_list:
	playbackdata = pd.read_csv("..//Data//"+file) 	
	readout = playbackdata.get("YawRate_seconds").values
	#yawrate_readout = yawrate_readout[mask]
	world_z = playbackdata.get("World_z").values
	readout = readout[world_z>= 16]
	sobol_YRs.append(readout)

onsets = sobol_condition_list.get('onsettime').values
autofile_idxs = sobol_condition_list.get('autofile_i').values
Cars_sobol = []

axes[1,0].axhline(y = -1.5, color = (.7, .7, .7), linestyle = '--')
axes[1,0].axhline(y = 1.5, color = (.7, .7, .7), linestyle = '--')


for sab, ons, afile in zip(sabs, onsets, autofile_idxs):
	yr_readout = sobol_YRs[int(afile)]
	Car, _ = tsim.runSimulation(course, yr_readout, 80, sab, ons-2)		
	Cars_sobol.append(Car)

for car in Cars_sobol:
	sb = car.error_history
	axes[1,0].plot((np.arange(0, len(sb)) / 60) + 2, sb, color = 'xkcd:light grey')

	axes[1,0].plot((len(sb)-1)/60 + 2, sb[-1], color = 'xkcd:light grey', marker = ".", markersize = 5)


	yr = car.yawrate_history

	axes[1,1].plot((np.arange(0, len(yr)) / 60) + 2, np.rad2deg(yr), color = 'xkcd:light grey')	

	axes[1,1].plot((len(yr)-1)/60 + 2, np.rad2deg(yr[-1]), color = 'xkcd:light grey', marker = ".", markersize = 5)


Cars_bal = simulate_offsets(course, yawrate_readout, 4, sabs = [-5.72957795, -1.19868047, -0.52191351, -0.3039716])

for car, col in zip(Cars_bal, rgbas):
	sb = car.error_history
	axes[1,0].plot((np.arange(0, len(sb)) / 60) + 2, sb, color = col)

	axes[1,0].plot((len(sb)-1)/60 + 2, sb[-1], color = col, marker = ".", markersize = 5)

	yr = car.yawrate_history
	axes[1,1].plot( (np.arange(0, len(yr)) / 60) + 2, np.rad2deg(yr), color = col)
	
	axes[1,1].plot((len(yr)-1)/60 + 2, np.rad2deg(yr[-1]), color = col, marker = ".", markersize = 5)

axes[1,0].set_xlabel('Time into Trial (s)', fontsize = 10)
axes[1,0].set_ylabel('Lane Position (m)', fontsize = 10)
axes[1,0].axvline(x = 15, color = 'k', linestyle = '-')
axes[1,0].set_xlim([2, 17])
axes[1,0].invert_yaxis()
axes[1,0].set_yticks([-1.5, 0, 1.5])  
axes[1,0].set_ylim([1.7,-1.7])  

axes[1,1].axhline(y =np.rad2deg(8.0 / 80), color = (.1, .1, .1), linestyle = ':')
axes[1,1].set_xlabel('Time into Trial (s)', fontsize = 10)
axes[1,1].set_ylabel('Yaw rate ($^\circ$/s)', fontsize = 10)
axes[1,1].axvline(x = 15, color = 'k', linestyle = '-')
axes[1,1].set_xlim([2, 17])

for a, label in zip(axes.flat, ['A','B','C','D']):
		a.text(0.05, 1.02, label, transform=a.transAxes,
      		fontsize=12, fontweight='bold', va='top')
		a.spines['right'].set_visible(False)
		a.spines['top'].set_visible(False)

for ax in axes.flat:
	for t in ax.get_xticklabels(): t.set_fontsize(8)		
	for t in ax.get_yticklabels(): t.set_fontsize(8)


legend_elements = [Line2D([0],[0],color = rgbas[0], lw = 2, label = "1.93 s"),
					Line2D([0],[0],color = rgbas[1], lw = 2, label = "4.21 s"),
					Line2D([0],[0],color = rgbas[2], lw = 2, label = "6.39 s"),
					Line2D([0],[0],color = rgbas[3], lw = 2, label = "8.37 s"),
				Line2D([0],[0], color = "xkcd:light grey", lw = 2, label = "Variable")] 


axes[0,1].legend(handles = legend_elements, loc = [.6,.02], fontsize = 8, frameon = False, title = r'$TTLC_F$', title_fontsize = 8, labelspacing = .4)
	

#annotations
axes[0,1].annotate('Bend YR', xy=(17, np.rad2deg(8.0 / 80)), xytext=(15, np.rad2deg(8.0 / 80)+1.5), fontsize = 8,
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

axes[1,1].annotate('Trial\nEnd', xy=(15, np.rad2deg(8.0 / 80)-4.5), xytext=(12, np.rad2deg(8.0 / 80)-6), fontsize = 8,
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

#axes[0,0].annotate("Understeering", xy=(15, -1.4), fontsize = 8)
#axes[0,0].annotate("Oversteering", xy=(15, 1.4), fontsize = 8)

#axes[1,0].annotate('Onset\nfor Repeated', xy=(6, .1), xytext=(3, 1), fontsize = 8,
#            arrowprops=dict(arrowstyle="->",
#                            connectionstyle="arc3"),
#            )

plt.savefig('failures_schematic.png', format='png', dpi=300, bbox_inches = "tight", facecolor=plt.gcf().get_facecolor(), edgecolor='none')
plt.savefig('failures_schematic.svg', format='svg', dpi=300, bbox_inches = "tight", facecolor=plt.gcf().get_facecolor(), edgecolor='none')
plt.savefig('failures_schematic.eps', format='eps', dpi=300, bbox_inches = "tight", facecolor=plt.gcf().get_facecolor(), edgecolor='none')
#	plt.savefig('Fits/linmix_sample_' + str(fit['ID']) +'_' + str(fit['drivingmode']) +'.svg', format='svg', dpi=300, bbox_inches = "tight", facecolor=plt.gcf().get_facecolor(), edgecolor='none')
plt.show()








"""
fig, ax = plt.subplots()
ax.plot(midline[:,0],midline[:,1])
for car in Cars:
	positions = np.array(car.pos_history)
	ax.plot(positions[:,0],positions[:,1])
plt.show()
"""

#in perspective. Doesn't really work as you cannot tell apart the faiures
"""
fig1, ax1 = plt.subplots()
ax1 = plot_road(ax1)
ax1 = plot_trajectories(Cars, ax1)
plt.show()
"""


"""
#steering bias
#can figure out how to do this in a smooth way. 
fig2, ax2 = plt.subplots()
ax2.axhline(y = -1.5, color = (.5, .5, .5), linestyle = '--')
ax2.axhline(y = 1.5, color = (.5, .5, .5), linestyle = '--')

for car in Cars:
	sb = car.error_history
	ax2.plot(np.arange(0, len(sb)) / 60, sb)
plt.show()

#yaw rate
fig3, ax3 = plt.subplots()
ax3.axhline(y =(8.0 / 80), color = (.5, .5, .5), linestyle = '--')
#ax3.axhline(y = 1.5, color = (.5, .5, .5), linestyle = '--')

for car in Cars:
	yr = car.yawrate_history
	ax3.plot(np.arange(0, len(yr)) / 60, yr)
plt.show()
"""


"""
plt.plot(inside_edge[0],inside_edge[1])
plt.plot(outside_edge[0],outside_edge[1])
for c in Cars:
	plt.plot(c[:,0],c[:,1])

plt.show()
"""
#need to simulate the trajectories, find the viewpoint that corresponds with one autofile.
#sobol sequence.
