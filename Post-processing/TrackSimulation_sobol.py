import sys
#rootpath = 'C:\\VENLAB data\\TrackMaker\\'
rootpath = 'C:\\git_repos\\TrackMaker\\'
sys.path.append(rootpath)

import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas as pd

import simTrackMaker

import sobol_seq

class vehicle:
    
    def __init__(self, initialyaw, speed, dt, yawrate_readout, Course):
            
        #self.pos = np.array([pos[0], pos[1]])
        self.yawrate_readout = yawrate_readout 
        self.playback_length = len(yawrate_readout)
        
        self.pos = np.array(Course[0])
        self.yaw = initialyaw #heading angle, radians
        self.speed = speed 
        self.dt = dt       
        self.midline = Course[1]
        self.trackorigin = Course[2]
        
        self.yawrate = 0
        
        self.pos_history = []
        self.yaw_history = []
        self.yawrate_history = []                
        self.error_history = []   
        self.closestpt_history = []         

        self.Course = Course      
        
        self.currenterror, self.closestpt = self.calculatebias()      

                # self.save_history()     
        

    def calculatebias(self):

        #TODO: cut down on processing but only selecting a window of points based on lastmidindex.
        midlinedist = np.sqrt(
            ((self.pos[0]-self.midline[:,0])**2)
            +((self.pos[1]-self.midline[:,1])**2)
            ) #get a 4000 array of distances from the midline
        idx = np.argmin(abs(midlinedist)) #find smallest difference. This is the closest index on the midline.	

        closestpt = self.midline[idx,:] #xy of closest point
        dist = midlinedist[idx] #distance from closest point				

        #Sign bias from assessing if the closest point on midline is closer to the track origin than the driver position. Since the track is an oval, closer = understeering, farther = oversteering.
        middist_from_origin = np.sqrt(
            ((closestpt[0]-self.trackorigin[0])**2)
            +((closestpt[1]-self.trackorigin[1])**2)
            )  #distance of midline to origin
        pos_from_trackorigin = np.sqrt(
            ((self.pos[0]-self.trackorigin[0])**2)
            +((self.pos[1]-self.trackorigin[1])**2)
            ) #distance of driver pos to origin
        distdiff = middist_from_origin - pos_from_trackorigin #if driver distance is greater than closest point distance, steering position should be understeering
        steeringbias = dist * np.sign(distdiff)     

        return steeringbias, closestpt


    def move_vehicle(self, newyawrate):           
        """update the position of the vehicle over timestep dt"""                        
                                 
        self.yawrate = newyawrate

        # self.yawrate = np.deg2rad(0.5) # np.random.normal(0, 0.001)

        maxheadingval = np.deg2rad(35.0) #in rads per second
        
        self.yawrate = np.clip(self.yawrate, -maxheadingval, maxheadingval)
        # print(self.yawrate)
        # self.yawrate = 0.0

        self.yaw = self.yaw + self.yawrate * self.dt  #+ np.random.normal(0, 0.005)
        
        #zrnew = znew*cos(omegaH) + xnew*sin(omegaH);
        #xrnew = xnew*cos(omegaH) - znew*sin(omegaH)

        x_change = self.speed * self.dt * np.sin(self.yaw)
        y_change = self.speed * self.dt * np.cos(self.yaw)
        
        self.pos = self.pos + np.array([x_change, y_change]) 

        self.currenterror, self.closestpt = self.calculatebias()
        
        self.save_history()
    
    def save_history(self):

        self.pos_history.append(self.pos)        
        self.yaw_history.append(self.yaw)
        self.yawrate_history.append(self.yawrate)
        self.error_history.append(self.currenterror)
        self.closestpt_history.append(self.closestpt)

    

def runSimulation(Course, yawrate_readout, myrads, yawrateoffset= 0, onsettime = 0):

    """run simulation and return RMS"""

    #Sim params
    fps = 60.0
    speed = 8.0
 
    yawrateoffset_rads = np.deg2rad(yawrateoffset)
   # print ("speed; ", speed)

    dt = 1.0 / fps
    run_time = 15 #seconds
    time = 0

    Car = vehicle(0.0, speed, dt, yawrate_readout, Course)

    i = 0
    
    crossed = False
    time_til_crossing = None

    f = lambda t: np.exp(-1/t)*(t > 0)
    smooth_step = lambda t: f(t)/(f(t) + f(1 - t))

    print ("playback lenght", Car.playback_length)
    while (time < run_time) and (crossed==False):

        #print i

        time += dt              


        if (i < Car.playback_length):
            newyawrate = np.deg2rad(Car.yawrate_readout[i])
        else:
            #if exceeding playback just put the bend yawrate in.
            newyawrate = speed / myrads


        if time > onsettime:
            time_after_onset = time - onsettime
            transition_duration = .5
            newyawrate += smooth_step(time_after_onset/transition_duration)*yawrateoffset_rads
            
        
        Car.move_vehicle(newyawrate)           
        
        if crossed == False and abs(Car.currenterror) > 1.5:
            time_til_crossing = time - onsettime
            crossed = True

        i += 1

    return Car, time_til_crossing
    
    #RMS = np.sqrt(np.mean(steeringbias**2))

    #print ("RMS: ", RMS)

     

def plotCar(plt, Car):
    """Plot results of simulations"""

    positions = np.array(Car.pos_history)
                        
    steeringbias = np.array(Car.error_history)

    if max(abs(steeringbias)) > 1.5:
        plt.plot(positions[:,0], positions[:,1], 'ro', markersize=.1)						
    else:
        plt.plot(positions[:,0], positions[:,1], 'go', markersize=.1)			            
    
if __name__ == '__main__':
    

    #create straight.
    L = 16#2sec.
    myStraight  = simTrackMaker.lineStraight(startpos = [0,0], length= 16)#, texturefile='strong_edge_soft.bmp')

    #Create Bend
    myrads = 80
    myBend = simTrackMaker.lineBend(startpos = myStraight.RoadEnd, rads = myrads, x_dir = 1, road_width=3.0) 

    #midline and edges
    Course_RoadStart = myStraight.RoadStart
    Course_midline = np.vstack((myStraight.midline, myBend.midline))
    Course_OutsideLine = np.vstack(
        (myStraight.LeftLine, myBend.OutsideLine)
        )
    Course_InsideLine = np.vstack((myStraight.RightLine, myBend.InsideLine))
    Course_CurveOrigin = myBend.CurveOrigin
    
    #Plot Bend
    plt.figure(1)
    plt.plot(Course_midline[:,0], Course_midline[:,1], '--k')
    
    #xlimits = Course_CurveOrigin[0]*2
        
    #plt.xlim([0-5, xlimits+5])
    
    #plt.ylim([-Course_CurveOrigin[0]-5, Course_CurveOrigin[1]*2 + Course_CurveOrigin[0]+5])
    plt.plot(Course_OutsideLine[:,0], Course_OutsideLine[:,1],'-k')
    plt.plot(Course_InsideLine[:,0], Course_InsideLine[:,1],'-k')
    plt.axis('equal')    
    plt.title("Sample Participant, Sobol selection")

    plt.xlim([0-5, 65])
    plt.ylim([20, Course_CurveOrigin[1]*2 + Course_CurveOrigin[0]+5])

    #Temp HACK to store in list while I improve trackmaker.
    Course = [
        Course_RoadStart,
        Course_midline, Course_CurveOrigin,
        Course_OutsideLine, Course_InsideLine
        ]
#

    #list of filenames
    if myrads == 40:
        #filename_list = ["Midline_40_0.csv","Midline_40_1.csv","Midline_40_2.csv","Midline_40_3.csv","Midline_40_4.csv","Midline_40_5.csv"]

        filename_list = ["Midline_40_0.csv"]
		
		
    elif myrads == 80:
        #removing _80_0 because of unusual yawrate changes.
        #removing _80_1 because the balanced portion of the experiment uses it.
        filename_list = ["Midline_80_2.csv","Midline_80_3.csv","Midline_80_4.csv","Midline_80_5.csv"]

        
    else:
        raise Exception('Unrecognised radius')

    Trials = 24
    sobol = sobol_seq.i4_sobol_generate(4, Trials) # 0,1 scale

    #print(sobol_3D)

    #rescale
    onset_sobol = sobol[:,2] * 4 + 5
    autofile_sobol = np.round(sobol[:,1] * 3,0) 

    ttlc_limit = 2
    ttlc_stay = 18

    ttlc_sobol = sobol[:,0] * (ttlc_stay-ttlc_limit) + ttlc_limit
    steer_sobol = sobol[:,2] #flag for understeering or oversteering
    
    #***** retrieve approximations of SAB ******
        
    #filename = "SimResults_onset_6_traj_80_1.csv"
    filename = "simulated_roadcrossing.csv"
    #columns are: yr_offset, file_i, onsettime, time_til_crossing
    balanced_results = np.genfromtxt(filename, delimiter=',')

    sab_sobol = np.ones(Trials)    

    balanced_results_notnan = balanced_results[~np.isnan(balanced_results[:,3])]
    balanced_results_understeer = balanced_results_notnan[balanced_results_notnan[:,0]<= 0]
    balanced_results_oversteer = balanced_results_notnan[balanced_results_notnan[:,0]>= 0]

    for i, ttlc in enumerate(ttlc_sobol):
        print(ttlc)
        steer = steer_sobol[i]
        if steer >= .7:
            sim = balanced_results_oversteer
        else:
            sim = balanced_results_understeer

        diffs = sim[:,3] - ttlc
        idx = np.argmin(abs(diffs)) 
        print(idx)
        sab = sim[idx, 0] #closest sab
        print("sab", sab)
        sab_sobol[i] = sab


    def map_ttlc_to_sab(ttlc, width = 1.5, vel = 8.0):

        sab = width / ((ttlc**2)*  vel)
        return (np.degrees(sab))

    print("2.23", map_ttlc_to_sab(2.23))
    sab_mapped = map_ttlc_to_sab(ttlc_sobol)
    print ("sab_sobol:", sab_sobol[1:5])
    print ("sab_mapped:", sab_mapped[1:5])

    #columns: yr_offset, file_i, onsettime, time_til_crossing
    totalrows = Trials
    
    simResults = np.empty([totalrows,4]) 
    


    #each run has pre-set parameters
    for i, sab in enumerate(sab_sobol):
        
        file_i = int(autofile_sobol[i])
        file = filename_list[file_i]        
        playbackdata = pd.read_csv("Data//"+file) 	
        yawrate_readout = playbackdata.get("YawRate_seconds")

        onset = onset_sobol[i]  
            

        Car, t = runSimulation(Course, yawrate_readout, myrads, sab, onset)
                #Plot Car
        plotCar(plt, Car)

        simResults[i] = [sab, file_i, onset, t]        

                #print(t)
        print ("SAB: ", sab, "Onset: ", onset, "Time til Crossing: ", t)

  #  plt.savefig('sample_participant.png', dpi=800)
    #plt.show()
    
    #np.savetxt("SimResults_OnsetTimes_"+str(myrads)+".csv", simResults, delimiter=",")

   # np.savetxt("SimResults_samplesobol_onsettimes.csv", simResults, delimiter=",")

    #***** plot against distribution of onset_times ******

    
    #balanced_ttlc =  [2.23333333,  4.68333333,  7.1       ,  9.5       , 12.15      ]
    
    balanced_ttlc =  [2.23333333,  4.68333333,  7.1,  9.5]
    balanced_sab = [-5.72957795, -1.19868047, -0.52191351, -0.3039716]
    
    #calculate proportion of takeovers
    def prop_stay(ttlcs, thresh = 9):
        stay_mask = ttlcs[ttlcs>=thresh]
        prop_stay = float(len(stay_mask)) / float(len(ttlcs))
        return prop_stay

    repetitions = 6
    balanced_ttlcs_reps = np.repeat(balanced_ttlc, 	repetitions)
    all_ttlcs = np.concatenate((ttlc_sobol,balanced_ttlcs_reps))
    print("all_ttlcs", all_ttlcs)

    proportion_stay = prop_stay(all_ttlcs)
    print("proportion_stay", proportion_stay)

    stay_sobol = prop_stay(ttlc_sobol)
    print("stay_sobol", stay_sobol)

    stay_balanced = prop_stay(balanced_ttlcs_reps)
    print("stay_balanced", stay_balanced)

    sab_rng = np.linspace(-6, 6, 1000)
    approx_ttlc = np.sqrt(3.0/(np.radians(np.abs(sab_rng))*8))

    plt.figure(2)
    plt.plot(balanced_results[:,0], balanced_results[:,3], 'k.', markersize=5, alpha = .2)
    plt.xlabel("Yaw Rate Offset (deg/s)")
    plt.ylabel("Time from Onset to Lane Crossing (s)")
    plt.plot(sab_sobol, ttlc_sobol, 'r.', markersize = 5)
    plt.plot(balanced_sab, balanced_ttlc, 'b.', markersize = 10)
    plt.title("Sample Participant, Proportion Takeover: " + str(1-proportion_stay))
    plt.axhline(y = 9)
    plt.ylim([0, 20])
    plt.plot(sab_rng, approx_ttlc)
    #plt.savefig('sample_sobol.png', dpi = 300)
    plt.show()



    #plot yr and time til crossing functions.
    # plt.figure(2)
    # plt.plot(simResults[:,0], simResults[:,1], 'k-')
    # plt.xlabel("Yaw Rate Offset (deg/s)")
    # plt.ylabel("Time from Onset to Lane Crossing (s)")
    # plt.title("Radius: " + str(myrads))
    # plt.savefig(str(myrads) + '_Sims_OnsetTimes.png', dpi = 300)
    # #plt.show()
    

    