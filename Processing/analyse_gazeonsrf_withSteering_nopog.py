###script to view and almagamate eye tracking files.
import numpy as np
import sys, os
from file_methods import *
import pickle
import matplotlib.pyplot as plt
import cv2
import csv
import pandas as pd
import math as mt 
from scipy.interpolate import interp1d
from nslr_hmm import *
from timeit import default_timer as timer

"""
Useful docs:
    
To understand marker transform:    
https://www.pyimagesearch.com/2014/08/25/4-point-opencv-getperspective-transform-example/

Marker centroids & vertices appear to be given in the camera's resolution.

Steps:
    1) Plot markers and gaze position per frame to explore data.
    2) Save in csv file with experiment and trial ID.
    3) Use matlab script -- converted into python -- to generate measures, such as angles etc. (or could do that in this script)
    
    #Screen co-ords.
    # +-----------+
    # |0,1     1,1|  ^
    # |           | / \
    # |           |  |  UP
    # |0,0     1,0|  |
    # +-----------+
    
Now I have the markers positions I need to do more with them than simple translate them into screen co-ords.
I also need to correct for tilt of the camera, and correct gaze position in relation to this. 
Look at how their screen tracker and surface transformations do this. 
    
"""


def correlateTimestamps(df):
   
    #TODO: use the eye_id to average over duplicate timestamps with multiple eyes. 
    #FIXME: This function is slow and doesn't work properly. Dealing with dataframes in this way seems slow. 
    
    #receives panda dataframe and crunches the multiple timestamps.
    world_ts = np.unique(df['world_timestamp'])
    data_by_ts = pd.DataFrame()    
        
    
    def avgAcrossTs(stamps):
        #receives a list of entries associated with that timestamp.
        #take the most confidence. If all equal, average.
        print(stamps)
        mx = stamps['confidence'].max() 
        gp= stamps.loc[stamps['confidence']==mx]
        avggp = gp.mean(0)        
        return avggp
    
    #loop through surface_timestamps and find closest frame index 
    
    for row in range(df.shape[0]):
        #print(row)
        datum = df.iloc[[row]] #needs list selection for dataframe.   
        surf_ts = df.iloc[row,1] #surface timestamp        
        if surf_ts >0 and surf_ts<20:    #check ts is during trial. If not don't process datum.
            tsdiff = world_ts - surf_ts 
            idx = np.argmin(abs(tsdiff))
            datum['frame_index'] = idx #Add new column.
            data_by_ts = pd.concat([data_by_ts,datum])
                
    #now you've correlated timestamps. Return average for duplicate entries. 
    #
    df_new = data_by_ts.groupby(['frame_index']).apply(avgAcrossTs)
   
#    for frame in range(len(world_ts)):
#        
#        #find all the entries.
#        ts_entries = []
#        for i in range(len(data_by_ts)):
#            dt = data_by_ts[i]
#            if dt['frame_index'] == frame:
#                ts_entries.append(i)
#                
#        if len(ts_entries) > 0:
#            Tlist = [data_by_ts[ts_entries[i]] for i in range(len(ts_entries))]    
#            Tdf = pd.DataFrame(Tlist)          
#            d= avgAcrossTs(Tdf)
#            df_new.append(d)            
    
    return df_new
    
def GazeAngles(df):
    
        
    def SurfaceToGazeAngle(gp):
        
        #proj_width = 1.965
        #proj_height = 1.115    
        
        #pixel size of markers of total white border marker is 118 x 107. But I think surface is up to black marker edge.
        #measurements of black marker edge in inkscape are ~75 x 72 pixels. 
        #NEED REAL_WORLD MEASUREMENTS OF SURFACE SIZE UP TO BLACK SQUARE.
        #AND HORIZON RELATIVE TO BOTTOM AND TOP OF SQUARE.
        #put the below into a function and call apply: 

        #need to check that surfaces is from edges, as in vizard the placement is the bottom left corner
        width = 1.656 #measured at 165.6 cm on 14/12/18 #real-world size of surface, in m.
        height = .634 #measured at 63.4 cm on 18/12/18
        #this should match the defined surface.
        
        
        centrex = .5 #horiz centre of surface is centre of experiment display.
        
        Horizon_relativeToSurfaceBottom = .455 #Horizon measured at 45.5 above bottom marker value . 45.5/63.5 = .7063

        #it is very important centrey is accurate. make sure you measure up to the true horizon, not the clipped horizon because this was overestimate how far people are looking
        #Measured at 46cm above the bottom marker. 46/60 is .7666667

        
        centrey = Horizon_relativeToSurfaceBottom / height #.7667 #placement of horizon in normalised surface units. Minus off norm_y to get gaze relative to horizon.


        #TODO: CHECK HORIZON MEASUREMENT AS SURFACE DOESN'T ALIGN WITH TOP OF MARKERS. NEED TO CORRECT FOR DISTORTION.
        screen_dist = 1.0 #in metres
        
        #convert the scale to real-distances from centre.
        x = gp['x_norm']
        y = gp['y_norm']
        real_h = (x-centrex)*width
        real_v = (y-centrey)*height
#	
#	
    	#calculate gaze angle
        hrad = mt.atan(real_h/screen_dist)
        vrad = mt.atan(real_v/screen_dist)
#	
    #	#convert to degrees
        hang = (hrad*180)/mt.pi
        vang= (vrad*180)/mt.pi
#	
        return (hang, vang) 
    
    
    df['hangle'], df['vangle'] = zip(*df.apply(SurfaceToGazeAngle,axis=1))
       
    return df	
    
def LoadSteering(path, trial, maxt, mint):
    #return steering data for trial. 
   # df_steer = pd.DataFrame() #trial frame by frame data
                        
    filename = path + str(trial) + '.csv' #formula for reading in steering file. 
    print("Steering Filename: ", filename)
    untrimmed_df_steer= pd.read_csv(filename) #read steering data                       
    #df_steer= pd.read_csv(filename) #read steering data                       
#    endtrial = len(importMat) #0 to 1 signifies end of trial.                    
#    trialdata = importMat.iloc[1:endtrial-1,:] #take all of the data, minus the first and last frame.
    
    #print(untrimmed_df_steer) #Not sure why I trimmed this.

    df_steer = untrimmed_df_steer.loc[:,['ppid','radius','yawrate_offset','trialn','timestamp_exp',
    'timestamp_trial','trialtype_signed','World_x','World_z','WorldYaw','SWA','YawRate_seconds','TurnAngle_frames','Distance_frames',
    'dt','WheelCorrection','SteeringBias','Closestpt','AutoFlag','AutoFile',
    'OnsetTime']].copy()
    
    #Since stitching together gaze and steering data relies on creating an interpolation function for the gaze data, then passing it the steering timestamps, we need the steering timestamps to be within the range of values for gaze timestamps.
        
    lower = df_steer['timestamp_exp']>mint
    upper = df_steer['timestamp_exp']<maxt
    
    df_steer = df_steer.loc[lower&upper, :].copy()
        
    return(df_steer)
    
def StitchGazeAndSteering(df, df_steer):
    
    #using gaze angles, linearly interpolate simulator angles. Frame rate is high enough that using the 'fake' gaze data isn't a problem.               
    #use all gaze angles to interpolate
    #This function relies on timestamps of simulator and pupil-labs to be synced at start of each trial.
    yinterpolater = interp1d(df['gaze_timestamp'].values,df['vangle'].values)
    y = yinterpolater(df_steer['timestamp_exp'].values)
                    
    xinterpolater = interp1d(df['gaze_timestamp'].values,df['hangle'].values)
    x = xinterpolater(df_steer['timestamp_exp'].values)

    df_steer['vangle'] = y
    df_steer['hangle'] = x    


    #Also add normed position on surface.
    yinterpolater_ynorm = interp1d(df['gaze_timestamp'].values,df['y_norm'].values)
    y_norm = yinterpolater_ynorm(df_steer['timestamp_exp'].values)
                    
    xinterpolater_xnorm = interp1d(df['gaze_timestamp'].values,df['x_norm'].values)
    x_norm = xinterpolater_xnorm(df_steer['timestamp_exp'].values)

    df_steer['y_norm'] = y_norm
    df_steer['x_norm'] = x_norm  

    # #interpolate confidence values.
    # confidenceinterpolater = interp1d(df['gaze_timestamp'].values,df['confidence'].values)
    # confidence = confidenceinterpolater(df_steer['currtime'].values)

    # df_steer['confidence'] = confidence

    return(df_steer)    
  
if __name__ == '__main__':
	#rootdir = sys.argv[1] 
    rootdir = "D:\\EyeTrike_Backup\\Recordings\\Orca_Copy"    
    savedir = "C:\\Users\\psccmo\\Orca18_Analysis\\"
    resave = False #boolean whether to move files to savedir
    
    steerdir = "E:\\ORCA_DATA\\Data\\"

    exp = 'Orca'
    
    master_stitch = pd.DataFrame() #master data for gaze and steering            

    gazedata_filename = '\\gaze_on_surface_Corrected.csv'
    #print(pfolder)
    #print ("here")
    for dirs in os.walk(rootdir):
        path = str(dirs[0]) 
        print (path)
        
        if os.path.exists(path +gazedata_filename):                    
            
            """                    
            ORCA PARTICULARS: There is a trial_timestamps file that contains the start and finish times (on the eyetrike)
            of the individual trials. Use this to partial out the main gaze and steering csv files.

            - For each row in eyetrack_timestamp, select the corresponding data within the start and finish range from the large gaze_df.
            - Then add the steering. Then stitch together.
            
            file has columns: 'world_timestamp', 'surface_timestamp', 'x_norm', 'y_norm', 'on_srf', 'confidence'
            There may be more than one estimate for each timestamp, and they may not be in order.
            Where they is more than one, either: Take the higher conf, or take the avg if conf is equal.         
            """                                                               
            
            begin = timer()

            allgaze_df = pd.read_csv(path+gazedata_filename, sep=',',header=0)                                        
            
            allgaze_df = GazeAngles(allgaze_df) #adds two new columns: hangle and vangle                     
            
            eyetrike_timestamps_df = pd.read_csv(path+'/Eyetrike_trialtimestamps.csv')
            
            ppid = eyetrike_timestamps_df['ppid'][0]
            
            donotprocess = ['P09']
            if ppid not in donotprocess:
            
                for index, row in eyetrike_timestamps_df.iterrows():
                                    
                    #t1 is start of trial.
                    #t2 is disengagement.
                    #t3 is end of trial.                     
                    mint = row['t1']
                    maxt = row['t3']
                
                    #print ("Mint", mint)
                    #print ("Maxt", maxt)                
                    trial = row['trialcode']
                    
                    print("Processing: ", trial)    
                    
                    #hack for pilot
                    #Use Eyetracking timestamps file to select the right portion of the eyetracking datafile.                 
                    lower = allgaze_df['gaze_timestamp'] > mint 
                    upper = allgaze_df['gaze_timestamp'] < maxt
                    trial_gazedata = allgaze_df.loc[lower & upper, :].copy() 
                    
                    if len(trial_gazedata) == 0:
                        print ("ZERO LENGTH:", trial)
                    else:

                        print("Gaze Data length: ", len(trial_gazedata))

                        #carry over section order.
                                            
                        #### Load Steering Timestamps of Same Trial ###########                        
                        #recalculate max and min so that gaze and steering data span same range for interpolation
                        mint = min(trial_gazedata['gaze_timestamp'].values)
                        maxt = max(trial_gazedata['gaze_timestamp'].values)
                        
                        df_steer = LoadSteering(steerdir,trial,maxt, mint) #loads steering file into a df.
                    #   
                        print("steering data length: ", len(df_steer))
                    
                        df_stitch = StitchGazeAndSteering(trial_gazedata,df_steer) #interpolates gaze with simulator timestamps.                        
                                        
                        print (list(df_stitch.columns.values))
                        
                        #Add trial identifiers to df_stitch.                
                        df_stitch['trialcode'] = trial
                        df_stitch['count'] = row['count'] 
                        df_stitch['cogload'] = row['cogload']        
                        df_stitch['ppid'] = row['ppid']
                        df_stitch['radii'] = row['radii']
                        df_stitch['block'] = row['block']     

                        print (list(df_stitch.columns.values))                         
                        
                        print ("added to master df")
                    
                        master_stitch = pd.concat([master_stitch,df_stitch])
            compute_time = timer()-begin
            print("gaze took %f seconds" % compute_time)
                
    
    #now you've built the master data of all trials, save it.
    master_stitch.to_csv(savedir + "GazeAndSteering_longFormat_080419.csv")



