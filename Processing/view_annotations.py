###script to view and almagamate eye tracking files.
import numpy as np
import sys, os
from file_methods import *
import pickle
import matplotlib.pyplot as plt
import cv2
import csv
import pandas as pd

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

if __name__ == '__main__':
	
    rootdir = "D:/EyeTrike_Backup/Recordings/Orca_Copy/Orca18_Easy_P01/000"
        
    for dirs in os.walk(rootdir):
        path = str(dirs[0]) + "/"
        print (path)
        if os.path.exists(path + 'pupil_data'):
        
            data = load_object(path + "/pupil_data")            

            trial_timestamps = pd.DataFrame(columns = ['trialcode','cogload', 'block', 'ppid', 'radii', 'count','t1','t2','t3']) 
            notes = data["notifications"]
            entry = 0
           
            for n in notes:
                if n['subject'] == 'annotation':
                    #collect trial type and timestamp        
                #### here you can crunch a trial to make sure that the correct timestamps are shown. 
                    t = n['timestamp']
                    print ("timestamp:", t)
                    label = n['label']
                    print ("label:", label)


                    #find reason
                    i1 = label.find('_')  #finds first underscore.
                    if i1 != -1: #if not the distractor screen 
                        reason = label[:i1]  
                        trialcode = label[i1+1:]

                        print ("trialcode:", trialcode)
                        print ("reason:", reason)                        
                        if ("Easy" in label) or ("Hard" in label):

                            exp_id, cogload, ppid, radii, count = trialcode.split('_')
                            block = 0
                        else:
                            exp_id, cogload, block, ppid, radii, count = trialcode.split('_')
                    
                        # check if the trialcode is already in the database. I want one row per trialcode. 
                        mask =  trial_timestamps[trial_timestamps['trialcode'].str.contains(trialcode)]
                        print("Mask: ", mask)
                        print("MaskEmpty: ", mask.empty)
                        if not mask.empty:
                            # if entry already exists, add the timestamp to the relevant column.
                            idx = trial_timestamps.index[trial_timestamps['trialcode']==trialcode]
                                
                            if "Sta" in reason:
                                trial_timestamps.loc[idx,'t1'] = t
                            elif "Dis" in reason:
                                trial_timestamps.loc[idx,'t2'] = t                            
                            elif "End" in reason:
                                trial_timestamps.loc[idx,'t3'] = t
                            
                            
                            
                        else: 
                            
                            # create new entry.                                         
                            if "Sta" in reason:
                                t1 = t
                                t2 = 0
                                t3 = 0
                            elif "Dis" in reason:
                                t1 = 0
                                t2 = t
                                t3 = 0
                            elif "End" in reason:
                                t1 = 0
                                t2 = 0
                                t3 = t
                                    
                            row = [trialcode, cogload, block, ppid, radii, count, t1, t2, t3]    
                            trial_timestamps.loc[entry,:] = row
                            entry += 1    
            print (trial_timestamps)

