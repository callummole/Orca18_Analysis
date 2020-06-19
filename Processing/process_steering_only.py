###script to view and almagamate eye tracking files.

import os
import pandas as pd
from glob import glob


if __name__ == '__main__':
	#rootdir = sys.argv[1] 
    
    #rootdir = "C:\\VENLAB data\\SP_18-19\\Data\\Orca19_Steering_Only"
    #rootdir = "C:\\VENLAB data\\SP_18-19\\Data\\Orca19_Steering_Only"
    #rootdir = "C:/git_repos/Orca18_Analysis/Data"
    rootdir = "D:/Orca19_FullSteeringDataset"
    output_filename = 'Orca19_collated_steering2.csv'

    

    master_stitch = pd.DataFrame() #master data for gaze and steering         

    glob_match1 = "Orca19_[MN]*[0-9]*_[0-9]*.csv" 
    

    for dirs in os.walk(rootdir): #does it for all dirs in that folder
        path = str(dirs[0]) + "/"
        print (path)
        for fn in glob(path + glob_match1):
            
            print(fn)
            
            if "old" in fn: continue
            trial_data = pd.read_csv(fn)
            #print(trial_data)
            
            trial_data["filename"] = fn
                
            master_stitch = pd.concat([master_stitch,trial_data])
            
        #now you've built the master data of all trials, save it.
    master_stitch.to_csv(rootdir + "/" + output_filename)



