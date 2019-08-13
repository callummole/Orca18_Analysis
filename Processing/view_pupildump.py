#file to play around with loading the dump.
    
import gzip
import json
import numpy as np
import os
import msgpack
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    
    #rootdir = "E:\\Trout_rerun_pupildumps"    
    rootdir = "E:\\Orca19_3pp_venlab\\P03"
        
    #fn = rootdir + "\\Trout_203_2.pupil.jsons.gz"
    fn = rootdir + "\\Orca19_None_Calibration_3.pupil.jsons.gz"   


    #json.dumps((topic, timestamper(), msgpack.loads(msg)))
    pupils = [],[]

    gaze_data = []
    timestamps = []

    """
    there is an error on one line of the pupil_log. Try skipping using the below code...

    """
    """

    """
    pupil_log = map(json.loads, gzip.open(fn))

    #for topic, *_ in pupil_log:
    #    if topic == "notify.recording.started":
    #        break
    
    lines = iter(pupil_log)
    #err = 0
    
    while lines:
    #for line in lines:
        
        try:   
            line = next(lines)            
            #print(line)
            topic, ts, data = line
            
            timestamps.append(ts)
            data['recv_ts'] = ts
            if topic == "pupil.0":
                pupils[0].append(data)
            if topic == "pupil.1":
                pupils[1].append(data)      
            
            gaze_data.append(data)
        except StopIteration:
            break
        except Exception as e:
            #err += 1
            print (e)
            pass
    
        
    print("length: ", (timestamps[-1] - timestamps[0]) / 60)

    num = 2000 #this will take ages.
    plt.figure(2)

    pupils = pupils[0][:num], pupils[1][:num]

    
    for pupil in pupils:
        for p in pupil:
            
            npos= p['norm_pos']
            t = p['recv_ts']

            #plt.plot(npos[0],npos[1], 'b.', alpha = .4)
            plt.plot(t,npos[0], 'b.', alpha = .4)
        break #just do the one pupil
     
    plt.show()