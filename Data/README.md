# Data Information and Dictionary

The Data folder contains examples of the raw steering data for one participant. The raw data for gaze is too large to upload an example.

The Midline_40_x and Midline_80_x files are the pool of pre-recorded trajectories that were selected from during the experiment. 

The data filenames are in different forms depending on the experimental block.

Practice files contain 'PRAC'. They will not be discussed.

## Count Task Only
_Distractor task without steering_

Two files:

_Orca18_Distractor_BlockNumber_PPID_EndofTrial.csv_ 

Contains the recorded count estimates. 

Columns are: [], ppid, targetoccurence, targetnumber, trialn, EoTScore1, TargetCount1, EoTScore2, TargetCount2, EoTScore3, TargetCount3

Most of these are self-explanatory. EoTScoreX is the recorded estimates. TargetCountX is the actual amount of times the target was presented.   
Empty cells are common if targetnumber is <3.



_Orca18_Distractor_BlockNumber_PPID_WithinTrial.csv_ 

Contains the button responses within a trial.

Columns: [], ppid, targetoccurence, targetnumber, trialn, CurrentAudio, RT, ResponseCategory, Target1, Target2,	Target3
CurrentAudio is target presented.  
RT is -1 if not responded  
ResponseCategory signifies whether it is a true positive, true negative, false positive, or false negative.  
Empty cells are common if targetnumber is <3.

## Driver Task Only
_Steering Task without Distractor_

One file per trial:

_Orca18_None_BlockNumber_PPID_Radii_TrialN.csv_

## Driver Task with Distractor

One file per trial:

_Orca18_CognitiveDifficulty_PPID_Radii_TrialN.csv_

Due to an error in the filename the distractor task data was over-written when there was a steering task. So we don't have this data.

Each trial file has the following columns:
ppid;	radius;	yawrate_offset;	trialn;	timestamp_exp;	timestamp_trial;	trialtype_signed;	World_x;	World_z;	WorldYaw;	SWA;	YawRate_seconds;	TurnAngle_frames;	Distance_frames;	dt;	WheelCorrection;	SteeringBias;	Closestpt;	AutoFlag;	AutoFile;	OnsetTime;


Wheel correction is the mismatch between virtual yaw-rate and the yaw rate as specified by the steering angle. 




#  Eyetracking Data

Eyetracking data was only collected for the steering task.

There is a separate eye recording, and therefore a separate calibration, for each steering experiment. 

Filenames are:
Orca18_None_BlockNumber_PPID (for steering without cognitive load)

Orca18_CognitiveDifficulty_PPID (for steering + cognitive load). 

## Annotations
"DistractorScreen"	- wait screen at the start of the experiment

'Start_' + self.Trial_SaveName - for the beginning of a trial

'Disengage_' + self.Trial_SaveName - when the automated vehicle was disengaged.

'Distractor_' + self.Trial_SaveName - when the count screen appears

'End_' + self.Trial_SaveName - for Distractor Screen end.
