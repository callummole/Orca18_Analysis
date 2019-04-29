library("tidyverse")
library("brms")
library(ggplot2)
#library(brmstools)
library(cowplot)
#library(dplyr)
#library(tidyr)

raw_eyedata = read_csv("Data/EyetrackingData/FramebyFrameData_LongFormat.csv")


#EYE HEADERS: # world_timestamp, world_frame_index, gaze_timestamp, ...
#norm_x, norm_y, on_srf, trialtype, ppcode, pp_id, exp_id, trialn, ...
#Bend, condition, hangle, vangle

eyedata = raw_eyedata
#Step1: Summarise Eye Data into individual trials
###mirror left bends hangle so all bends are on the same side.
eyedata[eyedata$Bend == "Left", "hangle"] <- eyedata[eyedata$Bend == "Left", "hangle"]*-1 

#remove  participant 1-4.
eyedata <- eyedata %>% filter(pp_id %in% c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

eyedata <- eyedata %>% mutate(Auto = ifelse(condition=='Free_Fix' | condition=='Free_Free','Free',
                              ifelse(condition=='Fix_Fix' | condition=='Fix_Free', 'Fix', NA)))
eyedata <- eyedata %>% mutate(Manual = ifelse(condition=='Fix_Free' | condition=='Free_Free','Free',
                              ifelse(condition=='Free_Fix' | condition=='Fix_Fix', 'Fix', NA)))


#Split Data by timestamp to get manual vs automation period
eyemanual <- eyedata %>% filter(world_timestamp > 10)
eyeauto <- eyedata %>% filter(world_timestamp < 10)

#use summarise to get trial measures and number of datapoints.
eyebytrial_manual <- eyemanual %>% group_by(pp_id, exp_id,trialn) %>% summarize(AvgHang = mean(hangle),
                                                                           StdHang = sd(hangle),
                                                                           MedHang = median(hangle),
                                                                           AvgVang = mean(vangle),
                                                                           StdVang = sd(vangle),
                                                                           MedVang = median(vangle),
                                                                           Datapoints = n(),
                                                                           IQRHang = IQR(hangle),
                                                                           IQRVang = IQR(vangle),
                                                                       nframes = n_distinct(world_frame_idx),
                                                                       Auto = Auto[1],
                                                                       Manual = Manual[1],
                                                                       cndt = condition[1])

#filter the eyedata by low nframes (at least 10s) and high std.
#eyebytrial_filtered <- eyebytrial %>% filter(nframes > 300, StdVang <30, StdHang<30)
#eyebytrial_unfiltered <- eyebytrial %>% filter(nframes < 300, StdVang >30, StdHang>30)
#eyebytrial_frames <- eyebytrial %>% filter(nframes < 300)
#eyebytrial_variable <- eyebytrial %>% filter(StdVang > 30, StdHang > 30)

###A VERY CONSERVATIVE FILTER OF EXTREME TRIALS###
eyebytrial_manual_filtered <- eyebytrial_manual %>% filter(nframes > 200, IQRHang < 20)



#now average by condition. Keep the number of trials going into each estimate.
eyebycondition_manual <- eyebytrial_manual_filtered %>% group_by(pp_id, exp_id, Auto, Manual) %>% summarize(AvgHang = mean(AvgHang),
                                                                            StdHang = mean(StdHang),
                                                                            MedHang = mean(MedHang),
                                                                            AvgVang = mean(AvgVang),
                                                                            StdVang = mean(StdVang),
                                                                            MedVang = mean(MedVang),
                                                                            Trialn = n(),
                                                                            IQRHang = mean(IQRHang),
                                                                            IQRVang = mean(IQRVang),
                                                                            cndt = cndt[1])



#use summarise to get trial measures and number of datapoints.
eyebytrial_auto <- eyeauto %>% group_by(pp_id, exp_id,trialn) %>% summarize(AvgHang = mean(hangle),
                                                                                StdHang = sd(hangle),
                                                                                MedHang = median(hangle),
                                                                                AvgVang = mean(vangle),
                                                                                StdVang = sd(vangle),
                                                                                MedVang = median(vangle),
                                                                                Datapoints = n(),
                                                                                IQRHang = IQR(hangle),
                                                                                IQRVang = IQR(vangle),
                                                                                nframes = n_distinct(world_frame_idx),
                                                                            Auto = Auto[1],
                                                                            Manual = Manual[1],
                                                                            cndt = condition[1])

#filter the eyedata by low nframes (at least 10s) and high std.
#eyebytrial_filtered <- eyebytrial %>% filter(nframes > 300, StdVang <30, StdHang<30)
#eyebytrial_unfiltered <- eyebytrial %>% filter(nframes < 300, StdVang >30, StdHang>30)
#eyebytrial_frames <- eyebytrial %>% filter(nframes < 300)
#eyebytrial_variable <- eyebytrial %>% filter(StdVang > 30, StdHang > 30)


###A VERY CONSERVATIVE FILTER OF EXTREME TRIALS###
eyebytrial_auto_filtered <- eyebytrial_auto %>% filter(nframes > 200, IQRHang < 20)


#now average by condition. Keep the number of trials going into each estimate.
eyebycondition_auto <- eyebytrial_auto_filtered %>% group_by(pp_id, exp_id, Auto, Manual) %>% summarize(AvgHang = mean(AvgHang),
                                                                                             StdHang = mean(StdHang),
                                                                                             MedHang = mean(MedHang),
                                                                                             AvgVang = mean(AvgVang),
                                                                                             StdVang = mean(StdVang),
                                                                                             MedVang = mean(MedVang),
                                                                                             Trialn = n(),
                                                                                             IQRHang = mean(IQRHang),
                                                                                             IQRVang = mean(IQRVang),
                                                                                             cndt = cndt[1])






##plot heatmaps for auto.
###PLOT HEATMAPPED GAZE POINTS POOLED ACROSS CONDITIONS
ggplot(eyeauto, mapping = aes(x=hangle, y=vangle)) +
  facet_grid(exp_id~.) +
  stat_density_2d(aes(fill=condition, alpha = ..level..), geom="polygon")+
  scale_fill_manual(values=c("Free_Fix"="#e6194b","Fix_Fix"="#3cb44b","Free_Free"="#0082c8","Fix_Free"="#ffe119")) +
  ylab("Angle from Horizon") + xlab("Angle From Centre") +
  coord_cartesian(ylim = c(-20,5), xlim = c(-10, 25)) + theme_light() + ggtitle("Automation Period")

ggplot(eyemanual, mapping = aes(x=hangle, y=vangle)) +
  facet_grid(exp_id~.) +
  stat_density_2d(aes(fill=condition, alpha = ..level..), geom="polygon")+
  scale_fill_manual(values=c("Free_Fix"="#e6194b","Fix_Fix"="#3cb44b","Free_Free"="#0082c8","Fix_Free"="#ffe119")) +
  ylab("Angle from Horizon") + xlab("Angle From Centre") +
  coord_cartesian(ylim = c(-20,5), xlim = c(-10, 25)) + theme_light() + ggtitle("Manual Period")





##save files as spreadsheet.
write_csv(eyebycondition_auto, path = "Data/EyetrackingData/AutomationPeriod_EyeTrackingbyCondition_CrowSparrow.csv")
write_csv(eyebycondition_manual, path = "Data/EyetrackingData/ManualPeriod_EyeTrackingbyCondition_CrowSparrow.csv")

#combine for plotting
eyebycondition_auto["Period"] = "Auto"
eyebycondition_manual["Period"] = "Manual"

eyebyC <- rbind(eyebycondition_auto, eyebycondition_manual)


#plot summary data.
ggplot(eyebycondition_auto, aes(x=cndt, y=StdHang, colour=exp_id, group=exp_id)) +# facet_wrap(~exp_id) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (sd(x)/sqrt(length(x))),
               fun.ymax = function(x) mean(x) + (sd(x)/sqrt(length(x))),
               geom = "pointrange",position = position_dodge(width=.5)) +
  ylab("StdHangle (degrees)") + theme_classic() + expand_limits(y=c(5,10)) + ggtitle("Automation Period")

ggplot(eyebycondition_manual, aes(x=cndt, y=StdHang, colour=exp_id, group=exp_id)) +# facet_wrap(~exp_id) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (sd(x)/sqrt(length(x))),
               fun.ymax = function(x) mean(x) + (sd(x)/sqrt(length(x))),
               geom = "pointrange",position = position_dodge(width=.5)) +
  ylab("StdHangle (degrees)") + theme_classic() + expand_limits(y=c(5,10)) + ggtitle("Manual Period")

ggplot(eyebyC, aes(x=cndt, y=IQRHang, colour=exp_id, group=exp_id)) + facet_wrap(~Period) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (sd(x)/sqrt(length(x))),
               fun.ymax = function(x) mean(x) + (sd(x)/sqrt(length(x))),
               geom = "pointrange",position = position_dodge(width=.5)) +
  ylab("IQRHang (degrees)") + theme_classic() + expand_limits(y=c(0,10)) + ggtitle("Driving Period")



#####STEERING DATA######

steer_sparrowdata = read_csv("Data/Sparrow17_010318.csv")
steer_NoLoad1data = read_csv("Data/NoLoad1_010318.csv")
steer_NoLoad2data = read_csv("Data/NoLoad2_010318.csv")
steer_crowdata = read_csv("Data/Crow17_010318.csv")

###These steering data are averages per trial.

#STEERING HEADERS: 'pp', 'tn_bycn','trialn', 'cn', 'exp_id', 'Fix', 'Bend', 'SB','RMS',...
#'AvgYaw','StdYaw','AvgSWA','StdSWA','AvgSWVel'

steerdata = rbind(steer_sparrowdata,steer_NoLoad1data,steer_NoLoad2data,steer_crowdata)

#remove pp1-4 as sparrow saving "-1".
steerdata <- steerdata %>% filter(pp %in% c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

steerdata <- steerdata %>% mutate(Auto = ifelse(Fix=='Free_Fix' | Fix=='Free_Free','Free',
                                            ifelse(Fix=='Fix_Fix' | Fix=='Fix_Free', 'Fix', NA)))
steerdata <- steerdata %>% mutate(Manual = ifelse(Fix=='Fix_Free' | Fix=='Free_Free','Free',
                                              ifelse(Fix=='Free_Fix' | Fix=='Fix_Fix', 'Fix', NA)))

#filter extreme.
steerdata <- steerdata %>% filter(abs(SB) <5)

steerbycondition <- steerdata %>% group_by(pp, exp_id, Auto, Manual) %>% summarize(AvgYaw = mean(AvgYaw),
                                                                                   StdYaw = mean(StdYaw),
                                                                                   AvgSWA = mean(AvgSWA),
                                                                                   StdSWA = mean(StdSWA),
                                                                                   AvgSWVel = mean(AvgSWVel),
                                                                                   SB = mean(SB),
                                                                                   RMS = mean(RMS),
                                                                                   cndt = Fix[1])
write_csv(steerbycondition, path = "Data/EyetrackingData/SteeringByCondition_CrowSparrow.csv")


#plot data.
ggplot(steerbycondition, aes(x=cndt, y=RMS, colour=exp_id)) +# facet_wrap(~exp_id) +
  stat_summary(fun.data=mean_se,position = position_dodge(width=.5)) +
  #geom_point(position=position_jitterdodge(), size=1) +
  ylab("RMS (m)") + theme_classic() + expand_limits(y=1)