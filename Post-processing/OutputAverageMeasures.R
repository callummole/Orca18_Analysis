library("tidyverse")

#FOr DA19 COgnitve load report. So that RW can do his own analysis of participant means.

###GAZE
steergazedata = read_csv("F:/TempStorageForMDrive/SparrowCrow/GazeAndSteeringData_SimErrorRemoved.csv")

#steergazedata$exp_id <- factor(steergazedata$exp_id, levels = c('NoLoad1','NoLoad2','Sparrow17','Crow17')) #reorder
#revalue(steergazedata$exp_id, c("NoLoad1"="NoLoad1", "NoLoad2"="NoLoad2", "Sparrow17" = "Count", "Crow17" = "Nback"))
#steergazedata$exp_id <- factor(steergazedata$exp_id, levels = c('NoLoad1','NoLoad2','Count','Nback')) #rename

#newgazedata <- read_csv("M:/SparrowCrow/EyetrackingData/GazeAndSteering_longFormat_160518.csv") #new gaze bias calculated
#newgazedata$trialid <- paste(newgazedata$exp_id, newgazedata$pp_id, newgazedata$trialn, sep="_")

#mutate condition codes into auto and manual
steergazedata <- steergazedata %>% mutate(Auto = case_when(condition=='Free_Fix' | condition=='Free_Free' ~ 'A_Free',
                                                           condition=='Fix_Fix' | condition=='Fix_Free' ~ 'A_Fix'))

steergazedata <- steergazedata %>% mutate(Manual = case_when(condition=='Fix_Free' | condition=='Free_Free' ~ 'M_Free',
                                                             condition=='Free_Fix' | condition=='Fix_Fix' ~ 'M_Fix'))

steergazedata[steergazedata$Bend == "Left", "hangle"] <- steergazedata[steergazedata$Bend == "Left", "hangle"]*-1 

steergazedata_classes <- steergazedata %>% filter(sample_class ==1 | sample_class == 4 ) #only include fixation & smooth pursuit segments.
steergazedata_classfilter <- steergazedata_classes %>% filter(vangle > -24 & vangle<0 & hangle >-39.5 & hangle<39.5 & lookahead < 60 & abs(gazebias) < 10)
steergazedata_classfilter <- steergazedata_classfilter %>% filter(!(pp_id %in% c(6,11,21))) #these three participant seem to be poorly calibrated.

SG_avgs <- steergazedata_classfilter %>% filter(timestamp > 1 & timestamp < 18) %>% group_by(trialid, autoflag) %>% summarize(gazebias = mean(gazebias), IQR_v = IQR(vangle),
                                                                                                                              IQR_h = IQR(hangle),
                                                                                                                              v_med = median(vangle),
                                                                                                                              h_med = median(hangle),
                                                                                                                              v_var = sd(vangle),
                                                                                                                              h_var = sd(hangle),
                                                                                                                              lookahead = mean(lookahead),
                                                                                                                              pp_id = pp_id[1],
                                                                                                                              Auto = Auto[1],
                                                                                                                              exp_id = exp_id[1],
                                                                                                                              Manual = Manual[1],
                                                                                                                              Bend = Bend[1],
                                                                                                                              trialn = trialn[1],
                                                                                                                              obs = n())

#set a limit of needing at least 3s of gaze data from automation/manual period = 180 obs.
SG_avgs_filter <- SG_avgs %>% filter(obs>120)

SG_MANAVGS <- filter(SG_avgs_filter, autoflag==0)

SG_AUTOVGS <- filter(SG_avgs_filter, autoflag==1)

pp10 <- filter(SG_avgs, pp_id == 10, exp_id == "NoLoad2", Auto == "A_Free")

SG_MAN_cndt <- SG_MANAVGS  %>% group_by(exp_id, Manual, Auto, pp_id) %>% summarize(gazebias = mean(gazebias), IQR_v = mean(IQR_v),
                                                                                                 IQR_h= mean(IQR_h),
                                                                                                 lookahead = mean(lookahead),
                                                                                                 h_med = mean(h_med),
                                                                                                 v_med = mean(v_med),
                                                                                                 h_var = mean(h_var),
                                                                                                 v_var = mean(v_var))

SG_AUTO_cndt <- SG_AUTOVGS  %>% group_by(exp_id, Auto, pp_id) %>% summarize(gazebias = mean(gazebias), IQR_v = mean(IQR_v),
                                                                                             IQR_h= mean(IQR_h),
                                                                                             lookahead = mean(lookahead),
                                                                                             h_med = mean(h_med),
                                                                                             v_med = mean(v_med),
                                                                                             h_var = mean(h_var),
                                                                                             v_var = mean(v_var))

SG_MAN_cndt$exp_id[SG_MAN_cndt$exp_id=="Crow17"] <- "Nback"
SG_MAN_cndt$exp_id[SG_MAN_cndt$exp_id=="Sparrow17"] <- "Count"
SG_MAN_cndt <- SG_MAN_cndt %>% unite(condition, exp_id, Auto, Manual) 


SG_AUTO_cndt$exp_id[SG_AUTO_cndt$exp_id=="Crow17"] <- "Nback"
SG_AUTO_cndt$exp_id[SG_AUTO_cndt$exp_id=="Sparrow17"] <- "Count"
SG_AUTO_cndt <- SG_AUTO_cndt %>% unite(condition, exp_id, Auto) 



#AUTO MEASURES

#just change the name to save out a different variable
SG_cndtavgs_AUTO_gazebias <- SG_AUTO_cndt %>% select(condition, pp_id, gazebias) %>%  spread(key = condition, value = gazebias)
write.csv(SG_cndtavgs_AUTO_gazebias, file = "gaze_averages_AUTO_gazebias.csv")

SG_cndtavgs_AUTO_h_var <- SG_AUTO_cndt %>% select(condition, pp_id, h_var) %>%  spread(key = condition, value = h_var)
write.csv(SG_cndtavgs_AUTO_h_var, file = "gaze_averages_AUTO_h_var.csv")

#just change the name to save out a different variable
SG_cndtavgs_MANUAL_gazebias <- SG_MAN_cndt %>% select(condition, pp_id, gazebias) %>%  spread(key = condition, value = gazebias)
write.csv(SG_cndtavgs_MANUAL_gazebias, file = "gaze_averages_MANUAL_gazebias.csv")

SG_cndtavgs_MANUAL_h_var <- SG_MAN_cndt %>% select(condition, pp_id, h_var) %>%  spread(key = condition, value = h_var)
write.csv(SG_cndtavgs_MANUAL_h_var, file = "gaze_averages_MANUAL_h_var.csv")






####Steering

steerdata = read_csv("F:/TempStorageForMDrive/SparrowCrow/FramebyFrameSteeringData_SimErrorRemoved.csv")

steeravgs <- steerdata %>% filter(autoflag==0, timestamp<18) %>% group_by(trialid) %>% summarize(OSBmn = mean(OSB), SWAmn = mean(SWA),
                                                                                                 SWAvar = sd(SWA),
                                                                                                 OSBvar = sd(OSB),
                                                                                                 RMS = sqrt(mean(OSB^2)),
                                                                                                 pp_id = pp_id[1],
                                                                                                 Auto = Auto[1],
                                                                                                 exp_id = exp_id[1],
                                                                                                 Manual = Manual[1],
                                                                                                 Bend = Bend[1],
                                                                                                 trialn = trialn[1])

#average over autotime
steercndts <- steeravgs %>% group_by(exp_id, Auto, Manual, pp_id) %>% summarize(OSBmn = mean(OSBmn), SWAmn = mean(SWAmn),
                                                                                SWAvar = mean(SWAvar),
                                                                                OSBvar = mean(OSBvar),
                                                                                RMS = mean(RMS))

###DA19### ####THIS SHOULD BE IT's OWN SCRIPT.

##to wide data and save for frequentist statistics.
head(steercndts)

#Change names.
steercndts_names <- steercndts
steercndts_names$exp_id[steercndts_names$exp_id=="Crow17"] <- "Nback"
steercndts_names$exp_id[steercndts_names$exp_id=="Sparrow17"] <- "Count"

steercndts_names <- steercndts_names %>% unite(condition, exp_id, Auto, Manual) 

#just change the name to save out a different variable
steercndts_wide_OSB <- steercndts_names %>% select(condition, pp_id, OSBmn) %>%  spread(key = condition, value = OSBmn)
write.csv(steercndts_wide_OSB, file = "steering_averages_OSB.csv")

steercndts_wide_RMS <- steercndts_names %>% select(condition, pp_id, RMS) %>%  spread(key = condition, value = RMS)
write.csv(steercndts_wide_RMS, file = "steering_averages_RMS.csv")

steercndts_wide_SWAmn <- steercndts_names %>% select(condition, pp_id, SWAmn) %>%  spread(key = condition, value = SWAmn)
write.csv(steercndts_wide_SWAmn, file = "steering_averages_SWAmn.csv")

steercndts_wide_SWAvar <- steercndts_names %>% select(condition, pp_id, SWAvar) %>%  spread(key = condition, value = SWAvar)
write.csv(steercndts_wide_SWAvar, file = "steering_averages_SWAvar.csv")
