library("tidyverse")


#for students for ORCA

###GAZE
steergazedata <- read_csv("C:/Users/psccmo/Orca18_Analysis/GazeAndSteering_longFormat_080419.csv")

#add left bend
steergazedata <- steergazedata %>% 
  mutate(bend = ifelse(trialtype_signed < 0, "left", "right"))

steergazedata[steergazedata$bend == "left", "hangle"] <- steergazedata[steergazedata$bend == "left", "hangle"]*-1 



steergazedata <- steergazedata %>% 
  rename(SWV = SWA) %>% 
  mutate(SWA = SWV * 90)

#create function retrieving RT.
disengage_RT <- function(onsettime, timestamp_trial, autoflag){
  
  #pick first frame where autoflag == false, then take the timestamp and minus the onset_time
  auto_false <- which(autoflag == "FALSE")
  disengage_index <- first(auto_false)
  disengage_trialtime <- timestamp_trial[disengage_index]
  onset_time <- first(onsettime)
  RT <- disengage_trialtime - onset_time #can be negative
  return(RT)
}


steergaze_trialavgs <- steergazedata  %>% 
  group_by(ppid, radius, yawrate_offset, cogload, block, count) %>% 
  summarize(RT = disengage_RT(OnsetTime, timestamp_trial, AutoFlag),
            SWA_var = sd(SWA),
            mean_SWA_vel = mean(diff(SWA)),
            SDLP = sd(SteeringBias),
            SB = mean(SteeringBias),
            RMS = mean(sqrt(SteeringBias^2)),
            hang_IQR = IQR(hangle),
            vang_IQR = IQR(vangle),
            YR_var = sd(YawRate_seconds),
            disengaged = ifelse(is.na(RT), 0, 1)
  )

steergaze_trialavgs <- steergaze_trialavgs %>% 
  group_by(ppid, radius, yawrate_offset, cogload) %>% 
  mutate(trialn_cndt = 1:n())

steergaze_trialavgs <- steergaze_trialavgs %>% 
  mutate(failure_type = yawrate_offset)

steergaze_trialavgs$failure_type <- as.factor(steergaze_trialavgs$failure_type)


levels(steergaze_trialavgs$failure_type)[1] = "Sudden"
levels(steergaze_trialavgs$failure_type)[2] = "Gradual"
levels(steergaze_trialavgs$failure_type)[3] = "None"
levels(steergaze_trialavgs$failure_type)[4] = "None"

steergaze_trialavgs$cogload <- factor(steergaze_trialavgs$cogload, levels = c("None", "Easy", "Hard")) #reorder

#plot some key variables
#RT
ggplot(filter(steergaze_trialavgs, failure_type == "Sudden", RT > 0), aes(x = factor(cogload), col = factor(cogload), y = RT)) + 
  geom_boxplot() +
  facet_grid(.~radius)



