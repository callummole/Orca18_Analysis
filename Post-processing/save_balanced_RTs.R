setwd("C:/git_repos/Orca18_Analysis/Post-Processing/")

#load steergaze data


#steergazedata <- read_csv("../Data/Orca19_collated_steering.csv")  
#saveRDS(steergazedata, "../Data/Orca19_collated_steering.rds")
steergazedata <- readRDS("../Data/Orca19_collated_steering.rds")  


steergazedata <- steergazedata %>% 
  rename(swv = swa) %>% 
  mutate(swa = swv * 90)

#mirror data
steergazedata <- steergazedata %>% 
  mutate(world_x_mirrored = if_else(bend == -1, world_x * -1, world_x),
         swa_mirrored = if_else(bend == -1, swa * -1, swa),
         swv_mirrored = if_else(bend == -1, swv * -1, swv),
         sb_mirrored = if_else(bend == -1, steeringbias * -1, steeringbias))


#unsophisticated calculation of TTC using rate of change in steering bias
calc_TTC <- function(sb_mirrored, sb_change){
  
  #need to incorporate the proximity to road edges and direction of travel.
  
  sb_mirrored <- array(sb_mirrored) #needs to be an array to use apply
  
  road_edges <- c(-1.5, 1.5)
  
  distance_to_edges <- apply(sb_mirrored, 1, function(x) road_edges - x)
  
  TTLCs <- apply(distance_to_edges,1, function(x) x / (sb_change * 60))
  
  TTLC <- apply(TTLCs, 1, max)
  
  return(TTLC)
  
}

#add ttlc
steergazedata <- steergazedata %>% 
  mutate(sb_change = prepend(diff(sb_mirrored), NA),
         TTLC = calc_TTC(sb_mirrored, sb_change)
  )


#add RT and disengage flag.
disengage_RT <- function(onsettime, timestamp_trial, autoflag){
  
  #pick first frame where autoflag == false, then take the timestamp and minus the onset_time
  auto_false <- which(autoflag == "FALSE")
  disengage_index <- first(auto_false)
  disengage_trialtime <- timestamp_trial[disengage_index]
  onset_time <- first(onsettime)
  RT <- disengage_trialtime - onset_time #can be negative
  return(RT)
  
}

#calculate RT
steergazedata <- steergazedata  %>% 
  group_by(ppid, sab, cogload, trialn) %>% 
  mutate(RT = disengage_RT(onsettime, timestamp_trial, autoflag),
         disengaged = ifelse(is.na(RT), 0, 1) #whether or not they actually took over.
  )

#create unique trial id
steergazedata <- steergazedata %>% 
  mutate(trialid = paste(ppid, cogload, trialn, sep = "_"))

steergazedata$cogload <- as.factor(steergazedata$cogload)


#rename cogload factors so that it ameks sense
steergazedata$cogload<- plyr::mapvalues(steergazedata$cogload, from = c("None", "Middle"), to = c("noload", "load"))

balanced_steerdata <- steergazedata %>% 
  filter(design=="balanced")

random_steerdata <- steergazedata %>% 
  filter(design=="random")

#refactor sab for 
balanced_steerdata <- balanced_steerdata %>% 
  ungroup() %>% 
  mutate (failure_type = case_when( sab == -0.3039716 ~ ".3",
                                    sab == -0.52191351 ~ ".5",
                                    sab == -1.19868047 ~ "1.19",
                                    sab == -5.72957795 ~ "5.72"
                                    )
  ) %>% 
  mutate(failure_type = as.factor(failure_type),
         cogload = as.factor(cogload),
         ppid = as.factor(ppid))



#we now have two data frames. "balanced_steerdata" and "random_steerdata"
steergaze_trialavgs <- steergazedata  %>% 
  ungroup() %>% 
  filter(design == "balanced") %>% 
  group_by(ppid, cogload, trialn) %>% 
  summarize(RT = first(RT),
            disengaged = first(disengaged), #whether or not they actually took over.
            premature = ifelse(RT <= 0, 1, 0),
            sab = first(sab),
            onsettime = first(onsettime),
            design = first(design),
            simTTLC = first(simulated_ttlc))


write.csv(steergaze_trialavgs, "balanced_RTs.csv")

#we now have two data frames. "balanced_steerdata" and "random_steerdata"
steergazedata  %>% 
  ungroup() %>% 
  filter(design == "random") %>% 
  group_by(ppid, cogload, trialn) %>% 
  summarize(RT = first(RT),
            disengaged = first(disengaged), #whether or not they actually took over.
            premature = ifelse(RT <= 0, 1, 0),
            sab = first(sab),
            onsettime = first(onsettime),
            design = first(design),
            simTTLC = first(simulated_ttlc)) %>% 
  write.csv(., "random_RTs.csv")


