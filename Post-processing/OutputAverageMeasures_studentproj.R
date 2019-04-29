library("tidyverse")
library(openxlsx)


#for students for ORCA

###GAZE
steergazedata <- read_csv("C:/Users/psccmo/Orca18_Analysis/GazeAndSteering_longFormat_080419.csv")
saveRDS(steergazedata, "orca_steergazedata.rds")

head(steergazedata)

#RT = need point at which autoflag == True minus Onset Time.


#add left bend
steergazedata <- steergazedata %>% 
  mutate(bend = ifelse(trialtype_signed < 0, "left", "right"))

steergazedata[steergazedata$bend == "left", "hangle"] <- steergazedata[steergazedata$bend == "left", "hangle"]*-1 



steergazedata <- steergazedata %>% 
  rename(SWV = SWA) %>% 
  mutate(SWA = SWV * 90)

#########
#yawrate offset wasn't taken into account when saving files, so data was over-written.
#so we essentially only have the last six trials per radii and cogload conditions.
#trialn is also the 1:n of the experiment. NOT the within condition exp.
#########

steergazedata <- steergazedata %>% 
  mutate(trialn_corrected = case_when(block == 0 ~ trialn,
                                      block == 1 ~ trialn,
                                      block == 2 ~ trialn+6))


#trial <- filter(steergazedata, trialcode == "Orca18_None_2_p30_40_3")


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


#across yawrates, radii, cogload.

#TRIAL_MEASURES:
#steering wheel variability
#steering wheel velocity
#RMS
#yaw rate variability.
#SDLP
#RT
#hangle IQR
#vangle IQR
#binary whether they disengaged

#create factors
steergazedata$ppid <- as.factor(steergazedata$ppid)
steergazedata$radius <- as.factor(steergazedata$radius)
steergazedata$yawrate_offset <- as.factor(steergazedata$yawrate_offset)
steergazedata$cogload <- as.factor(steergazedata$cogload)


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
 
ggplot(steergaze_trialavgs, aes(x = yawrate_offset, y= SWA_var, col = cogload)) +
  geom_jitter(alpha = .2)

ggplot(steergaze_trialavgs, aes(x = cogload, y= RT, col = radius)) +
  geom_jitter(alpha = .2)

#calculate grand means of the main measures for imputing.
grandmeans <- steergaze_trialavgs %>% 
  ungroup() %>% 
  filter(RT > 0) %>% 
  group_by(radius, yawrate_offset, cogload) %>%
  summarise(mn_RT = mean(RT, na.rm= T),
            sd_RT = sd(RT, na.rm = T),
            mn_SWA_var = mean(SWA_var),
            sd_SWA_var = sd(SWA_var),
            mn_SWA_vel = mean(mean_SWA_vel),
            sd_SWA_vel = sd(mean_SWA_vel),
            mn_SDLP = mean(SDLP),
            sd_SDLP = sd(SDLP),
            mn_SB = mean(SB),
            sd_SB = sd(SB),
            mn_RMS = mean(RMS),
            sd_RMS = sd(RMS),
            mn_hang_IQR = mean(hang_IQR),
            sd_hang_IQR = sd(hang_IQR),
            mn_vang_IQR = mean(vang_IQR),
            sd_vang_IQR = sd(vang_IQR),
            mn_disengaged = mean(disengaged),
            sd_disengaged = sd(disengaged))

#expand so there are entries for combination of condition
steergaze_trialavgs_expanded <- steergaze_trialavgs %>% 
  complete(ppid, radius, yawrate_offset, cogload, 
           fill = list(RT = 99,
                       SWA_var = 99,
                       mean_SWA_vel = 99,
                       SDLP = 99,
                       SB = 99,
                       RMS = 99,
                       hang_IQR = 99,
                       vang_IQR = 99,
                       YR_var = 99,
                       disengaged = 99)
           )


#add a monotonically increasing trial number.
steergaze_trialavgs_expanded <- steergaze_trialavgs_expanded %>% 
  group_by(ppid, radius, yawrate_offset, cogload) %>% 
  mutate(trialn_cndt = 1:n())
  

#join grandmeans and steergaze_trialavgs_expanded.
joined <- left_join(steergaze_trialavgs_expanded, grandmeans, by = c("radius","yawrate_offset","cogload"))


#impute means from a normal distribution if the value is 99
joined_imputed <- joined %>% 
  ungroup() %>% 
  group_by(ppid, radius, yawrate_offset, cogload, trialn_cndt) %>% 
  mutate(RT = ifelse(RT == 99, rnorm(1, mn_RT, sd_RT/4), RT),
         SWA_var = ifelse(SWA_var == 99, rnorm(1, mn_SWA_var, sd_SWA_var/2), SWA_var),
         mean_SWA_vel = ifelse(mean_SWA_vel == 99, rnorm(1, mn_SWA_vel, sd_SWA_vel/2), mean_SWA_vel),
         SDLP = ifelse(SDLP == 99, rnorm(1, mn_SDLP, sd_SDLP/2), SDLP),
         SB = ifelse(SB == 99, rnorm(1, mn_SB, sd_SB/2), SB),
         RMS = ifelse(RMS == 99, rnorm(1, mn_RMS, sd_RMS/2), RMS),
         hang_IQR = ifelse(hang_IQR == 99, rnorm(1, mn_hang_IQR, sd_hang_IQR/2), hang_IQR),
         vang_IQR = ifelse(vang_IQR == 99, rnorm(1, mn_vang_IQR, sd_vang_IQR/2), vang_IQR),
         disengaged = ifelse(disengaged == 99, rnorm(1, mn_disengaged, sd_disengaged/2), disengaged)
  )




#refactor yawrate_offset
steergaze_trialavgs <- steergaze_trialavgs %>% 
  mutate(failure_type = case_when(yawrate_offset %in% c(-.2, .15) ~ "None",
                                   yawrate_offset == -9 ~ "Sudden",
                                   yawrate_offset == -1.5 ~ "Gradual"
                                   ))

steergazedata <- steergazedata %>% 
  mutate(failure_type = case_when(yawrate_offset %in% c(-.2, .15) ~ "None",
                                  yawrate_offset == -9 ~ "Sudden",
                                  yawrate_offset == -1.5 ~ "Gradual"
  ))


joined_imputed <- joined_imputed %>% 
  mutate(failure_type = yawrate_offset)

joined_imputed$failure_type <- as.double(joined_imputed$failure_type)

levels(joined_imputed$failure_type)[1] = "Sudden"
levels(joined_imputed$failure_type)[2] = "Gradual"
levels(joined_imputed$failure_type)[3] = "None"
levels(joined_imputed$failure_type)[4] = "None"

#condition averages
#mean RT
#median RT
#mean SWA var
#mean SWA vel
steergaze_condtavgs <- joined_imputed  %>% 
  group_by(ppid, radius, failure_type, cogload) %>% 
  filter(RT > 0) %>% 
  summarize(mn_RT = mean(RT, na.rm = T),
            med_RT = median(RT, na.rm = T),
            mn_SWA_var = mean(SWA_var),
            mn_SWA_vel = mean(mean_SWA_vel),
            mn_SDLP = mean(SDLP),
            mn_SB = mean(SB),
            mn_RMS = mean(RMS),
            mn_hang_IQR = mean(hang_IQR),
            mn_vang_IQR = mean(vang_IQR),
            perc_takeover = sum(disengaged)/n()
            )

head(steergaze_condtavgs)

write_csv(steergaze_condtavgs, "orca_conditionaverages_longformat.csv")

#wide format. 
steergaze_united <- steergaze_condtavgs %>% 
  unite("condition", 
        failure_type, cogload, radius,
        sep = '_')

#need to save separate measures
OUT <- openxlsx::createWorkbook()
addsheet <- function(OUT, varname){
 
  var = as.symbol(varname)
  steercndts_wide <- steergaze_united %>% 
    select(condition, ppid, !!var) %>%  
    spread(key = condition, value = !!var)

  addWorksheet(OUT, varname)
  writeData(OUT, sheet=varname, x = steercndts_wide)  
  
}

addsheet(OUT, "mn_RT")
addsheet(OUT, "med_RT")
addsheet(OUT, "mn_SWA_var")
addsheet(OUT, "mn_SWA_vel")
addsheet(OUT, "mn_SDLP")
addsheet(OUT, "mn_SB")
addsheet(OUT, "mn_RMS")
addsheet(OUT, "mn_hang_IQR")
addsheet(OUT, "mn_vang_IQR")
addsheet(OUT, "perc_takeover")

openxlsx::saveWorkbook(OUT, "orca_conditionaverages_wideformat2.xlsx")
