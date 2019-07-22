library("tidyverse")
library(openxlsx)
library(ggplot2)
library(wesanderson)


# Output measures Orca19

###STEERING
steerdata <- read_csv("C:/VENLAB data/SP_18-19/Data/Orca19_Steering_Only/collated_steering.csv")
saveRDS(steerdata, "orca19_steerdata.rds")

head(steerdata)

#RT = need point at which autoflag == True minus Onset Time.


#add left bend
steerdata <- steerdata %>% 
  mutate(bend = ifelse(bend < 0, "left", "right"))

#Filter random data
balanced_steerdata <- steerdata %>%
  filter(design == "balanced")


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
#binary whether they disengaged

#refactor sab
balanced_steerdata <- balanced_steerdata %>% 
  ungroup() %>% 
  mutate (failure_type = case_when( sab == -0.20073596 ~ ".2",
                                    sab == -0.3039716 ~ ".3",
                                    sab == -0.52191351 ~ ".5",
                                    sab == -1.19868047 ~ "1.19",
                                    sab == -5.72957795 ~ "5.72"
  )) %>% 
  mutate(failure_type = as.factor(failure_type),
         cogload = as.factor(cogload),
         ppid = as.factor(ppid))


#create factors


# create averages for trials 
balanced_steer_trialavgs <- balanced_steerdata  %>% 
  group_by(ppid, cogload, failure_type, trialn) %>% 
  summarize(RT = disengage_RT(onsettime, timestamp_trial, autoflag),
            SWA_var = sd(swa),
            mean_SWA_vel = mean(diff(swa)),
            SDLP = sd(steeringbias),
            SB = mean(steeringbias),
            RMS = mean(sqrt(steeringbias^2)),
            YR_var = sd(yawrate_seconds),
            disengaged = ifelse(is.na(RT), 0, 1),
            sab = first(sab)
  )

  



ggplot(balanced_steer_trialavgs, aes(x = sab, y = RT)) +
  geom_point(alpha = .5)

#create average data points
average_balanced_steer_trialavgs <- balanced_steer_trialavgs %>%
  group_by(sab) %>%
  summarise(meanRT = mean(RT),
            meanRMS = mean(RMS)
)


# Plot average data
#### Not sure why this does not work!!
ungroup_df <- ungroup(average_balanced_steer_trialavgs)

head(ungroup_df)
ggplot(ungroup_df, aes(x = as.numeric(sab), y = meanRT)) +
  geom_line()

# RT by failure typpe and participant
ggplot(balanced_steer_trialavgs, aes(x = as.numeric(sab), y = RT, col = ppid)) +
  geom_point()

ggplot(balanced_steer_trialavgs, aes(x = sab, y = RMS, col = ppid)) +
  geom_point(alpha = 0.6) +
  scale_fill_manual(values = wes_palette(n = 2, name="Royal1"))


# Steering error by failure type and participant
ggplot(balanced_steer_trialavgs, aes(x = sab, y = RMS, col = ppid)) +
  geom_col()



# Box plots


ggplot(balanced_steer_trialavgs, aes(x = cogload, y = RMS, color = cogload)) +
  geom_boxplot(outlier.colour="red", outlier.shape = 9,
                outlier.size = 8) +
  stat_summary(fun.y = mean, geom="point", shape = 10, size = 4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = .5)



ggplot(balanced_steer_trialavgs, aes(x = cogload, y = RMS, color = cogload, fill = sab)) +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette(n = 4, name="GrandBudapest2"))


ggplot(balanced_steer_trialavgs, aes(x = cogload, y = mean_SWa_vel, color = cogload, fill = sab)) +
  geom_boxplot() +
  scale_fill_manual(values = wes_palette(n = 4, name="GrandBudapest2"))
 
