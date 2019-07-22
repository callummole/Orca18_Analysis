library("tidyverse")
library("brms")
library(ggplot2)
#library(brmstools)
library(cowplot)
#library(dplyr)
#library(tidyr)


#####STEERING DATA######

steer_Orca19_data = read_csv("C:/VENLAB data/SP_18-19/DataOrca19_Steering_Only/collated_steering.csv")
#steer_NoLoad1data = read_csv("Data/NoLoad1_010318.csv")
#steer_NoLoad2data = read_csv("Data/NoLoad2_010318.csv")
#steer_crowdata = read_csv("Data/Crow17_010318.csv")

###These steering data are averages per trial.

#STEERING HEADERS: 'pp', 'tn_bycn','trialn', 'cn', 'exp_id', 'Fix', 'Bend', 'SB','RMS',...
#'AvgYaw','StdYaw','AvgSWA','StdSWA','AvgSWVel'

#steerdata = rbind(steer_sparrowdata,steer_NoLoad1data,steer_NoLoad2data,steer_crowdata)

#remove pp1-4 as sparrow saving "-1".
#steerdata <- steerdata %>% filter(pp %in% c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

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