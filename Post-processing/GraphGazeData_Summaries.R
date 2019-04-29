library("tidyverse")
library(ggExtra)
library(ggridges)
library(magrittr) #for extra pipe functions
library(zoo) #for rollmean moving average function
library(cowplot)
library(osfr) #Package for interfacing with OSF. https://rdrr.io/github/CenterForOpenScience/osfr/

#For Data Dictionary see: https://github.com/callummole/Trout18_Analysis/tree/master/Data

#theme for plots on TRANSITION grant.
theme_transition <- theme_classic() +
  theme(strip.background = element_rect(fill=NA,color=NA), 
        strip.text = element_text(face="bold",colour="black",size="8"), 
        axis.title = element_text(face="bold",colour="black",size="8"),
        axis.text.x = element_text(vjust=-.5),
        axis.text.y = element_text(vjust=.5),
        axis.text = element_text(face="plain",colour="black",size="7"),
        legend.text = element_text(face="plain",colour="black",size="7"),
        legend.title = element_text(face="bold",colour="black",size="8"),
        legend.key = element_blank(),
        panel.grid.major.y = element_line(color="grey85",size=.2, linetype = 2))

mypalette = c("#FF9966","#33CC33","#FF0000","#0000FF")


# Load data file. Is a large file, >1GB.
#only load the file if it isn't already loaded.
if (!exists("steergazedata")){
  steergazedata <- read_csv("C:/Users/psccmo/Trout18_Analysis/Data/GazeAndSteering_longFormat_25PPs_181218.csv")  
  
  #Append a column (currtimezero) which gives the time starting at zero
  steergazedata <- steergazedata %>% group_by(trialcode) %>% mutate(currtimezero= currtime - currtime[1] )
  
  ### Remove backup trials (obstacleoffset = .5 or condition == 99), or surplus manual trials (count > 100).
  steergazedata <-steergazedata  %>% filter(condition != 99,count < 100)
  
  #Make factors.
  steergazedata$sectiontype <- as.factor(steergazedata_onsrf$sectiontype)
  steergazedata$condition <- as.factor(steergazedata_onsrf$condition)
  
  
  #Trial code is a unique identifier
  #Block_PP_sectiontype_condition_count
}


#For the first step, we want to assess at the highest level variability across different conditions and blocks.
# The measures we are interested in at a first pass are the variability and the number of observations.

#For Diagnostics, let's look at active control for bends only. 
Bends_data <- steergazedata %>% filter( (posz > 50 | posz < -50), sectiontype %in% c(0))


#### TRIAL SUMMARIES #####
TrialData = GazeSummaries(Bends_data, 'trial')
TrialData_onlyFullTrials <- filter(TrialData, sectiontype %in% c(0,1,2))

### Let's assess whether we have a consistent number of datapoints ###
DataPoints_ppMean <- TrialData_onlyFullTrials %>% group_by(ID) %>% summarise(number_of_obs = mean(obs, na.rm=TRUE))

ggplot(TrialData_onlyFullTrials, aes(y =obs, x=ID)) + geom_point(alpha = .2) + ylab("Number of Obs") +facet_grid(block~.) +
  geom_point(data=DataPoints_ppMean, aes(y=number_of_obs, x=ID), colour = "Red", shape=17, size=2 )


### ARE THERE ANY Participants that are way different to the rest? #######
#some key metrics. 1) Lookahead
LookaheadMean_ppMean <- TrialData_onlyFullTrials %>% group_by(ID) %>% summarise(lookahead_mean = mean(lookahead_mean, na.rm=TRUE))
ggplot(TrialData_onlyFullTrials, aes(y =lookahead_mean, x=ID)) + geom_point(alpha = .2) + ylab("Lookahead_mean") +
  geom_point(data=LookaheadMean_ppMean, aes(y=lookahead_mean, x=ID), colour = "Red", shape=17, size=2 )

LookaheadMed_ppMean <- TrialData_onlyFullTrials %>% group_by(ID) %>% summarise(lookahead_med = mean(lookahead_med, na.rm=TRUE))
ggplot(TrialData_onlyFullTrials, aes(y =lookahead_med, x=ID)) + geom_point(alpha = .2) + ylab("Lookahead_med") + facet_grid(block~.) +
  geom_point(data=LookaheadMed_ppMean, aes(y=lookahead_med, x=ID), colour = "Red", shape=17, size=2 ) + ggtitle("Active Control, On Bends")

# V_Var
v_var_ppMean <- TrialData_onlyFullTrials %>% group_by(ID) %>% summarise(v_var = mean(v_var, na.rm=TRUE))
ggplot(TrialData_onlyFullTrials, aes(y =v_var, x=ID)) + geom_point(alpha = .2) + ylab("v_var") +facet_grid(block~.) +
  geom_point(data=v_var_ppMean, aes(y=v_var, x=ID), colour = "Red", shape=17, size=2 ) + ggtitle("Active Control, On Bends")

# V_med
ManualTrials <- filter(TrialData_onlyFullTrials, sectiontype == 0)
v_med_ppMean <- ManualTrials %>% group_by(ID) %>% filter(sectiontype == 0) %>%  summarise(v_med = mean(v_med, na.rm=TRUE))
ggplot(ManualTrials, aes(y =v_med, x=ID)) + geom_point(alpha = .2) + ylab("v_med") +facet_grid(block~.) +
  geom_point(data=v_med_ppMean, aes(y=v_med, x=ID), colour = "Red", shape=17, size=2 )+ ggtitle("Active Control, On Bends")

#density plot of v_med.
times = array(1:4)
distances = times * 8.0
vrad = atan(-1.2 / distances)
vangle = vrad * 180 / pi
ggplot(Bends_data, aes(vangle)) + geom_density() + coord_cartesian(xlim = c(-15, 10)) + geom_vline(xintercept = vangle, color = "red")



# h_med
ManualTrials <- filter(TrialData_onlyFullTrials, sectiontype == 0)
h_med_ppMean <- TrialData_onlyFullTrials %>% group_by(ID) %>% filter(sectiontype == 0) %>%  summarise(h_med = mean(h_med, na.rm=TRUE))
ggplot(TrialData_onlyFullTrials, aes(y =h_med, x=ID)) + geom_point(alpha = .2) + ylab("h_med") +facet_grid(block~.) +
  geom_point(data=h_med_ppMean, aes(y=h_med, x=ID), colour = "Red", shape=17, size=2 )+ ggtitle("Active Control, On Bends")


# gaze_bias
ManualTrials <- filter(TrialData_onlyFullTrials, sectiontype == 0)
gazebias_ppMean <- ManualTrials %>% group_by(ID) %>% filter(sectiontype == 0) %>%  summarise(gazebias = mean(gazebias, na.rm=TRUE))
ggplot(ManualTrials, aes(y =gazebias, x=ID)) + geom_point(alpha = .2) + ylab("gazebias") +facet_grid(block~.) +
  geom_point(data=gazebias_ppMean, aes(y=gazebias, x=ID), colour = "Red", shape=17, size=2 )+ ggtitle("Active Control, On Bends")

min(Bends_data$gazebias)
ggplot(filter(Bends_data, abs(gazedistance) < 10), aes(gazedistance)) + geom_density() + coord_cartesian(xlim = c(-10, 10)) 

###### Retrieve Summaries for Plotting ###########
PPCndtData = GazeSummaries(steergazedata, 'participant')




GroupCndtData = GazeSummaries(steergazedata, 'group')
#Filter out PID and Interp
GroupCndtData <- filter(GroupCndtData, sectiontype %in% c(0,1,2))

### Explore Trial Data for crazy values. ####

ggplot(data=filter(TrialData, sectiontype %in% c(0,1,2)), aes(y=v_med, x=trialcode, col = factor(sectiontype))) + geom_point() + 
  ylab("v_med") +
  theme_transition + theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
  
  




####### Plot Point Ranges for Whole Trial, for each Section ##########

####IQR Ranges #######
#Doesn't seem to be many trends across the entire trial

vplot <- PlotPointRange(GroupCndtData,'IQR_v','Vert IQR (degs)', c(-2,8))
hplot <- PlotPointRange(GroupCndtData, 'IQR_h', 'Horiz IQR (degs)', c(3,25))
IQRplot <- plot_grid(vplot,hplot, nrow=2)
add_cowplottitle(IQRplot, "Whole Trial")

vplot <- PlotPointRange(GroupCndtData, 'v_var','Vert Std (degs)')
hplot <- PlotPointRange(GroupCndtData, 'h_var','Horiz Std (degs)')
plot_grid(vplot,hplot, nrow=2)

lookaheadplot <- PlotPointRange(GroupCndtData, 'lookahead','Lookahead (m)')
lookaheadplot 


### SPLIT steergazedata_onsrf into sections

##BENDS
Bends_onsrf <- steergazedata_onsrf %>% filter( (posz > 50 | posz < -50), sectiontype %in% c(0,1,2))

Bends_Trials <- GazeSummaries(Bends_onsrf, 'trial')
ggplot(data=filter(Bends_Trials, sectiontype %in% c(0,1,2)), aes(y=v_med, x=trialcode, col = factor(sectiontype))) + geom_point() + 
  ylab("v_med") +
  theme_transition + theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 

Bends_Groups <- GazeSummaries(Bends_onsrf, 'group')
vplot <- PlotPointRange(Bends_Groups,'IQR_v','Vert IQR (degs)', c(-2,8))
hplot <- PlotPointRange(Bends_Groups, 'IQR_h', 'Horiz IQR (degs)', c(3,25))
IQRplot <- plot_grid(vplot,hplot, nrow=2)

#add title
add_cowplottitle(IQRplot, "Bends")

vplot <- PlotPointRange(Bends_Groups, 'v_var','Vert Std (degs)')
hplot <- PlotPointRange(Bends_Groups, 'h_var','Horiz Std (degs)')
varplot <- plot_grid(vplot,hplot, nrow=2)
add_cowplottitle(varplot, "Bends")


##STRAIGHTS
Straights_onsrf <- steergazedata_onsrf %>% filter( (posz > -50 & posz < 50), sectiontype %in% c(0,1,2))

Straights_Trials <- GazeSummaries(Straights_onsrf, 'trial')
ggplot(data=filter(Straights_Trials, sectiontype %in% c(0,1,2)), aes(y=IQR_v, x=trialcode, col = factor(sectiontype))) + geom_point() + 
  ylab("IQR_v") +
  theme_transition + theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 

Straights_Groups <- GazeSummaries(Straights_onsrf, 'group')
vplot <- PlotPointRange(Straights_Groups,'IQR_v','Vert IQR (degs)', c(-2,8))
hplot <- PlotPointRange(Straights_Groups, 'IQR_h', 'Horiz IQR (degs)', c(3,25))
IQRplot <- plot_grid(vplot,hplot, nrow=2)
#add title
add_cowplottitle(IQRplot, "Straights")

vplot <- PlotPointRange(Straights_Groups, 'v_var','Vert Std (degs)')
hplot <- PlotPointRange(Straights_Groups, 'h_var','Horiz Std (degs)')
varplot <- plot_grid(vplot,hplot, nrow=2)

#add title
add_cowplottitle(varplot, "Straights")


###### MY FUNCTIONS ########


# Create function that returns summary dataframe for a variable given as a string.
#https://stackoverflow.com/questions/49371260/using-variables-as-arguments-in-summarize

#function uses as.symbol to change the name of a variable
make_name <- function(var, suffix)(as.symbol(paste(var,suffix, sep=".")))


CalculateSummaries <- function(data, var, drop){
  # returns a dataframe with summaries for var
  # var should be a string that corresponds to a column in PPCndt_avgs, above
  # if drop = True, remove the group_by columns and only keep the columns added in summarize()
  
  #create symbols.
  var = as.symbol(var)
  #print(var)
  mymean = make_name(var,"mean")
  #print(mymean)
  myCIupper = make_name(var,"CIupper")
  myCIlower = make_name(var,"CIlower")
  mySEupper = make_name(var,"SEupper")
  mySElower = make_name(var,"SElower")
  #mySDupper = make_name(var,"SDupper")
  #mySlower = make_name(var,"SDlower")
  
  #:= is essential for any !! use, otherwise there is a parse error. Dunno why
  #!! means use the name/symbol of the variable.
  out <- data %>% group_by(obstacleoffset, obstaclecolour, sectiontype) %>% summarize(!!mymean := mean(!!var), #lookahead
                                                                                    !!myCIupper := !!mymean + (qt(0.975,df=n()-1)*sd(!!var)/sqrt(n())),
                                                                                    !!myCIlower := !!mymean - (qt(0.975,df=n()-1)*sd(!!var)/sqrt(n())),
                                                                                    #!!myCIupper := !!mymean + (qt(0.975,df=n()-1)*sd(!!var)),
                                                                                    #!!myCIlower := !!mymean - (qt(0.975,df=n()-1)*sd(!!var)),
                                                                                    !!mySEupper := !!mymean + (sd(!!var)/sqrt(n())),
                                                                                    !!mySElower := !!mymean  - (sd(!!var)/sqrt(n())),
                                                                                    #!!mySEupper = make_name(var,"SEupper")
                                                                                    #mySElower = make_name(var,"SElower")
                                                                                    condition = condition[1])
  if (drop) out[, c('obstacleoffset','obstaclecolour','sectiontype','condition')] <- list(NULL)
  return(out)        
}


GazeSummaries <- function(data, level = 'trial') {
  #Takes the raw frame by frame data and returns trial, participant, or group summaries.
  #level accepts 'trial', 'participant','group'.
  #'trial' returns trial averages. participant returns condition averages for each participant. 'group' returns condition averages averaged over participant.
  # Data should be the original data frame with all variables.
  
  L <- 99
  if (level == 'trial') L = 1
  if (level == 'participant') L = 2
  if (level == 'group') L = 3
  
  if( L == 99) stop('level must be trial, participant, or group')
  
  if (L >= 1){
    WholeTrial_avgs <- data %>% group_by(trialcode) %>% summarize(gazebias = mean(gazebias), #avg gazebias
                                                                  IQR_v = IQR(vangle), IQR_h = IQR(hangle), #IQR of gaze
                                                                  v_med = median(vangle), h_med = median(hangle), #median of gaze angle
                                                                  v_var = sd(vangle), h_var = sd(hangle), #SDs of gaze
                                                                  lookahead_mean = mean(lookahead), #lookahead
                                                                  lookahead_med = median(lookahead), #lookahead median
                                                                  obs = n(), #amount of data
                                                                  perc_on_srf = (sum(on_srf) / obs) * 100,
                                                                  ID = unique(ID),
                                                                  condition = unique(condition), count = unique(count), block = unique(block), #condition & count & block
                                                                  sectionorder = sectionorder[2], sectiontype = sectiontype[2], #sectiondata  
                                                                  obstacleoffset = obstacleoffset[2], obstaclecolour = obstaclecolour[2]) #obstacle info
  
    
    output <- WholeTrial_avgs
    if (L >= 2){
      #do participant condition avgs
      
      PPCndt_avgs <- WholeTrial_avgs  %>% group_by(condition, sectiontype, ID, block) %>% summarize(gazebias = mean(gazebias), 
                                                                                IQR_v = mean(IQR_v), IQR_h= mean(IQR_h), #mean IQRs
                                                                                lookahead_mean = mean(lookahead), #mean lookahead
                                                                                lookahead_med = mean(lookahead), #mean lookahead
                                                                                h_med = mean(h_med), v_med = mean(v_med), #mean of gaze medians
                                                                                h_var = mean(h_var), v_var = mean(v_var), #mean of gaze SDs.
                                                                                perc_on_srf = mean(perc_on_srf), #percentage on srf
                                                                                obstacleoffset = unique(obstacleoffset), obstaclecolour = unique(obstaclecolour)) #carry through obstacle info
      
      output <- PPCndt_avgs
                
      if (L == 3){
        #do group averages with CIs
        #is there are way to put this into a loop, with different variable names using paste and assign() 
        
        lookahead_mean <- CalculateSummaries(PPCndt_avgs, 'lookahead_mean', drop = FALSE)
        lookahead_med <- CalculateSummaries(PPCndt_avgs, 'lookahead_med', drop = FALSE)
        IQR_v <- CalculateSummaries(PPCndt_avgs, 'IQR_v', drop = TRUE)
        IQR_h <- CalculateSummaries(PPCndt_avgs, 'IQR_h', drop = TRUE)
        v_var <- CalculateSummaries(PPCndt_avgs, 'v_var', drop = TRUE)
        h_var <- CalculateSummaries(PPCndt_avgs, 'h_var', drop = TRUE)
        h_var <- CalculateSummaries(PPCndt_avgs, 'perc_on_srf', drop = TRUE)
        
        GroupCndts_avgs <- cbind(lookahead, IQR_v, IQR_h, v_var, h_var, perc_on_srf)
        
        output = GroupCndts_avgs
      }
    }
  }
  return(output)
}

make_colname <- function(var, suffix)(paste(var,suffix, sep="."))

PlotPointRange <- function(data, var, ylabel, ylims = NULL){
  # Plot range plot depending on var
  #ylabel = ylab
  
  #make strings
  mymean = make_colname(var,"mean")
  myCIupper = make_colname(var,"CIupper")
  myCIlower = make_colname(var,"CIlower")
  
  #mySEupper = make_colname(var,"SEupper")
  #mySElower = make_colname(var,"SElower")
  
  #make plot
  myplot <- ggplot(data=data, aes_string(x='sectiontype', y=mymean, ymin=myCIlower, ymax=myCIupper, col='condition')) +
    geom_pointrange(position = position_dodge(width=.3), fatten=2, size=.8) + 
    xlab("Section Type") + ylab(ylabel) + 
    theme_transition + 
    scale_colour_manual(values = mypalette, labels = c("0"="Attract_Narrow", "1"="Attract_Wide", "2"="Avoid_Narrow", "3"="Avoid_Wide"), name = "Condition") +
    theme(panel.grid.major.x = element_line(color="grey85",size=.2, linetype = 2), strip.text = element_text(face="plain",colour="black",size="7"), plot.title= element_text(face="bold",colour="black",size="8", hjust=.5)) +
    #coord_cartesian(ylim = c(-3,2)) +
    scale_x_discrete(labels = c("Manual","Playback","Stock","PID","Interp"))
  
  if (!is.null(ylims)){
    print("True")
    myplot <- myplot + coord_cartesian(ylim = ylims)
  } 
  
  if (is.null(ylims)) print("NULL")
  
  return(myplot)
  
  
}

add_cowplottitle <- function(cp, titlestr){
  #takes a cowplot plot 'cp', adds a title 'titlestr'
  title <- ggdraw() + draw_label(titlestr, fontface='bold')
  plot_grid(title, cp, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
}
