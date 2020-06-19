data <- read.csv("D:/Orca19_FullSteeringDataset/Orca19_collated_steering2.csv")
setwd("C:/git_repos/Orca18_Analysis/Post-Processing/")
saveRDS(data, "../Data/Orca19_collated_steering2.rds")
