setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv("../Data/collated_steering.csv")
saveRDS(data, "../Data/collated_steering.rds")
