library(tidyverse)
library(reshape2)
setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/05_Drift")
dataA <- read.table("A_drift.txt",header=F,stringsAsFactors = F)
dataA$group <- "A"
dataB <- read.table("B_drift.txt",header=F,stringsAsFactors = F)
dataB$group <- "B"
dataD <- read.table("D_drift.txt",header=F,stringsAsFactors = F)
dataD$group <- "D"
all<- as.data.frame(rbind(dataA,dataB,dataD))
all$id <- c(1:64)
colnames(all) <- c("individual","value","group","id")
data<-all
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
pi <- read.table("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/06_PI/Pi_stat.txt",header=T,stringsAsFactors = F,sep="\t")
data$pi <- pi$value
data$sd <- pi$sd
label_data$value[1] <- label_data$value[1]
label_data$value[26] <- label_data$value[26]
label_data$value[27] <- label_data$value[27]
pd <- position_dodge(0.1)
data$value <- data$value/10
data$pi <- -data$pi
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.012, xend = start, yend = 0.012), colour = "grey", alpha=1, size=0.8 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.009, xend = start, yend = 0.009), colour = "grey", alpha=1, size=0.8, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.006, xend = start, yend = 0.006), colour = "grey", alpha=1, size=0.8, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.003, xend = start, yend = 0.003), colour = "grey", alpha=1, size=0.8 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.003, 0.006, 0.009, 0.012), label = c("0.003", "0.006", "0.009", "0.012") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) + 
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-0.015,0.02) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value/10+0.001, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=label_data, aes(x=id, y=value/10+0.001, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  #geom_text(data=base_data, aes(x = title, y = -0.005, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)+
  geom_errorbar(aes(ymin=pi-sd, ymax=pi+sd), width=.1, position=pd) +
  geom_point(data=data, aes(x = id, y = pi), size=2)
p



