if (!require("RColorBrewer")) { install.packages("RColorBrewer") }
library(RColorBrewer)
library(tidyverse)
setwd("~/mastersThesis/FST/")
data1 <- read.table("pairwisePopulations.FST", header = TRUE, sep = "\t")



mycolors = brewer.pal(name="Dark2", n = 8)

q <- ggplot2::ggplot() +
  geom_boxplot(data = data1, aes(x = reorder(data1[,1], data1[,5]), y = data1[,5], color = data1[,8]))



q + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

######Valley1 subset##########

Valley1 <- data1[data1$region1=="Valley1" | data1$region2=="Valley1",]

valley1plot <- ggplot2::ggplot() +
  geom_boxplot(data = Valley1, aes(x = reorder(Valley1[,1], Valley1[,5]), y = Valley1[,5], color = Valley1[,8])) +
  scale_color_manual(values = mycolors)

valley1plot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

############Valley2 subset#############

Valley2 <- data1[data1$V6=="Valley2" | data1$V7=="Valley2",]

Valley2plot <- ggplot2::ggplot() +
  geom_boxplot(data = Valley2, aes(x = reorder(Valley2[,1], Valley2[,5]), y = Valley2[,5], color = Valley2[,8])) +
  scale_color_manual(values = mycolors)

Valley2plot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

############Valley3 subset#############

Valley3 <- data1[data1$V6=="Valley3" | data1$V7=="Valley3",]

Valley3plot <- ggplot2::ggplot() +
  geom_boxplot(data = Valley3, aes(x = reorder(Valley3[,1], Valley3[,5]), y = Valley3[,5], color = Valley3[,8])) +
  scale_color_manual(values = mycolors)

Valley3plot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

############Valley4 subset#############

Valley4 <- data1[data1$V6=="Valley4" | data1$V7=="Valley4",]

Valley4plot <- ggplot2::ggplot() +
  geom_boxplot(data = Valley4, aes(x = reorder(Valley4[,1], Valley4[,5]), y = Valley4[,5], color = Valley4[,8])) +
  scale_color_manual(values = mycolors)

Valley4plot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

#########Within groups subset###############

withinGroup <- data1[data1$grouping=="Within group",]

withinGroupPlot <- ggplot2::ggplot() +
  geom_boxplot(data = withinGroup, aes(x = reorder(withinGroup[,1], withinGroup[,6]), y = withinGroup[,5], color = withinGroup[,6])) 


withinGroupPlot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

anova <- aov(withinGroup$V5 ~ withinGroup$V6)
summary(anova)

tukey <- TukeyHSD(anova)
tukey

############Within group excluded###############

outsideGroup <- data1[data1$V8!="Within group",]

outsideGroupPlot <- ggplot2::ggplot(data = outsideGroup, aes(x = reorder(outsideGroup[,1], outsideGroup[,8]), y = outsideGroup[,5], color = outsideGroup[,8])) +
  geom_boxplot() +
  scale_color_manual(values = mycolors)


outsideGroupPlot + theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 5))

anova <- aov(outsideGroup$V5 ~ outsideGroup$V8)
summary(anova)

par(mfrow=c(2,2))
plot(anova)
par(mfrow=c(1,1))
