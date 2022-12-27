library(adegenet)
library("mmod")
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

setwd("~/mastersThesis/")

obj2 <- read.genepop("intergenicSNPs.noSingletons.noOutgroup.gen", ncode = 3)

popName <- as.data.frame(as.character(obj2@pop))
popCoordinates <- cbind(popName, "X", "Y")
str(popCoordinates)

popInfo <- as.data.frame(read.table(file = "selected_pops_genepopNames.txt", sep = "\t"))

for (x in popCoordinates[1:158,1]) {
  
  spatialCoordinates <- popInfo[which(popInfo$V1==x),2:3]
  popCoordinates[which(popCoordinates==x),2:3] <- spatialCoordinates
  
  
}
str(unique(popCoordinates))
str(popCoordinates)
obj2@other$xy <- as.matrix(unique(popCoordinates))

obj2Subset <- subset(obj2, obj2@pop!="LA4339_NO1")
obj2Subset2 <- subset(obj2Subset, obj2Subset@pop!="LA2880_N17")

obj2Subset2.smry <- summary(obj2Subset2)

plot(obj2Subset2.smry$Hexp, obj2Subset2.smry$Hobs, main = "Observed vs expected heterozygosity")
abline(0, 1, col = "red")

t.test(obj2Subset2.smry$Hexp, obj2Subset2.smry$Hobs, paired = TRUE, var.equal = TRUE)

toto <- genind2genpop(obj2Subset2)
Dgen <- dist.genpop(toto, method = 2)
Dgeo <- dist(obj2Subset2$other$xy[1:16,], method = "euclidian")
dim(as.matrix(Dgen))
dim(as.matrix(Dgeo))
str(Dgeo)
str(Dgen)
ibd <- mantel.randtest(Dgen, Dgeo)
plot(ibd)

plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red", lty=2)

library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")

library(graph4lg)

pairwiseFST <- mat_pw_fst(obj2Subset2)
str(pairwiseFST)
summary(pairwiseFST)
boxplot(pairwiseFST)
boxplot(pairwiseFSTSI)
