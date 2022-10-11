library(adegenet)
library("mmod")
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

setwd("~/mastersThesis/")

objSI2 <- read.genepop("intergenicSNPS.noOutgroup.rename.noOutliers.gen", ncode = 3)

popName <- as.data.frame(as.character(objSI@pop))
popCoordinates <- cbind(popName, "X", "Y")
str(popCoordinates)

popInfo <- as.data.frame(read.table(file = "selected_pops_genepopNames.txt", sep = "\t"))

for (x in popCoordinates[,1]) {
  
  spatialCoordinates <- popInfo[which(popInfo$V1==x),2:3]
  popCoordinates[which(popCoordinates==x),2:3] <- spatialCoordinates
  
  
}
str(unique(popCoordinates))
str(popCoordinates)
objSI@other$xy <- as.matrix(unique(popCoordinates))

objSISubset <- subset(objSI, objSI@pop!="LA4339_NO1")
objSISubset2 <- subset(objSISubset, objSISubset@pop!="LA2880_N17")

objSISubset2.smry <- summary(objSISubset2)

plot(objSISubset2.smry$Hexp, objSISubset2.smry$Hobs, main = "Observed vs expected heterozygosity")
abline(0, 1, col = "red")

t.test(objSISubset2.smry$Hexp, objSISubset2.smry$Hobs, paired = TRUE, var.equal = TRUE)

totoSI <- genind2genpop(objSISubset2)
DgenSI <- dist.genpop(totoSI, method = 2)
DgeoSI <- dist(objSISubset2$other$xy[1:16,], method = "euclidian")
dim(as.matrix(Dgen))
dim(as.matrix(Dgeo))
str(Dgeo)
str(Dgen)
ibdSI <- mantel.randtest(DgenSI, DgeoSI)
plot(ibdSI)

plot(DgeoSI, DgenSI)
abline(lm(DgenSI~DgeoSI), col="red", lty=2)

library(MASS)
dens <- kde2d(DgeoSI,DgenSI, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(DgeoSI, DgenSI, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(DgenSI~DgeoSI))
title("Isolation by distance plot")

pairwiseFSTSI <- mat_pw_fst(objSISubset2)
str(pairwiseFSTSI)
summary(pairwiseFSTSI)
boxplot(pairwiseFSTSI)
