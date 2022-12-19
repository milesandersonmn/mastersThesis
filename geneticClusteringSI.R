library(adegenet)
library("mmod")
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

setwd("~/mastersThesis/")

objSI <- read.genepop("pruned.SIgenesSNPs.noOutgroup.rename.gen", ncode = 3)

grpSI <- find.clusters(objSI, method = "ward", max.n.clust = 40, choose.n.clust =TRUE, criterion = "diffNgroup" )
200
5

dapcSI <- dapc(objSI, grpSI$grp)
21
5


summary(dapcSI)
assignplot(dapcSI)

scatter(dapcSI)

attr(dapcSI$assign[dapcSI$assign==5], which = "names")
dapcSI$assign

cluster1individualsSI <- attr(dapcSI$grp[dapcSI$grp==1], which = "names")
cluster2individualsSI <- attr(dapcSI$grp[dapcSI$grp==2], which = "names")
cluster3individualsSI <- attr(dapcSI$grp[dapcSI$grp==3], which = "names")
cluster4individualsSI <- attr(dapcSI$grp[dapcSI$grp==4], which = "names")
cluster5individualsSI <- attr(dapcSI$grp[dapcSI$grp==5], which = "names")




sampleClusters <- objSI[c(cluster4individuals,cluster5individuals)]

distance <- provesti.dist(sampleClusters)
theTree <- distance %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.
?aboot
aboot(sampleClusters, dist = nei.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

optim.a.score(dapcSI)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapcSI, lab="", posi=list(x=12,y=-.01), cleg=.7)
temp <- which(apply(dapcSI$posterior,1, function(e) all(e<0.9)))
temp

loadingplot(dapcSI$var.contr,axis =1,
            thres=.004, lab.jitter=1)

contrib1 <- attr(dapcSI$var.contr[,1][dapcSI$var.contr[,2] > 0.004], which = "names")

str(cluster2individuals)
