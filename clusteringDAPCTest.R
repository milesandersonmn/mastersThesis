library(adegenet)
library("mmod")
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

setwd("~/mastersThesis/")
?read.genepop

obj1 <- read.genepop("pruned.thinnedSNPs-1.gen", ncode = 3)
?find.clusters
grp <- find.clusters(obj1, method = "ward", max.n.clust = 40, choose.n.clust =TRUE, criterion = "diffNgroup" )
200
7
grp$Kstat
?Gst_Hedrick
Gst_Hedrick(obj1)

sampleClusters <- obj1[c(cluster3individuals,cluster2individuals)]

distance <- provesti.dist(sampleClusters)
theTree <- distance %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

aboot(sampleClusters, dist = nei.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

dapc1 <- dapc(obj1, grp$grp)
21
6

cluster7individuals <- attr(dapc1$grp[dapc1$grp==7], which = "names")
cluster2individuals <- attr(dapc1$grp[dapc1$grp==2], which = "names")
cluster3individuals <- attr(dapc1$grp[dapc1$grp==3], which = "names")
cluster4individuals <- attr(dapc1$grp[dapc1$grp==4], which = "names")
obj1@tab[c(cluster7individuals,"LA4329_NO2_LA4329_NO2"),contrib1]

scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")

summary(dapc1)
assignplot(dapc1, subset = 175:200)
dapc1$grp==3
subset(dapc1, dapc1$grp==3)
dapc1[dapc1$grp == levels(dapc1$grp)[3]]
attr(dapc1$grp[dapc1$grp==3], which = "names")
dapc1$posterior

?dapc()

dapc2 <- dapc(obj1, n.da = 50, n.pca = 100)

optim.a.score(dapc1)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc1, lab="", posi=list(x=12,y=-.01), cleg=.7)
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

loadingplot(dapc1$var.contr,axis =1,
            thres=.004, lab.jitter=1)

contrib1 <- attr(dapc1$var.contr[,1][dapc1$var.contr[,2] > 0.004], which = "names")

str(cluster2individuals)
