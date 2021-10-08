
setwd("../Becher2022data/perSuperClusterJoins/")

# This script requires functions defined and objects computed in 01analyses.R!

# read in data ####
clJoinFiles <- dir(pattern = "*/*Counts", recursive = T)

clJoins = list()
length(clJoins)
str(clJoins)

for(i in clJoinFiles){
  if(length(readLines(i)) == 0){
    clJoins[[substr(i, 6, 13)]] <- data.frame(0,0,0)
  } else {
  clJoins[[substr(i, 6, 13)]] <- read.table(i, header = F)
  }
  names(clJoins[[substr(i, 6, 13)]]) <- c("count", "E030", substr(i, 6, 13))
}


clJMats <- list()
for(i in names(clJoins)){
  clJMats[[i]] <- matrix(0, nrow=max(clJoins[[i]][,2]) + 1, ncol=max(clJoins[[i]][,3]) + 1)

  for(j in 1:nrow(clJoins[[i]])){
    clJMats[[i]][clJoins[[i]][j, 2] + 1, clJoins[[i]][j, 3] + 1] <- clJoins[[i]][j, 1]
  }
}
# plotMatCont(clJMats[[1]])
# plotMatCont(clJMats[[2]])
# plotMatCont(clJMats[[3]])
# plotMatCont(clJMats[[4]])
# plotMatCont(clJMats[[5]])
# plotMatCont(clJMats[[51]])
# plotMatCont(clJMats[[2]])
# plotMatCont(clJMats[[3]])
# plotMatCont(clJMats[[4]])
# plotMatCont(clJMats[[155]]) # an2
# plotMatCont(clJMats[[5]]) # vi
# plotMatCont(clJMats[[104]])
# plotMatCont(clJMats[[255]])
# plotMatCont(clJMats[[55]]) # ri1
# plotMatCont(clJMats[[205]]) # ri3
# plotMatCont(clJMats[[255]]) # ri3
# plotMatCont(clJMats[[105]]) # ro
# 
# plotMatCont(clJMats[[204]]) # ri2
# 
# 
# plotMatCont(clJMats[[101]]) # Angela
# plotMatCont(clJMats[[52]]) # sat
# plotMatCont(clJMats[[104]]) # 45S rDNA
# 
# 
# 
# plotMatCont(clJMats[[57]])


# contribution of clusters to a sample's GS

# gsFromBinJoin(clJMats[[1]])
# gsFromBinJoin(clJMats[[101]])
# gsFromBinJoin(clJMats[[2]])
# gsFromBinJoin(clJMats[[3]])
# gsFromBinJoin(clJMats[[4]])
# gsFromBinJoin(clJMats[[5]])
# gsFromBinJoin(clJMats[[7]])
# gsFromBinJoin(clJMats[[57]])
# gsFromBinJoin(clJMats[[2]])
# gsFromBinJoin(clJMats[[102]])
# gsFromBinJoin(clJMats[[29]])

clDiffs <- sapply(clJMats, function(x) {
  a <- gsFromBinJoin(x)
  a[1]-a[2]
})


# contribution to GS diff for the top 50 repeat super clusters
plot(clDiffs[1:50], ylim=c(-75,75)) # Vi-An1
points(clDiffs[51:100], col=2)      # Ri1-An1
points(clDiffs[101:150], col=3)     # Ro-An1
points(clDiffs[151:200], col=4)     # An2-An1
points(clDiffs[201:250], col=5)     # Ri2-An1
points(clDiffs[251:300], col=6)     # Ri3-An1
abline(v=(1:10)*5)

diffSuVi <- sum(clDiffs[1:50], na.rm = T) # Vi
diffSuRi1 <-sum(clDiffs[51:100], na.rm = T) # Ri1
diffSuRo <-sum(clDiffs[101:150], na.rm = T) #Ro
diffSuAn2 <-sum(clDiffs[151:200], na.rm = T) #An2
diffSuRi2 <-sum(clDiffs[201:250], na.rm = T) # Ri2
diffSuRi3 <-sum(clDiffs[251:300], na.rm = T) # Ri3


clCOntribs <- sapply(clJMats, function(x) {
  a <- gsFromBinJoin(x)
  a[1]
})
plot(clCOntribs)
plot(clCOntribsE030)
clCOntribsE030 <- sapply(clJMats[1:50], function(x) {
  a <- gsFromBinJoin(x)
  a[2]
})
# plot(clCOntribs)
# plot(clCOntribsE030)

# Which super clusters show the largest variation in their contribution to GS?
varCoef <- function(x) sd(x)/mean(x) # coefficient of variation, can be high if difference are small in magnitude
varCoefs <- apply(matrix(c(clCOntribsE030, clCOntribs), ncol=7), 1, varCoef)
sds <- apply(matrix(c(clCOntribsE030, clCOntribs), ncol=7), 1, sd) # standard dev
plot(varCoefs)
abline(v=(1:10)*5)
plot(sds)
abline(v=(1:10)*5)

plot(apply(sapply((0:5) * 50, function(x) (clDiffs[(1:50) + x])), 1, sd))
suSd <- apply(sapply((0:5) * 50, function(x) (clDiffs[(1:50) + x])), 1, sd)
suSd[order(suSd, decreasing = T)[1:3]]


# difference between size difference and difference in contibution of first 50 super clusters
cumsumDiff <- function(jm){
  -sum(jm[,4] * bin2lin(jm[,1]))/1000000
}

# requires object compDFs from 01analyses.R
gsDiffsCumsum <- sapply(c("E031", "E032", "E040", "E065", "E068", "E073"), function(x) cumsumDiff(compDfs[[x]]))
c(diffSuVi, diffSuRi1, diffSuRo, diffSuAn2, diffSuRi2, diffSuRi3) / -gsDiffsCumsum




# contribution of super cluster to GS in copy number range ####
compDfs <- lapply(jMats, specComp)

contByOom(specComp(clJMats[[101]]))

diffsByOom(specComp(clJMats[[101]]))


contByOom(specComp(clJMats[[102]]))
diffsByOom(specComp(clJMats[[102]]))

contByOom(specComp(clJMats[[104]]))
diffsByOom(specComp(clJMats[[104]]))

contOomSu <- function(su){
  a <- sapply(c(0:5)*50 + su, function(x) contByOom(specComp(clJMats[[x]])))
  dff <- specComp(clJMats[[su]])
  dff[,2] <- dff[,3]
  nms <- c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  br <- lin2bin(c(0.47, 10, 100, 1000, 10000, 100000, 1000000, 10000000))
  nums030 <- sapply(1:(length(br)-1), function(x){
      #print(x)
      a <- dff[(dff[,1] > br[x]) & (dff[,1] <= br[x+1]), ]

      b <- sum(bin2lin(a[,1]) * a[,2])/1000000
      names(b) <- paste0(nms[x],"-", nms[x+1])
      b
    })
  c <- rbind(nums030,
             t(a))
  #rownames(c) <- c("E030", "E031", "E032", "E040", "E065", "E068", "E073")
  rownames(c) <- c("An1", "Vi", "Ri1", "Ro", "An2", "Ri2", "Ri3")
  c
  }
contOomSu(1)

#pdf("ContOrdersOfMagSu.pdf", width=6, height=4)
barplot(t(contPerOom[c(1, 5, 2, 4, 3, 6, 7),]),
        beside =T,
        names.arg = c("An1", "An2", "Vi", "Ro", "Ri1", "Ri2", "Ri3"),
        ylab = "Contibution to genome size difference (Mbp)",
        xlab="Sample")
nms <- c("0", "10", "100", "1k", "10k", "100k", "1M", "10M")
legend("topleft",
       fill=grey.colors(7),
       legend=sapply(1:7, function(x) paste0(nms[x],"-", nms[x+1])),
       title="Copy number range"
)

#t(sapply(c(0,0:5)*50 + 1, function(x) contByOom(specComp(clJMats[[x]]))))
barplot(t(contOomSu(1)[c(1, 5, 2, 4, 3, 6, 7),]),
        beside = T,
        add = T,
        col=2)
barplot(t(contOomSu(4)[c(1, 5, 2, 4, 3, 6, 7),]),
        beside = T,
        add = T,
        col=3)
barplot(t(contOomSu(2)[c(1, 5, 2, 4, 3, 6, 7),]),
        beside = T,
        add = T,
        col=4)

#dev.off()


#pdf("OrdersOfMag.pdf", width=6, height=4)
diffOomSu <- function(su){
  a <- t(sapply(c(0:5)*50 + su, function(x) diffsByOom(specComp(clJMats[[x]]))))

  rownames(a) <- c("Vi", "Ri1", "Ro", "An2", "Ri2", "Ri3")
  a
}
diffOomSu(1)
diffOomSu(4)


barplot(t(perOom[c(4, 1, 3, 2, 5, 6),]),
        beside =T,
        names.arg = c("An2", "Vi", "Ro", "Ri1", "Ri2", "Ri3"),
        ylab = "Contibution to genome size difference (Mbp)",
        xlab="Sample")
nms <- c("0", "10", "100", "1k", "10k", "100k", "1M", "10M")
legend("topleft",
       fill=grey.colors(7),
       legend=sapply(1:7, function(x) paste0(nms[x],"-", nms[x+1])),
       title="Copy number range"
)
barplot(t(diffOomSu(1)[c(4, 1, 3, 2, 5, 6),]),
        beside = T,
        add = T,
        col = 2)
barplot(t(diffOomSu(4)[c(4, 1, 3, 2, 5, 6),]),
        beside = T,
        add = T,
        col = 3)
barplot(t(diffOomSu(2)[c(4, 1, 3, 2, 5, 6),]),
        beside = T,
        add = T,
        col = 4)
#dev.off()
