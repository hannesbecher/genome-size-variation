# data ####
library(Tetmer)

setwd("~/Dropbox/manuscripts/2108_kmers_gs/data/analyseSpectraFastp/")
samples <- dir(pattern = "PlMt.hist")
samplesO <- dir(pattern = "full.hist")
spectra = list()
spectraO = list()
for(i in samples){
  spectra[[substr(i,1,4)]] <- read.spectrum(i, substr(i,1,4), 21)
}

for(i in samplesO){
  spectraO[[substr(i,1,4)]] <- read.spectrum(i, substr(i,1,4), 21)
}


names(spectra)
str(spectra['E030'])

# Tetmer ####

tetmer(spectra[["E030"]])
tetmer(spectra[["E031"]])
tetmer(spectra[["E032"]])
tetmer(spectra[["E040"]])
tetmer(spectra[["E065"]])
tetmer(spectra[["E068"]])
tetmer(spectra[["E073"]])


# GS from uncropped spectrum ####
gsFromSpectrum <- function(spec, dep){
  s <- spec@data[spec@data[,1] > dep*0.55,]
  sum(s[,1] * s[,2])/dep/1000000
}
head(spectra[["E030"]]@data)
gsE030a <- gsFromSpectrum(spectra[["E030"]], 54)
gsE031a <- gsFromSpectrum(spectra[["E031"]], 42.4)
gsE032a <- gsFromSpectrum(spectra[["E032"]], 35)
gsE040a <- gsFromSpectrum(spectra[["E040"]], 67.4)
gsE065a <- gsFromSpectrum(spectra[["E065"]], 28.5)
gsE068a <- gsFromSpectrum(spectra[["E068"]], 25.5)
gsE073a <- gsFromSpectrum(spectra[["E073"]], 20.8)



gsFromSpectrum(spectraO[["E030"]], 54)/gsFromSpectrum(spectra[["E030"]], 54)
gsFromSpectrum(spectraO[["E031"]], 42.4)/gsFromSpectrum(spectra[["E031"]], 42.4)
gsFromSpectrum(spectraO[["E032"]], 35)/gsFromSpectrum(spectra[["E032"]], 35)
gsFromSpectrum(spectraO[["E040"]], 67.4)/gsFromSpectrum(spectra[["E040"]], 67.4)
gsFromSpectrum(spectraO[["E065"]], 28.5)/gsFromSpectrum(spectra[["E065"]], 28.5)
gsFromSpectrum(spectraO[["E068"]], 25.5)/gsFromSpectrum(spectra[["E068"]], 25.5)
gsFromSpectrum(spectraO[["E073"]], 20.8)/gsFromSpectrum(spectra[["E073"]], 20.8)


# joint binned spectra
fJoins = dir(pattern = "trim")

bin2lin <- function(x, binwidth=1.1){
  a <- which(x == 0)
  b <- 0.5 * binwidth^(x-1) * 1.05
  b[a] <- 0
  return(b)
}


lin2bin <- function(x){
  return(floor(log(x/0.5)/log(1.1)) + 1)
}

joins = list()
for(i in fJoins){
  joins[[substr(i, 5, 8)]] <- read.table(i, header = F)
  names(joins[[substr(i, 5, 8)]]) <- c("count", "E030", substr(i, 5, 8))
}
head(joins[["E031"]])

jMats <- list()
for(i in names(joins)){
  jMats[[i]] <- matrix(0, nrow=max(joins[[i]][,2]) + 1, ncol=max(joins[[i]][,3]) + 1)

  for(j in 1:nrow(joins[[i]])){
    jMats[[i]][joins[[i]][j, 2] + 1, joins[[i]][j, 3] + 1] <- joins[[i]][j, 1]
  }
}

plotJMat <- function(jMat, ...){
  image(0:(nrow(jMat)-1), 0:(ncol(jMat)-1), jMat, asp=1, ...)
}
jMats[["E031"]]

contJMat <- function(jMat, ...){
  contour(0:(nrow(jMat)-1), 0:(ncol(jMat)-1), jMat, asp=1, ...)
}
plotJMat(log(jMats[["E032"]]), col = (heat.colors(50)))
abline(0,1)

plotJMat(log(jMats[["E040"]]), col = (rainbow(30)), xaxt='n', yaxt='n', ylab="E. rostkoviana (Ro)", xlab="E. anglica (An1)")
contJMat(log(jMats[["E040"]]), add=T)
contJMat(log(jMats[["E040"]]))

abline(v=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
abline(h=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
plotJMat(log(jMats[["E040"]]), col = (rainbow(30)), xaxt='n', yaxt='n', add=T)
abline(0,1)
abline(v=lin2bin(c(1, 2, 4)), col="grey", lty=2)
abline(h=lin2bin(c(1, 2, 4)), col="grey", lty=2)
ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
axis(2, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))



# example plot joint
#pdf("jointRo.pdf", width=5, height=5)
plotJMat(log(jMats[["E040"]]),
         col = "#FFFFFF",
         xaxt='n',
         yaxt='n',
         xlab="Copy number bins An",
         ylab="Copy number bins Ro")
abline(v=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
abline(h=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
plotJMat(log(jMats[["E040"]]), col = (rainbow(30)), add=T)
abline(0,1)
abline(v=lin2bin(c(1, 2, 4)), col="grey", lty=2)
abline(h=lin2bin(c(1, 2, 4)), col="grey", lty=2)
contJMat(log(jMats[["E040"]]),
         xaxt='n',
         yaxt='n',
         add=T,
         drawlabels=F)
ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
axis(2, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
#dev.off()
plotMatCont <- function(jMat,...){
  plotJMat(log(jMat),
           col = "#FFFFFF",
           xaxt='n',
           yaxt='n',
           xlab="Copy number bins An",
           ...
           )
  abline(v=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
  abline(h=lin2bin(c(10, 100, 1000, 10000, 100000, 1000000)), col="grey")
  plotJMat(log(jMat), col = (rainbow(30)), add=T)
  abline(0,1)
  abline(v=lin2bin(c(1, 2, 4)), col="grey", lty=2)
  abline(h=lin2bin(c(1, 2, 4)), col="grey", lty=2)
  contJMat(log10(jMat),
           xaxt='n',
           yaxt='n',
           add=T,
           drawlabels=F
           )
  ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
  axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
  axis(2, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
}
#pdf("jointVi.pdf", width=5, height=5)
plotMatCont(jMats[["E031"]], ylab="Copy number bins in Vi")
#dev.off()

#pdf("jointRi1.pdf", width=5, height=5)
plotMatCont(jMats[["E032"]], ylab="Copy number bins in Ri1")
#dev.off()

#pdf("jointRo.pdf", width=5, height=5)
plotMatCont(jMats[["E040"]], ylab="Copy number bins in Ro")
#dev.off()

#pdf("jointAn2.pdf", width=5, height=5)
plotMatCont(jMats[["E065"]], ylab="Copy number bins in An2")
#dev.off()

#pdf("jointRi2.pdf", width=5, height=5)
plotMatCont(jMats[["E068"]], ylab="Copy number bins in Ri2")
#dev.off()

#pdf("jointRi3.pdf", width=5, height=5)
plotMatCont(jMats[["E073"]], ylab="Copy number bins in Ri3")
#dev.off()

dim(jMats[["E040"]])

plotJMat(log(jMats[["E065"]]))


plotJMat(log(jMats[["E065"]]))

plotJMat(log(jMats[["E068"]]))


plotJMat(log(jMats[["E073"]]))
abline(0, 1)
abline(v=lin2bin(c(1, 2, 4)))
abline(h=lin2bin(c(1, 2, 4)))
gsFromBinJoin <- function(bMat){
  cs = ncol(bMat)
  rs = nrow(bMat)
  return(c(cols = sum(colSums(bMat)[2:cs] * bin2lin(1:(cs-1))) / 1000000,
           rows = sum(rowSums(bMat)[2:rs] * bin2lin(1:(rs-1))) / 1000000))
}

gsE030b <- gsFromBinJoin(jMats[["E031"]])[2]
gsE031b <- gsFromBinJoin(jMats[["E031"]])[1]
gsE032b <- gsFromBinJoin(jMats[["E032"]])[1]
gsE040b <- gsFromBinJoin(jMats[["E040"]])[1]
gsE065b <- gsFromBinJoin(jMats[["E065"]])[1]
gsE068b <- gsFromBinJoin(jMats[["E068"]])[1]
gsE073b <- gsFromBinJoin(jMats[["E073"]])[1]

plot(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a),
     c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b)
     )
c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a)/0.978
text(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a),
     c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b),
     labels = c("E030", "E031", "E032", "E040", "sE065", "E068", "E073"))
summary(lm(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a)~
           c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b),
        offset=c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b)
))
abline(0,1)
lin2bin(1)
lin2bin(0.5)
lin2bin(0.5775)

lin2bin(0.51)

bin2lin(0)
bin2lin(1)

gsDiffs <- gsE030a - c(gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a)
gsDiffsB <- gsE030b - c(gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b)
# comparing binned spectra starting from bin joins ####

#  function to produce a df from a matrix
specComp <- function(bMat){
  cs <- dim(bMat)[2]
  rs <- dim(bMat)[1]
  cDf <- data.frame(bins = 0:(cs-1), counts1 = colSums(bMat))
  rDf <- data.frame(bins = 0:(rs-1), counts2 = rowSums(bMat))
  jDf <- merge(cDf, rDf, by="bins", all=T)
  jDf[is.na(jDf)] <- 0
  jDf$diff <- jDf$counts1 - jDf$counts2
  jDf
}
sCE040 <- specComp(jMats[["E040"]])

compDfs <- lapply(jMats, specComp)
str(compDfs)

compCumsumPlot <- function(jm, main=""){
  plot(jm[,1],
       cumsum(jm[,4] * bin2lin(jm[,1]))/1000000,
       xlab="Genomic copy number",
       ylab="Cumulative GS difference (Mbp)",
       main = main,
       ylim=c(-15, 250),
       xaxt = 'n'

       )
  ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
  axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
  abline(v=lin2bin(ticks), lty=2, col="grey")
  abline(h=0, col = 'grey')
  d = sum(jm[,4] * bin2lin(jm[,1]))/1000000
  abline(h = d,  col = "grey")
  axis(2, at=d, labels=round(d,0), las=1)


}

compCumsumPlotInset <- function(jm, main=""){
  plot(jm[,1],
       cumsum(jm[,4] * bin2lin(jm[,1]))/1000000,
       xlab="Genomic copy number",
       ylab="Cumulative GS difference (Mbp)",
       main = main,
       ylim=c(-15, 50),
       xlim=c(0, lin2bin(10)),
       xaxt = 'n'

  )
  ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
  axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
  abline(v=lin2bin(ticks), lty=2, col="grey")
  abline(h=0, col = 'grey')
  d = sum(jm[,4] * bin2lin(jm[,1]))/1000000
  abline(h = d,  col = "grey")
  axis(2, at=d, labels=round(d,0), las=1)


}

cumsumDiff <- function(jm){
  -sum(jm[,4] * bin2lin(jm[,1]))/1000000
}

compPlot <- function(jm, main = ""){
  plot(jm[,1],
       (jm[,4] * bin2lin(jm[,1]))/1000000,
       xlab="Genomic copy number",
       ylab="GS difference per bin (Mbp)",
       xaxt = 'n',
       main = main)
  ticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
  axis(1, at = lin2bin(ticks), labels = c("1", "2", "4", "10", "100", "1000", "10,000", "100,000", "1,000,000"))
  abline(v=lin2bin(ticks), lty=2, col="grey")
  abline(h=0, col = 'grey')

}





compPlot2 <- function(jm, main=""){
  plot(jm[,1], (jm[,2]/jm[,3]), xlab="Bin", ylab="Prop difference", log="y", main = main)
  abline(v=lin2bin(c(1,2)), lty=2, col="grey")
  abline(h=1, col = 'grey')
  abline(h=c(1.1, 1/1.1), col = 'grey', lty=2)

}

#pdf("cumsumE031.pdf", width=6, height = 4)
pdf("cumsumE031norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E031"]], main = "Vi - An1")
dev.off()
compCumsumPlotInset(compDfs[["E031"]], main = "Vi - An1")
compPlot(compDfs[["E031"]], main = "Vi - An1")

#pdf("cumsumE032.pdf", width=6, height = 4)
pdf("cumsumE032norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E032"]], main = "Ri1 - An1")
dev.off()
compPlot(compDfs[["E032"]], main = "Ri1 - An1")

#pdf("cumsumE040.pdf", width=6, height = 4)
pdf("cumsumE040norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E040"]], main = "Ro - An1")
dev.off()
compCumsumPlotInset(compDfs[["E040"]], main = "Ro - An1")
compPlot(compDfs[["E040"]], main = "Ro - An1")

#pdf("cumsumE065.pdf", width=6, height = 4)
pdf("cumsumE065norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E065"]], main = "An2 - An1")
dev.off()
compPlot(compDfs[["E065"]], main = "An2 - An1")

#pdf("cumsumE068.pdf", width=6, height = 4)
pdf("cumsumE068norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E068"]], main = "Ri2 - An1")
dev.off()
compPlot(compDfs[["E068"]], main = "Ri2 - An1")

#pdf("cumsumE073.pdf", width=6, height = 4)
pdf("cumsumE073norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E073"]], main = "Ri3 - An1")
dev.off()
compPlot(compDfs[["E073"]], main = "Ri3 - An1")


compPlot(compDfs[["E031"]])
plot(specComp(jMats[["E031"]])[2:30,1:2])
points(specComp(jMats[["E031"]])[2:30,c(1,3)], pch=3)


compCumsumPlot(compDfs[["E032"]])

plot(specComp(jMats[["E032"]])[2:30,1:2])
points(specComp(jMats[["E032"]])[2:30,c(1,3)], pch=3)


compPlot(compDfs[["E040"]])
compCumsumPlot(compDfs[["E040"]], main = "RO - AN1")

compPlot(compDfs[["E040"]])

plot(specComp(jMats[["E040"]])[2:30,1:2])
points(specComp(jMats[["E040"]])[2:30,c(1,3)], pch=3)




compPlot(compDfs[["E065"]])
compPlot2(compDfs[["E040"]])
compPlot2(compDfs[["E065"]])

compCumsumPlot(compDfs[["E065"]])
compCumsumPlot(compDfs[["E068"]])
compCumsumPlot(compDfs[["E073"]])

plot(sCE040[,1], cumsum(sCE040[,3]))
plot(sCE040[,1], cumsum(sCE040[,4] * bin2lin(sCE040[,1])))
cumsum(sCE040[,4] * bin2lin(sCE040[,1]))
gsE040a - gsE030a
gsE040b - gsE030b
abline(v=lin2bin(c(1,2, 4, 8)))


gsDiffs # from full spectra
gsDiffsCumsum <- sapply(c("E031", "E032", "E040", "E065", "E068", "E073"), function(x) cumsumDiff(compDfs[[x]]))
plot(gsDiffs, gsDiffsCumsum)
abline(0, 1)
points(gsDiffs, gsDiffsB, col = 2)
points(gsDiffsCumsum, gsDiffsB, col = 3)





# orders of magnitude ####
# contributions
contByOom <- function(dff){
  nms <- c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  br <- lin2bin(c(0.47, 10, 100, 1000, 10000, 100000, 1000000, 10000000))
  sapply(1:(length(br)-1), function(x){
    #print(x)
    a <- dff[(dff[,1] > br[x]) & (dff[,1] <= br[x+1]), ]

    b <- sum(bin2lin(a[,1]) * a[,2])/1000000
    names(b) <- paste0(nms[x],"-", nms[x+1])
    b
  })
}
numsE030 <- {
  dff <- compDfs[[1]]
  dff[,2] <- dff[,3]
  nms <- c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  br <- lin2bin(c(0.47, 10, 100, 1000, 10000, 100000, 1000000, 10000000))
  sapply(1:(length(br)-1), function(x){
    #print(x)
    a <- dff[(dff[,1] > br[x]) & (dff[,1] <= br[x+1]), ]

    b <- sum(bin2lin(a[,1]) * a[,2])/1000000
    names(b) <- paste0(nms[x],"-", nms[x+1])
    b
  })
}
contPerOom <- t(sapply(c("E031", "E032", "E040", "E065", "E068", "E073"), function(x) contByOom(compDfs[[x]])))
contPerOom <- rbind(numsE030, contPerOom)
rownames(contPerOom)[1] <- "E030"
rowSums(contPerOom)
#pdf("ContOrdersOfMag.pdf", width=6, height=4)
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
#dev.off()

# differences
lin2bin(c(1, 10, 100, 1000, 10000, 100000, 1000000))
compDfs[["E040"]]
compDfs[["E031"]]
diffsByOom <- function(dff){
  nms <- c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  br <- lin2bin(c(0.47, 10, 100, 1000, 10000, 100000, 1000000, 10000000))
  sapply(1:(length(br)-1), function(x){
    #print(x)
    a <- dff[(dff[,1] > br[x]) & (dff[,1] <= br[x+1]), ]

    b <- sum(bin2lin(a[,1]) * a[,4])/1000000
    names(b) <- paste0(nms[x],"-", nms[x+1])
    b
  })
}
perOom <- t(sapply(c("E031", "E032", "E040", "E065", "E068", "E073"), function(x) diffsByOom(compDfs[[x]])))
#barplot(perOom, beside =T)
#pdf("OrdersOfMag.pdf", width=6, height=4)
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
#dev.off()
# comparing samples in a join ####
sliceOfJoin <- function(j, l){
  a <- j[j[,2]==l,][,c(3,1)]
  a[order(a[,1]),]
}

sliceOfJoinMean <- function(j, l){
  lin2bin(weighted.mean(x=bin2lin(j[j[,2]==l,][,3]), w=j[j[,2]==l,][,1]))

}
sliceOfJoinMean(joins[["E040"]], 50)

normalisedSlicePoints <- function(j, l){
  a <- j[j[,2]==l,][,c(3,1)]
  a[,2] <- a[,2]/sum(a[,2])
  a[order(a[,1]),]
}

plot(1:140, sapply(1:140, function(x) sliceOfJoinMean(joins[["E040"]], x)))
abline(0, 1)

plot(normalisedSlicePoints(joins[["E040"]], 5), type="l", xlim=c(0,150))
sapply(seq(10, 140, 5), function(x) points(normalisedSlicePoints(joins[["E040"]], x), type="l"))

head(joins[["E040"]])
str(joins[["E040"]])







# partitioning gs variance ####


gsDat <- data.frame(sp=c("An","An","Vi","Ro","Ri","Ri","Ri"),
                    gs=c(999.98, 989.23, 1055.93, 1227.92, 1126.64, 1096.44, 1104.84)
)
summary(lm(gs ~ sp, gsDat))
var(gsDat[,2])
a <- anova(lm(gs ~ sp, gsDat))
str(a)
a$`Pr(>F)`

resample <- function(dff, n=999){
  sapply(1:n, function(x){
    a <- dff
  a[,2] <- dff[sample(1:7),2]
  anova(lm(gs ~ sp, a))$`Pr(>F)`[1]
  })
}
ps <- resample(gsDat)
hist(ps)
quantile(ps, 0.05)

# cumsum example #####
a <- dnorm(seq(-5, 5, 0.1), 0, 0.5)
b <- dnorm(seq(-5, 5, 0.1), 0, 1)
plot(a, type = "p")
points(b)

plot(cumsum(a-b))
points(a-b)


# synthetic data ####
# generated following the Tetmer logic
probs <- expression(rbind(
  dnbinom(txmin:txmax, tkcov/tbias*1, mu = tkcov * 1),
  dnbinom(txmin:txmax, tkcov/tbias*2, mu = tkcov * 2)
))

factors <- expression(
  c(
    2*tth/(1+tth),
    1/(1+tth)
  )
)

binSpectrum <- function(sp){
  bins <- lin2bin((1:150)/25)
  bins[bins < 0] <- 0
  tapply(sp, bins, sum)
}

# 0. different heterozygosity
spec01 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1)) *
                           eval(factors, envir=list(tth=0.01)))*200*1000000

spec02 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1)) *
                    eval(factors, envir=list(tth=0.1)))*200*1000000

#pdf("schematic01.pdf", width=5, height=4)
plot(0:27,
     binSpectrum(spec01)/1000000,
     xaxt='n',
     xlab="Copy number bins",
     ylab="Number of k-mers per bin (millions)")
points(0:27, binSpectrum(spec02)/1000000, pch=3)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
legend("left",
       pch=c(3, 1),
       legend=c("high", "low "),
       title="Heterozygosity"
)
#dev.off()

#pdf("schematic02.pdf", width=5, height=4)
plot(0:27,
     cumsum((binSpectrum(spec01) - binSpectrum(spec02))*bin2lin(0:27))/1000000,
     xaxt='n',
     ylab="Difference (millions)",
     xlab="Copy number bins",
     ylim=c(-35,30),
     pch=8
)
points(0:27,
       (binSpectrum(spec01) - binSpectrum(spec02))*bin2lin(0:27)/1000000,
       pch=2)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
abline(h=0, col="grey")
legend("topleft",
       pch=c(2, 8),
       legend=c("as is", "cumulative"),
       title="Weighted difference"
)
#dev.off()


# 1. different peak width
spec11 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1)) *
                    eval(factors, envir=list(tth=0.1)))*200*1000000
spec12 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1.5)) *
                    eval(factors, envir=list(tth=0.1)))*200*1000000

#pdf("schematic11.pdf", width=5, height=4)
plot(0:27,
     binSpectrum(spec11)/1000000,
     xaxt='n',
     xlab="Copy number bins",
     ylab="Number of k-mers per bin (millions)")
points(0:27, binSpectrum(spec12)/1000000, pch=3)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
legend("left",
       pch=c(3, 1),
       legend=c("high", "low "),
       title="Peak width"
)
#dev.off()

#pdf("schematic12.pdf", width=5, height=4)
plot(0:27,
     cumsum((binSpectrum(spec11) - binSpectrum(spec12))*bin2lin(0:27))/1000000,
     xaxt='n',
     ylab="Difference (millions)",
     xlab="Copy number bins",
     ylim=c(-10,15),
     pch=8
)
points(0:27,
       (binSpectrum(spec11) - binSpectrum(spec12))*bin2lin(0:27)/1000000,
       pch=2)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
abline(h=0, col="grey")
legend("topleft",
       pch=c(2, 8),
       legend=c("as is", "cumulative"),
       title="Weighted difference"
)
#dev.off()

# 2. partial duplication
spec21 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1)) *
                    eval(factors, envir=list(tth=0.01)))*200*1000000 +
  colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=50, tbias=1)) *
            eval(factors, envir=list(tth=0)))*20*1000000
spec22 <- colSums(eval(probs, envir = list(txmin= 1, txmax=150, tkcov=25, tbias=1)) *
                    eval(factors, envir=list(tth=0.01)))*200*1000000
#pdf("schematic21.pdf", width=5, height=4)
plot(0:27,
     binSpectrum(spec21)/1000000,
     xaxt='n',
     xlab="Copy number bins",
     ylab="Number of k-mers per bin (millions)")
points(0:27, binSpectrum(spec22)/1000000, pch=3)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
legend("left",
      pch=c(1, 3),
      legend=c("With duplication", "Without duplication")
      )
#dev.off()

#pdf("schematic22.pdf", width=5, height=4)
plot(0:27,
     cumsum((binSpectrum(spec21) - binSpectrum(spec22))*bin2lin(0:27))/1000000,
     xaxt='n',
     ylab="Difference (millions)",
     xlab="Copy number bins",
     pch=8
     )
points(0:27,
       (binSpectrum(spec21) - binSpectrum(spec22))*bin2lin(0:27)/1000000,
       pch=2)
axis(1, at = lin2bin(xticks), labels = xticks)
abline(v=lin2bin(c(1, 2, 4)), lty=2, col='grey')
abline(h=80, col="grey")
legend("left",
       pch=c(2, 8),
       legend=c("as is", "cumulative"),
       title="Weighted difference"
)
#dev.off()
# two diploids, different heterozygosity


