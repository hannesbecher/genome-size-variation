# Prepare ####
# download and install Tetmer:
# https://github.com/hannesbecher/shiny-k-mers/releases/download/v2.0.0-beta/Tetmer_2.0.0.tar.gz
# install via: install.packages("path/to/Tetmer_2.0.0.tar.gz", repo=NULL)

# This file uses the lab-internal E-number notation.
# The corresponding names in Becher et al. 2022 are
# E030 - An1
# E031 - Vi
# E032 - Ri1
# E040 - Ro
# E065 - An2
# E068 - Ri2
# E073 - Ri3


# data ####
library(Tetmer)

# set path to directory with your k-mer spectra (and binned joint dumps)
setwd("../Becher2022data/")

samples <- dir(pattern = "PlMt.hist") # plastid/mito seqs removed
samplesO <- dir(pattern = "full.hist") # nothing removed
spectra = list()
spectraO = list()

# Loop over names and use the Tetmer function "read.spectrum" to read the k-mer spectra
# Try out ?read.spectrum ! You can specify a sample name and a k-mer size
for(i in samples){
  spectra[[substr(i,1,4)]] <- read.spectrum(i, substr(i,1,4), 21)
}

for(i in samplesO){
  spectraO[[substr(i,1,4)]] <- read.spectrum(i, substr(i,1,4), 21)
}


names(spectra)
str(spectra['E030'])

# Tetmer ####
# Use Tetmer to estimate monoploid k-mer coverage values:
# tetmer(spectra[["E030"]])
# tetmer(spectra[["E031"]])
# tetmer(spectra[["E032"]])
# tetmer(spectra[["E040"]])
# tetmer(spectra[["E065"]])
# tetmer(spectra[["E068"]])
# tetmer(spectra[["E073"]])


# Estimate GS ####
# From uncropped spectra 

# Function to return GS in Mbp. Takes a Tetmer spectrum object (as created by
# read.spectrum()) `spec` and a monoploid multiplicity value `dep`.
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



# compute size ratios between clean and unfiltered spectra
gsFromSpectrum(spectraO[["E030"]], 54)/gsFromSpectrum(spectra[["E030"]], 54)
gsFromSpectrum(spectraO[["E031"]], 42.4)/gsFromSpectrum(spectra[["E031"]], 42.4)
gsFromSpectrum(spectraO[["E032"]], 35)/gsFromSpectrum(spectra[["E032"]], 35)
gsFromSpectrum(spectraO[["E040"]], 67.4)/gsFromSpectrum(spectra[["E040"]], 67.4)
gsFromSpectrum(spectraO[["E065"]], 28.5)/gsFromSpectrum(spectra[["E065"]], 28.5)
gsFromSpectrum(spectraO[["E068"]], 25.5)/gsFromSpectrum(spectra[["E068"]], 25.5)
gsFromSpectrum(spectraO[["E073"]], 20.8)/gsFromSpectrum(spectra[["E073"]], 20.8)
# Dirty estimates tend to be 3.8-7.2% higher than clean ones!

# joint binned spectra
fJoins = dir(pattern = "trim")

# expected multiplicity form bin number
bin2lin <- function(x, binwidth=1.1){
  a <- which(x == 0)
  b <- 0.5 * binwidth^(x-1) * 1.05
  b[a] <- 0
  return(b)
}
# bin2lin(8)
# bin2lin(100)

# turn multiplicity into bin
lin2bin <- function(x){
  return(floor(log(x/0.5)/log(1.1)) + 1)
}
# lin2bin(1)
# lin2bin(100)

# loop over binned joint k-mer dump files and import
joins = list()
for(i in fJoins){
  joins[[substr(i, 5, 8)]] <- read.table(i, header = F)
  names(joins[[substr(i, 5, 8)]]) <- c("count", "E030", substr(i, 5, 8))
}
#head(joins[["E031"]])

# loop ofer imported joins and convert to joint binned spectra (matrices)
jMats <- list()
for(i in names(joins)){
  jMats[[i]] <- matrix(0, nrow=max(joins[[i]][,2]) + 1, ncol=max(joins[[i]][,3]) + 1)

  for(j in 1:nrow(joins[[i]])){
    jMats[[i]][joins[[i]][j, 2] + 1, joins[[i]][j, 3] + 1] <- joins[[i]][j, 1]
  }
}

# a function to nicely plot a joint spectrum
plotJMat <- function(jMat, ...){
  image(0:(nrow(jMat)-1), 0:(ncol(jMat)-1), jMat, asp=1, ...)
}
#jMats[["E031"]]
#plotJMat(log10(jMats[["E031"]]))

# a function to generate contours to be added to a plot
contJMat <- function(jMat, ...){
  contour(0:(nrow(jMat)-1), 0:(ncol(jMat)-1), jMat, asp=1, ...)
}

# A function to plot a joint spectrum together with contour lines, log axis ticks,
# and other things
plotMatCont <- function(jMat,...){
  plotJMat(log(jMat),
           col = "#FFFFFF",
           xaxt='n',
           yaxt='n',
           xlab="Copy number bins An1",
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

# a function to compute GS in Mbp from a joint spectrum
# returns two values because each joint spectrum contains data of two individuals
gsFromBinJoin <- function(bMat){
  cs = ncol(bMat)
  rs = nrow(bMat)
  return(c(cols = sum(colSums(bMat)[2:cs] * bin2lin(1:(cs-1))) / 1000000,
           rows = sum(rowSums(bMat)[2:rs] * bin2lin(1:(rs-1))) / 1000000))
}
# GS estimates from joined (binned) spectra
gsE030b <- gsFromBinJoin(jMats[["E031"]])[2]
gsE031b <- gsFromBinJoin(jMats[["E031"]])[1]
gsE032b <- gsFromBinJoin(jMats[["E032"]])[1]
gsE040b <- gsFromBinJoin(jMats[["E040"]])[1]
gsE065b <- gsFromBinJoin(jMats[["E065"]])[1]
gsE068b <- gsFromBinJoin(jMats[["E068"]])[1]
gsE073b <- gsFromBinJoin(jMats[["E073"]])[1]

# The GS estimtes from un-binned and binned spectra match quite well
plot(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a),
     c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b),
     xlab="Un-binned conventional spectra",
     ylab="Computed fro binned joint spectra"
     )

text(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a),
     c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b) + 10,
     labels = c("E030", "E031", "E032", "E040", "sE065", "E068", "E073"))
# Neither slope nor intercept are significant:
summary(lm(c(gsE030a, gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a)~
           c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b),
        offset=c(gsE030b, gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b)
))


# Difference to ref individual (E030, An1)
gsDiffs <- gsE030a - c(gsE031a, gsE032a, gsE040a, gsE065a, gsE068a, gsE073a)
gsDiffsB <- gsE030b - c(gsE031b, gsE032b, gsE040b, gsE065b, gsE068b, gsE073b)
# comparing binned spectra starting from bin joins ####

#  function to produce a k-mer difference data frame from a matrix
# This is going a step back and discarding the information in a joint spectrum.
# The df cols are the sums for each bin over the matrix rows and columns. There
# is also a difference column. This is used to generate the difference graphs.
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



#pdf("cumsumE031.pdf", width=6, height = 4)
#pdf("cumsumE031norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E031"]], main = "Vi - An1")
#dev.off()

compPlot(compDfs[["E031"]], main = "Vi - An1")

#pdf("cumsumE032.pdf", width=6, height = 4)
#pdf("cumsumE032norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E032"]], main = "Ri1 - An1")
#dev.off()
compPlot(compDfs[["E032"]], main = "Ri1 - An1")

#pdf("cumsumE040.pdf", width=6, height = 4)
#pdf("cumsumE040norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E040"]], main = "Ro - An1")
#dev.off()
compPlot(compDfs[["E040"]], main = "Ro - An1")

#pdf("cumsumE065.pdf", width=6, height = 4)
#pdf("cumsumE065norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E065"]], main = "An2 - An1")
#dev.off()
compPlot(compDfs[["E065"]], main = "An2 - An1")

#pdf("cumsumE068.pdf", width=6, height = 4)
#pdf("cumsumE068norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E068"]], main = "Ri2 - An1")
#dev.off()
compPlot(compDfs[["E068"]], main = "Ri2 - An1")

#pdf("cumsumE073.pdf", width=6, height = 4)
#pdf("cumsumE073norm.pdf", width=6, height = 4)
compCumsumPlot(compDfs[["E073"]], main = "Ri3 - An1")
#dev.off()
compPlot(compDfs[["E073"]], main = "Ri3 - An1")



gsDiffs # from full spectra
# Function to compute GS difference from difference data frame:
cumsumDiff <- function(jm){
  -sum(jm[,4] * bin2lin(jm[,1]))/1000000
}

gsDiffsCumsum <- sapply(c("E031", "E032", "E040", "E065", "E068", "E073"), function(x) cumsumDiff(compDfs[[x]]))
# bith are very similar
plot(gsDiffs, gsDiffsCumsum)
abline(0, 1)


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
        ylab = "Contibution to genome size (Mbp)",
        xlab="Sample")
nms <- c("0", "10", "100", "1k", "10k", "100k", "1M", "10M")
legend("topleft",
       fill=grey.colors(7),
       legend=sapply(1:7, function(x) paste0(nms[x],"-", nms[x+1])),
       title="Copy number range"
)
#dev.off()

# differences
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
        ylab = "Contribution to genome size difference (Mbp)",
        xlab="Sample")
nms <- c("0", "10", "100", "1k", "10k", "100k", "1M", "10M")
legend("topleft",
       fill=grey.colors(7),
       legend=sapply(1:7, function(x) paste0(nms[x],"-", nms[x+1])),
       title="Copy number range"
       )
#dev.off()

# partitioning gs variance ####
gsDat <- data.frame(spec=c("An","An","Vi","Ro","Ri","Ri","Ri"),
                    genSiz=c(999.98, 989.23, 1055.93, 1227.92, 1126.64, 1096.44, 1104.84)
)
summary(lm(genSiz ~ spec, gsDat))
var(gsDat[,2])
anova(lm(genSiz ~ spec, gsDat))
a <- anova(lm(genSiz ~ spec, gsDat))
str(a)
a$`Pr(>F)`


# This p-value is based on very few samples. How often do we get low p-values
# is we permute the data?
set.seed(123345)
resample <- function(dff, n=999){
  sapply(1:n, function(x){
    a <- dff
  a[,2] <- dff[sample(1:7),2]
  anova(lm(genSiz ~ spec, a))$`Pr(>F)`[1]
  })
}
pvals <- resample(gsDat)
hist(pvals)
quantile(pvals, 0.05)
# Rarely.

# synthetic data ####
# For schematic figure 1

# generated  using Tetmer code
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

xticks = c(1, 2, 4, 10, 100, 1000, 10000, 100000, 1000000)
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

