# Write out a 2D histogram of a transformed k-mer dump

# Character vector of command line arguments
args = commandArgs(trailingOnly=TRUE)

# Read in transformed dump:
td = read.table(args[1], header=F, sep=" ")

# Decide on dimensions for a count matrix
xdim = max(td[,2]) + 1 # number of matrix rows
ydim = max(td[,3]) + 1 # number of matrix cols
xbins = 1:xdim
ybins = 1:ydim

# Initialise matrix
countMat = matrix(0, nrow=xdim, ncol = ydim)

# Fill matrix by going through dump line by line
for(i in 1:nrow(td)){
    countMat[td[i, 2], td[i, 3]] <- countMat[td[i, 2], td[i, 3]] + td[i, 1]
}
png("hist2d.png")
image(xbins, ybins, countMat, asp=1, main="True counts", xlab="Bins, Sample 1", ylab = "Bins, Sample 2")
abline(0, 1)
dev.off()

png("hist2dlog.png")
image(xbins, ybins, log10(countMat), asp=1, main="Log counts", xlab="Bins, Sample 1", ylab = "Bins, Sample 2")
abline(0, 1)
dev.off()
