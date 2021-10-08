setwd("~/Dropbox/manuscripts/2108_kmers_gs/data/RE1813/")


# read file of super clusters ####
# this table was manually truncated, check example file

sc <- read.table("SUPERCLUSTER_TABLE_trunc.csv", header=T, stringsAsFactors = F)
str(sc)
# all lines in order? yes!
all(sc[,1] == 1:50)

# Are all clusters uniquely assigned to one super cluster?
length(do.call(c, sapply(1:50, function(ln){as.numeric(strsplit(sc[ln,2], ",")[[1]])})))
length(unique(do.call(c, sapply(1:50, function(ln){as.numeric(strsplit(sc[ln,2], ",")[[1]])}))))
# Yes!

##############################
# ADUST these for your needs:
smallestCluster <- 201 # Depends on the specific RE data set! Pick bottom cluster from cluster summary HTML.
outPath = "~/temp/2108_gs_kmers/"
rePref = "/mnt/1TBSSD/RE_Euphrasia_GS_ms/seqclust/clustering/clusters/" # the directory that contains the RE cluster directories

# Function to paste all read files of a super cluster's member clusters
pasteReadFiles <- function(ln){

  print(paste0("SU: ", ln))
  membs <- as.numeric(strsplit(sc[ln,2], ",")[[1]])
  if (file.exists(paste0(outPath, sprintf("su%02d", ln)))){
    #Delete file if it exists
    file.remove(paste0(outPath, sprintf("su%02d", ln)))
  }
  for(i in membs){
    print(i)
    if( i <= smallestCluster){
      cat(readLines(paste0(rePref, sprintf("dir_CL%04d/reads.fasta", i))),
          file=paste0(outPath, sprintf("su%02d", ln)),
          sep="\n",
          append=T)
    }
  }

}
for(su in 1:50) pasteReadFiles(su)

# then, in a terminal, run: 
# touch allSuReadsConcat.fa; for i in su??; do echo ">$i" >> allSuReadsConcat.fa; grep -v ">" $i >> allSuReadsConcat.fa; done
# and then run UniqueKMER (https://github.com/OpenGene/UniqueKMER/archive/refs/tags/v0.1.0.tar.gz) on the concatenated file
# Finally, generate KMC data bases for each of the ca. 50 sets of unique k-mers. There might be fewer than 50 if some super clusters did not contain unique k-mers.