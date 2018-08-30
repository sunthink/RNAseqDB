library(Rsubread)
library(limma)
library(edgeR)


args <- commandArgs( trailingOnly = TRUE )
bamFile <- args[1]
gtfFile <- args[2]
nthreads <- args[3]
outFilePref <- args[4]

outStatsFilePath  <- paste(outFilePref, '.stat',  sep = '');
outCountsFilePath <- paste(outFilePref, '.count', sep = '');
outFpkmFilePath   <- paste(outFilePref, '.fpkm',  sep = '');
outTpmFilePath    <- paste(outFilePref, '.tpm',   sep = '');

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, nthreads=nthreads, isPairedEnd=TRUE)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = t(t(fpkm)*1e6/(colSums(fpkm)))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(cbind(fCountsList$annotation[,1], fpkm), outFpkmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#write.table(cbind(fCountsList$annotation[,1], log2(fpkm + 1)), outFpkmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], tpm), outTpmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#write.table(cbind(fCountsList$annotation[,1], log2(tpm + 1)), outTpmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
