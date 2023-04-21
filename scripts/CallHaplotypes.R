
library("HaplotypR")
library("ShortRead")

args <- commandArgs(trailingOnly=TRUE)

marker = args[1]
merged_file = args[2]
markerTab <- read.csv2(args[3])
outputDir <- args[4]

refSeq = DNAStringSet(markerTab$ReferenceSequence[markerTab$MarkerID == marker])
names(refSeq) <- marker
dir.create(outputDir)

if (file.info(merged_file)$size == 0) {
    quit()
}
################################################################################
#Calculate mismatch rate and call SNPs
################################################################################
# Options
minMMrate <- 0.05
minOccGen <- 1

# Calculate mismatch rate
seqErrLst <- calculateMismatchFrequencies(merged_file, 
                                        refSeq, 
                                        method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                        minCoverage=1L)

seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s.txt", marker)), sep="\t", row.names=F)

# Call SNPs
potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
snpRef <- unlist(lapply(potSNP, function(snp){
as.character(subseq(refSeq[marker], start=snp, width=1))
}))
if (length(potSNP) == 0){
    snpRef <- as.character(subseq(refSeq[marker], start=1, width=1))
    snps <- data.frame(Chr=marker, Pos=1, Ref=snpRef, Alt=snpRef, stringsAsFactors=F)
}else{
snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
}
write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s.txt", 
                                                    minMMrate*100, minOccGen, marker)), 
            row.names=F, col.names=T, sep="\t", quote=F)


################################################################################
# call haplotypes
################################################################################
# call haplotype options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2
minCovSample <- 25

sampleTab = data.frame(ReadFile = merged_file, MarkerID = marker, SampleID = 0, SampleName = '')

snpList = list(snps)
names(snpList) <- marker


# call final haplotypes
finalTab <- createFinalHaplotypTable(
  outputDir = outputDir, sampleTable = sampleTab, markerTable = sampleTab, referenceSeq = refSeq,
  snpList = snpList, postfix = '', minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)

