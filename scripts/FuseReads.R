#!/usr/bin/env Rscript
# usage Rscript --vanilla FuseReads.R [outputDir] [FileFW] [FileRV]

library("HaplotypR")
library("ShortRead")

args <- commandArgs(trailingOnly=TRUE)

################################################################################
##Fuse paired reads
################################################################################
outProcFiles <- args[1]
FileR1  <- args[2]
FileR2  <- args[3]

dir.create(outProcFiles)

procReadsMerge <- mergeAmpliconReads(as.character(FileR1), as.character(FileR2), outProcFiles, method='vsearch',mergePrefix="", trimFilenameExt="_R1\\.fastq.gz$")
