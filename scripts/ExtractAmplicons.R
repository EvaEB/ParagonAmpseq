######################################
##
## Rscript ExtractAmplicons.R [inputfile] [outputfile] [marker] [primerfile] [primerVariation]
##
## extract amplicons and remove primers from fused reads 
## required arguments:
## * inputfile - input .fastq file with reads aligning to the amplicon region in question
## * outputfile - output .fastq, will contain filtered & trimmed reads
## * marker - name of the marker, as in [primerfile]
## * primerfile - tab-separated file containing the primer sequences with (at least) columns "MarkerID", "Forward", "Reverse"
## * primervariation  - allowed fraction of variation in the primer sequence
##
## Author: Dr. Eva Bons
##
## Date Created: 2023-02-21
##
## Copyright (c) Eva Bons, 2023
## Email: eva.bons@gmail.com
##
######################################


library("ShortRead")
library("Biostrings")
library("tools")

args <- commandArgs(trailingOnly=TRUE)
inputfile = args[1]
outputfile = args[2]

marker = args[3]
primerfile = args[4]
primer_variation = as.numeric(args[5])

######################################
# INPUT CHECKS
#####################################
## inputfile must be fastq, and must exist
if (file_ext(inputfile) != "fastq") {
  stop("input file must be .fastq")
}

if (!file.exists(inputfile)) {
  stop(paste("input file '",inputfile, "' does not exist",sep=""))
}


## primerfile must exist
if (!file.exists(primerfile)) {
  stop(paste("primer file '",primerfile, "' does not exist",sep=""))
}
markerTab <- read.csv2(primerfile)

## marker must be present in primerfile
if (!any(markerTab$MarkerID == marker)){
  stop(paste("marker",marker,"not found in primerfile"))
}
  
## primer variation must be float, <1, >=0
if (primer_variation >= 1 | primer_variation < 0){
  stop(paste("primer variation (", primer_variation, ") must be a float, <1, >=0",sep=""))
}

## if not provided, it defaults to 0.2
##TODO
##primer_variation = 0.2
  
###########################################
# read in files & extract relevant info
###########################################
# Read fastq input file
fastq <- readFastq(inputfile,withIds = TRUE)

# Define primer strings
fw_primer <- DNAString(markerTab$Forward[markerTab$MarkerID == marker])
rv_primer <- DNAString(markerTab$Reverse[markerTab$MarkerID == marker])

############################################
# locate primers & trim
###########################################
# fw primer
maxmismatchfw = round(length(fw_primer)*primer_variation) 
fw_matches = vmatchPattern(fw_primer, sread(fastq),max.mismatch=maxmismatchfw)
start_trim = end(fw_matches) + 1

# rv primer
maxmismatchrv = round(length(rv_primer)*primer_variation)
rv_matches = vmatchPattern(rv_primer, sread(fastq),max.mismatch=maxmismatchrv)
end_trim = start(rv_matches) - 1

#for amplicons on the reverse strand. This is a little hacky but it works
if (all(all(start_trim > end_trim))){ #it found the forward primer in the end and vice versa
  start_trim = end(rv_matches) + 1
  end_trim = start(fw_matches) - 1
}



# find unmatched reads to filter out
filter = (which(sapply(start_trim, isEmpty) | sapply(end_trim, isEmpty)))
write(sprintf("sequences with missing primers: %d (%.2f %%)", length(filter), 100*length(filter)/ length(fastq)),stdout())

# filter out these reads
start_trim = start_trim[-filter]
end_trim = end_trim[-filter]
fastq_filtered = fastq[-filter]

start_trim = unlist(lapply(start_trim, as.list), recursive=TRUE)
end_trim = unlist(lapply(end_trim, as.list), recursive=TRUE)

# find invalid matches (start > end)
filter_valid = start_trim < end_trim
number_invalid = length(start_trim)-sum(filter_valid)
write(sprintf("sequences with invalid primer positions: %d (%.2f %%)",number_invalid, 100*number_invalid/ length(fastq)),stdout())

start_trim = start_trim[filter_valid]
end_trim = end_trim[filter_valid]
fastq_filtered = fastq_filtered[filter_valid]

# trim reads
trimmed = lapply(1:length(fastq_filtered), function(i){ subseq(sread(fastq_filtered[i])[[1]],start_trim[[i]], end_trim[[i]])})
trimmed = DNAStringSet(trimmed)

#trim qualities
quals = quality(fastq)[-filter][filter_valid]
trimmed_quals = lapply(1:length(fastq_filtered), function(i){ subseq(quals[[i]],start_trim[[i]], end_trim[[i]])})

#save result
names_fastq = id(fastq)[-filter][filter_valid]
names(trimmed) = names_fastq

writeXStringSet(x=trimmed, qualities=BStringSet(trimmed_quals) , outputfile, format="fastq")


