# script to extract primer and reference sequences per amplicon from a reference genome given a file with amplicon positions
# 
# usage:
# bash amplicon_positions_to_primer.sh amplicon_positions.csv reference.fasta
#   amplicon_positions.csv: file with amplicon info: must be formatted as (make sure all primer positions are refering to the + strand!!!):
#     `amplicon_name chromosome  fw_primer_start fw_primer_end   rv_primer_start rv_primer_end`
echo "MarkerID;Forward;Reverse;ReferenceSequence" #header

{
    read
    while read -r line; do 
        echo $line | awk '{printf $1";"}' #marker name
        fw=$(echo $line | awk '{printf $2":"$3"-"$4}') #fw details
        rv=$(echo $line | awk '{printf $2":"$5"-"$6}') #rv details
        ref=$(echo $line | awk '{printf $2":"$4+1"-"$5-1}') #reference details

        samtools faidx -n 1000 $2 $fw | tail -1 | tr '\n' ';' #extract forward sequence
        samtools faidx -n 1000 $2 $rv | tail -1 | tr '\n' ';' #extract reverse sequence
        samtools faidx -n 1000 $2 $ref | tail -1 #extract amplicon reference sequence

    done 
} < $1
