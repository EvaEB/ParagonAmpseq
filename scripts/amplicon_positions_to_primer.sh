echo "MarkerID;Forward;Reverse;ReferenceSequence"

{
    read
    while read -r line; do 
        echo $line | awk '{printf $1";"}'
        fw=$(echo $line | awk '{printf $2":"$3"-"$4}')
        rv=$(echo $line | awk '{printf $2":"$5"-"$6}')
        ref=$(echo $line | awk '{printf $2":"$4+1"-"$5-1}')

        samtools faidx -n 1000 $2 $fw | tail -1 | tr '\n' ';'
        samtools faidx -n 1000 $2 $rv | tail -1 | tr '\n' ';'
        samtools faidx -n 1000 $2 $ref | tail -1

    done 
} < $1

Amplicon;>chrom:left_start-left_end;>chrom:right_start-right_end;>chrom:1--1
