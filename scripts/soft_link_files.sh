for file in /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/202305*/*.fastq.gz; 
do 
    #name=$(echo $file | sed -n "s/.*L29P4_1_\(.*\)_[ATGC]\+_[ATGC]\+.*\(R[12]\).*/\1_\2.fastq.gz/p")
    name=$(echo $file | sed -n "s/.*L29P4_1_\(.*\)_S.*\(R[12]\).*/\1_\2.fastq.gz/p")
    ln -s $file $name
done

