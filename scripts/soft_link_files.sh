for file in ../../../../../../../projects/openbis/userstore/uni_basel_stph_nsanzabana/2022*/*.fastq.gz; 
do 
    name=$(echo $file | sed -n "s/.*KMM6N_1_\(.*\)_[ATGC]\+_[ATGC]\+.*\(R[12]\).*/\1_\2.fastq.gz/p")
    ln -s $file $name
done

