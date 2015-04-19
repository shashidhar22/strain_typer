#! /bin/bash
FILES=/data/assembly/*-SPADES-contig.fa
for file in $FILES
do 
    echo "k-merizing $file"
    fasta=${file##*/}
    ./count -f fasta -d 4 -l 4 -k 8 -o kmer_out/${fasta%.*}".kc" $file
done
