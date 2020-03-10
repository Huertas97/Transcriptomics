#!/bin/bash

file_in=$1
file_out=$2

bedtools getfasta -fi hg19.fa -bed $file_in -fo $file_out 

# ./getFastaFromBed_MEME.sh bed_E1_filtrado_0.7/intersect_E1_f0.7.bed bed_E1_filtrado_0.7/intersect_E1_f0.7.fasta