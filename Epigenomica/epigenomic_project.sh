#!/bin/bash
mkdir -p "resultados" # no error if existing, make parent directories as needed

# -wb to report the original “B” feature when an overlap is found
bedtools intersect -wb -a Monocyte1_11_Master_11_segments.bed -b Monocyte2_11_Master_11_segments.bed -sorted > ./resultados/temp_E1.aux
# awk will select only overlapped segments belonging to E1 in both files
awk '{OFS="\t"; if ($4 == $8 && $4 == "E1" || $4 == $8 && $4 == "E2") print $1, $2, $3, $4}' ./resultados/temp_E1.aux > ./resultados/E1_E2.bed

rm -r ./resultados/*.aux 


# comprobacion
# grep -cE "^chr.*E1[^0-9].*E1$" temp_E1.aux 
# grep -cE "^chr.*E2.*E2$" temp_E1.aux