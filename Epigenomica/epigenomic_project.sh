#!/bin/bash
mkdir -p "resultados" # no error if existing, make parent directories as needed

# -wb to report the original “B” feature when an overlap is found
bedtools intersect -wb -a Monocyte1_11_Master_11_segments.bed -b Monocyte2_11_Master_11_segments.bed -sorted > ./resultados/temp_E1.aux
# awk will select only overlapped segments belonging to E1 in both files
# +1 para que este 1 based
awk '{OFS="\t"; if ($4 == $8 && $4 == "E1") print $1, $2+1, $3+1}' ./resultados/temp_E1.aux > ./resultados/E1.bed

# The same for E2
bedtools intersect -wb -a Monocyte1_11_Master_11_segments.bed -b Monocyte2_11_Master_11_segments.bed -sorted > ./resultados/temp_E2.aux
# awk will select only overlapped segments belonging to E1 in both files
# +1 para que este 1 based
awk '{OFS="\t"; if ($4 == $8 && $4 == "E2") print $1, $2+1, $3+1}' ./resultados/temp_E2.aux > ./resultados/E2.bed

rm -r ./resultados/*.aux 

# Añadimos "chr" a la columna correspondiente al cromosoma en GFF
# awk '{if !($1 == "#") print "chr"$4}' GRCh37_latest_genomic.gff | head