python3 FileCreator.py

bedtools intersect -a ./Monocyte1/States_E1.bed -b ./Monocyte2/States_E1.bed > Monocytes_E1.bed
bedtools intersect -a ./Monocyte1/States_E2.bed -b ./Monocyte2/States_E2.bed > Monocytes_E2.bed



