#!/bin/bash

# Pipeline sin threshold

python3 FileCreator.py ./Monocyte1 Monocyte1_States_0 0.0
python3 FileCreator.py ./Monocyte2 Monocyte2_States_0 0.0

mv ./Monocyte1/Monocyte1_States_0_E1.bed -t ./
mv ./Monocyte1/Monocyte1_States_0_E2.bed -t ./
mv ./Monocyte2/Monocyte2_States_0_E1.bed -t ./
mv ./Monocyte2/Monocyte2_States_0_E2.bed -t ./

echo "Creando los archivos de intersección para los estados E1 y E2"
echo " "

bedtools intersect -a ./Monocyte1_States_0_E1.bed -b ./Monocyte2_States_0_E1.bed > Monocytes_0_E1.bed
bedtools intersect -a ./Monocyte1_States_0_E2.bed -b ./Monocyte2_States_0_E2.bed > Monocytes_0_E2.bed

echo "Creando el directorio de resultados"
echo " "

mkdir Results_0

echo "Moviendo todos los ficheros creados al directorio de resultados"
echo " "

mv $(ls | grep -w Monocyte._States_0_.*) -t ./Results_0
mv $(ls | grep -w Monocytes_0_E.*) -t ./Results_0

# Pipeline con threshold de 0.95

echo "#######################################################################################################################################"
echo " "

python3 FileCreator.py ./Monocyte1 Monocyte1_States_0.95 0.95
python3 FileCreator.py ./Monocyte2 Monocyte2_States_0.95 0.95

mv ./Monocyte1/Monocyte1_States_0.95_E1.bed -t ./
mv ./Monocyte1/Monocyte1_States_0.95_E2.bed -t ./
mv ./Monocyte2/Monocyte2_States_0.95_E1.bed -t ./
mv ./Monocyte2/Monocyte2_States_0.95_E2.bed -t ./

echo "Creando los archivos de intersección para los estados E1 y E2"
echo " "

bedtools intersect -a ./Monocyte1_States_0.95_E1.bed -b ./Monocyte2_States_0.95_E1.bed > Monocytes_0.95_E1.bed
bedtools intersect -a ./Monocyte1_States_0.95_E2.bed -b ./Monocyte2_States_0.95_E2.bed > Monocytes_0.95_E2.bed

echo "Creando el directorio de resultados"
echo " "

mkdir Results_0.95

echo "Moviendo todos los ficheros creados al directorio de resultados"
echo " "

mv $(ls | grep -w Monocyte._States_0.95_.*) -t ./Results_0.95
mv $(ls | grep -w Monocytes_0.95_E.*) -t ./Results_0.95
