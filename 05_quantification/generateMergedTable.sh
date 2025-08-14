#!/bin/bash

# script to generate a merge gene counts file
# usage: bash generateMergedTable.sh

# set the inputs path
inputsPath="/scratch365/ebrooks5/D_melanica_UV_exposure/short_read_data_processed_EGAPx/Pfrender_MP-3533_250512_CMG/counted"

# initialize the merged counts file
echo "gene" > $inputsPath"/counts_merged.csv"
cat $inputsPath"/MUV10_S10_L004/counts.txt" | cut -f1 >> $inputsPath"/counts_merged.csv"

# merge counts for each sample
for i in $inputsPath"/"*"/"; do 
	newName=$(basename $i | sed "s/_S.*_L004//g")
	echo $newName > $inputsPath"/counts_merged.tmp.csv"
	cat $i/counts.txt | cut -f2 >> $inputsPath"/counts_merged.tmp.csv"
	paste -d, $inputsPath"/counts_merged.csv" $inputsPath"/counts_merged.tmp.csv" >> $inputsPath"/counts_merged.tmp.csv"
done

# clean up
rm $inputsPath"/counts_merged.tmp.csv"
