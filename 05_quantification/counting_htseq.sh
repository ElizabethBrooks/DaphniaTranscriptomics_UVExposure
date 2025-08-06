#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput

# script to perform htseq-count counting of trimmed, aligned, then name sorted
# paired end reads
# Usage: qsub counting_htseq.sh sortedFolder
# Usage Ex: qsub counting_htseq.sh sorted_name

#Required modules for ND CRC servers
module load bio/3.0
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2

# retrieve input folder name
sortedFolder=$1

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../"inputData/shortReads/inputPaths_D_melanica.txt" | tr -d " " | sed "s/genomeFeatures://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

# setup the inputs path
inputsPath=$outputsPath"/"$sortedFolder

#Move to outputs directory
cd "$outputsPath"

# create outputs directory
outputFolder="counted"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsPath"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using htseq-count
	echo "Sample $curSampleNoPath is being counted..."
	#Determine which flags to use based on sorting method
	if [[ "$1"  == *name* ]]; then
		#Use name sorted flag (default)
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
	elif [[ "$1"  == *coordinate* ]]; then
		#Use coordinate sorted flag
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
	else
		echo "ERROR: The bam file "$f1" was not found... exiting"
		exit 1
	fi
	echo "Sample $curSampleNoPath has been counted!"
done
#Copy previous summaries
cp "$inputsDir"/*.txt "$outputFolder"

# clean up
rm -r $inputsDir
