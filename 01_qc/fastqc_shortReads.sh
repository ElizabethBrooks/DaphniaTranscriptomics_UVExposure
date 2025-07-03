#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N fastqc_MultiGenome_jobOutput

# Script to perform fastqc quality control of paired end reads
# Usage: qsub fastqc_shortReads.sh inputsType
# Usage Ex: qsub fastqc_shortReads.sh raw
## job 
# Usage Ex: qsub fastqc_shortReads.sh trimmed
## job 

# Required modules for ND CRC servers
module load bio

# retrieve input arguments
inputsType=$1

# Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
# Retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_D_melanica.txt" | tr -d " " | sed "s/outputs://g")

# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir
#mkdir $outputsPath

# check inputs type
if [[ $inputsType == "trimmed" ]]; then
	# set the directory for inputs
	readPath=$outputsPath"/trimmed"
	# set the directory for analysis
	qcOut=$outputsPath"/qc_trimmed"
elif [[ $inputsType == "raw" ]]; then
	# set the directory for analysis
	qcOut=$outputsPath"/qc_raw"
fi

# Make a new directory for analysis
mkdir $qcOut
# Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $qcOut directory already exsists... please remove before proceeding."
	exit 1
fi
# Move to the new directory
cd $qcOut

# Name version output file
versionFile=$qcOut"/software_version_summary.txt"

# Report software version
fastqc -version > $versionFile

# perform QC
fastqc $readPath"/"*\.f*q.gz -o $qcOut

# run multiqc
multiqc $qcOut -o $qcOut -n "multiqc_raw"

# Print status message
echo "Analysis complete!"
