#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_test_jobOutput
#$ -pe smp 4

# Script to perform hisat2 alignment of trimmed
# paired end reads
# Note that a hisat2 genome refernce build folder needs to be generated first
# usage: qsub realignment_hisat2_test.sh
# un-conc and al-conc
# test
## job 2073938
# EGAPx test
## job 

#Required modules for ND CRC servers
module load bio/2.0
#module load bio/hisat2/2.1.0

#Retrieve genome reference absolute path for alignment
#buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
buildFile=$(grep "genomeReference:" ../inputData/shortReads/inputPaths_EGAPx_D_pulex.txt | tr -d " " | sed "s/genomeReference://g")
# Retrieve analysis outputs absolute path
#outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_D_pulex.txt" | tr -d " " | sed "s/outputs://g")
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_EGAPx_D_pulex.txt" | tr -d " " | sed "s/outputs://g")
# Retrieve paired reads absolute path for alignment
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_pulex.txt" | tr -d " " | sed "s/pairedReads://g")

# Make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

# set inputs absolute path
inputsFolder=$outputsPath"/aligned_conc"

# move to outputs directory
cd "$outputsPath"

# set output directory name
outputFolder=$outputsPath"/realigned_conc"
# create output directory
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
#Add software version to output summary file
hisat2 --version > $inputOutFile

#Build output directory for Hisat reference
buildOut="$outputsPath"/"reference_hisat2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//' | sed 's/\.fna//' | sed 's/\.fa//')

#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in $inputsFolder"/"*/; do
	# status message
	echo "Processing directory $f1 ..."
	curSample=$f1
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Run hisat2 with default settings
	echo "Sample $curSampleNoPath is being aligned and converted..."
	hisat2 -p 4 -q -x $buildOut"/"$buildFileNoEx -1 $curSample"/un_conc.fq.1.gz" -2 $curSample"/un_conc.fq.2.gz" -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam \
	--un-conc-gz "$outputFolder"/"$curSampleNoPath"/un_conc.fq.gz --al-conc-gz "$outputFolder"/"$curSampleNoPath"/al_conc.fq.gz --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
	#Remove the now converted .sam file
	rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	# status message
	echo "Sample $curSampleNoPath has been aligned and converted!"
done
