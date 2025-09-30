#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N BHB_blastp_jobOutput
#$ -pe smp 8

# script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
# usage: qsub BHB_blastp.sh
## genomeReference
## job 2083470
## proteins
## job

# load necessary modules for ND CRC servers
module load bio/2.0

# retrieve genome reference absolute path for alignment
#inputDB=$(grep "proteins:" ../inputData/shortReads/inputPaths_ZQ_D_melanica.txt | tr -d " " | sed "s/proteins://g")
inputDB=$(grep "proteins:" ../inputData/shortReads/inputPaths_EGAPx_D_melanica.txt | tr -d " " | sed "s/proteins://g")
# retrieve genome reference absolute path for alignment
#inputQuery=$(grep "proteins:" ../inputData/shortReads/inputPaths_D_pulex.txt | tr -d " " | sed "s/proteins://g")
inputQuery=$(grep "proteins:" ../inputData/shortReads/inputPaths_EGAPx_D_pulex.txt | tr -d " " | sed "s/proteins://g")

# retrieve analysis outputs absolute path
#outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
outputsPath=$(grep "outputs:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/outputs://g")
# retrieve paired reads absolute path for alignment
#readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_ZQ_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")
readPath=$(grep "pairedReads:" ../"inputData/shortReads/inputPaths_EGAPx_D_melanica.txt" | tr -d " " | sed "s/pairedReads://g")

# make a new directory for project analysis
projectDir=$(basename $readPath)
outputsPath=$outputsPath"/"$projectDir

# move to outputs directory
cd "$outputsPath"

# set output directory name
outputFolder=$outputsPath"/RBHB_proteins"
# create output directory
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists..."
fi

# name output file of inputs
inputOutFile=$outputFolder"/software_summary.txt"
# add software version to output summary file
blastp -version > $inputOutFile

# check if first DB of seqs exsists
directoryDB=$(dirname "$inputDB")
if [ -f "$inputDB".phr ]; then
	# status message
	echo "Using exsisting "$inputDB".phr DB..."
else # make blastable DB of transcriptome
	# status message
	echo "Creating "$inputDB".phr DB..."
	# determine current script location
	currLoc=$(echo $PWD)
	# move to DB directory
	cd $directoryDB
	makeblastdb -in $inputDB -dbtype prot
	# move back to script location
	cd $currLoc
fi

# check if second DB of seqs exsists
directoryQuery=$(dirname "$inputQuery")
if [ -f "$inputQuery".phr ]; then
	# status message
	echo "Using exsisting "$inputQuery".phr DB..."
else # make blastable DB of seqs
	# status message
	echo "Creating "$inputQuery".phr DB..."
	# move to DB directory
	cd $directoryQuery
	makeblastdb -in $inputQuery -dbtype prot
fi

# use blastp to search a database
# output with outfmt6 header:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# switch query and search paths for reciprocal search
echo "Beginning blastp search..."
blastp -num_threads 8 -query "$inputQuery" -db "$inputDB" -outfmt 6 -evalue 0.01 -num_alignments 1 > "$outputFolder"/RBHB_pulex_melanica.outfmt6
echo "Finished reciprocal blastp database search!"

# status message
echo "Analysis complete!"
