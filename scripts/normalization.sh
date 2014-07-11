#!/bin/bash
source /home/arun/.bash_profile
#this script uses bbnorm for normalization rather than subsampling
#sample reads, place in /temp/[]x/
#map reads to refseq, end up with a sorted bam in results
#use an Rscript to plot the coverage

set -e
set -u
suffix='.norm'
reads='../data/trimmed-reads/3163_Y920_NoIndex_L002_R1_001trimmed.fasta'
refseq='../data/refseq/citrus-yellow-vein-associated-virus.fasta'
covs=('5x' '10x' '20x' '30x' '40x' '50x')
karray=(2 4 8 12 16 20)

sample_reads() {
	bash bbnorm.sh in=$reads out=../temp/$1/normalized_reads.fasta target=$2 -Xmx2g
}

map () {
	bwa index -a bwtsw $refseq
	bwa bwasw $refseq ../temp/$1/normalized_reads.fasta > ../temp/$1/alignment.sam
	samtools view -bS -F 4 ../temp/$1/alignment.sam > ../temp/$1/mapped.alignment.bam
	samtools sort ../temp/$1/mapped.alignment.bam ../results/$1/mapped.sorted.alignment
	samtools index ../results/$1/mapped.sorted.alignment.bam
}

plot() {
	genomeCoverageBed -ibam ../results/$1/mapped.sorted.alignment.bam -g $refseq -d > ../results/$1/coverageHist$suffix.txt
	cat ../results/$1/coverageHist$suffix.txt | ./rplot.r ../results/$1/plot$suffix.png
}

for index in ${!covs[*]}
do
	printf "%4d\n" $index
	echo "###################"
	echo "sampling reads..."
	echo "###################"
	sample_reads ${covs[$index]} ${karray[$index]}
	echo "###################"
	echo "mapping..."
	echo "###################"
	map ${covs[$index]}
	echo "###################"
	echo "plotting..."
	echo "###################"
	plot ${covs[$index]}
done
echo "DONE"
