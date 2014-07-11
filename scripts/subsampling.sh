#!/bin/bash
source /home/arun/.bash_profile
#sample reads, place in /temp/[]x/
#map reads to refseq, end up with a sorted bam in results
#use an Rscript to plot the coverage

set -e
set -u
suffix='.subsamp'
reads='../data/trimmed-reads/3163_Y920_NoIndex_L002_R1_001trimmed.fasta'
refseq='../data/refseq/citrus-yellow-vein-associated-virus.fasta'
covs=('5x' '10x' '20x' '30x' '40x' '50x')
readarray=(260 520 1040 1560 2080 2600)

sample_reads() {
	bioawk -c fastx -v k=$2 '{y=x++<k?x-1:int(rand()*x);if(y<k)a[y]=">"$name"\n"$seq}END{for(z in a)print a[z]}' $reads > ../temp/$1/sampled_reads.fasta
}

map () {
	bwa index -a bwtsw $refseq
	bwa bwasw $refseq ../temp/$1/sampled_reads.fasta > ../temp/$1/alignment.sam
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
	sample_reads ${covs[$index]} ${readarray[$index]}
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
