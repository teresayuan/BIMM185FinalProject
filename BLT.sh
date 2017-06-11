#!/bin/bash

# 1 --> reference
# 2 --> 

DIR=/home/bimm185

if [ "$1" == "-h" ]
then
	if [ "$1" != "-h" ] 
	then
		echo "wrong number of files"
	fi
	echo "scriptname reference sample1forward sample1reverse sample2forward sample2reverse"

	exit 1
fi

echo "hello world"

mkdir $DIR/working_test
mkdir $DIR/working_test/fastqc

# for f in $(ls $DIR/data/*.fastq); do echo $f; fastqc -o $DIR/working_test/fastqc $f ; done

for f in $(ls $DIR/data/*.fastq); do echo $f; python modfastq.py $f ${f/data/working_test}; done



#bwa index -p $f{${f/reference/working_test}%%.fna}; done
cp $DIR/reference/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/.
for f in $(ls $DIR/working_test/*.fna); do echo $f; bwa index $f; done

bwa mem $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC8_1.fastq $DIR/working_test/UC_2.fastq > $DIR/working_test/UC8.sam

samtools view -o $DIR/working_test/UC8.bam -b $DIR/working_test/UC8.sam

samtools flagstat $DIR/working_test/UC8.bam > $DIR/working_test/UC8_align_stats.txt

samtools sort $DIR/working_test/UC8.bam -o $DIR/working_test/UC8_sorted.bam 

samtools mpileup -f $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC8_sorted.bam > $DIR/working_test/UC8.mpileup

java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp $DIR/working_test/UC8.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > $DIR/working_test/UC8.vcf


