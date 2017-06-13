#!/bin/bash

# 1 --> reference
# 2 --> forward strand of sample 1
# 3 --> reverse strand of sample 1
# 4 --> forward strand of sample 2
# 5 --> reverse strand of sample 2


if [ "$1" == "-h" ] || [ "$#" -ne 5 ];
then
	if [ "$1" != "-h" ];
	then
		echo "wrong number of files"
	fi
	echo "BLT.sh reference sample1forward sample1reverse sample2forward sample2reverse"

	exit 1
fi

echo "Welcome to BLT!"

a=$(basename $1)
working_test=${a%%.*}

mkdir $working_test
# cd $working_test
mkdir $working_test/fast_qc

cp $1 $working_test
cp $2 $working_test
cp $3 $working_test
cp $4 $working_test
cp $5 $working_test

fastqc -o $working_test/fast_qc $2
fastqc -o $working_test/fast_qc $3
fastqc -o $working_test/fast_qc $4
fastqc -o $working_test/fast_qc $5

# for f in $(ls $DIR/$working_test/*.fna); do echo $f; bwa index $f; done

bwa index $1

s1=$(basename $2 .fastq)
s2=$(basename $4 .fastq)

bwa mem $1 $2 $3 > $working_test/$s1.sam
bwa mem $1 $4 $5 > $working_test/$s2.sam

samtools view -o $working_test/$s1.bam -b $working_test/$s1.sam
samtools view -o $working_test/$s2.bam -b $working_test/$s2.sam

samtools flagstat $working_test/$s1.bam > $working_test/$s1_align_stats.txt
samtools flagstat $working_test/$s2.bam > $working_test/$s2_align_stats.txt

samtools sort $working_test/$s1.bam -o $working_test/$s1_sorted.bam 
samtools sort $working_test/$s2.bam -o $working_test/$s2_sorted.bam 

samtools mpileup -f $1 $working_test/$s1_sorted.bam > $working_test/$s1.mpileup
samtools mpileup -f $1 $working_test/$s2_sorted.bam > $working_test/$s2.mpileup

java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp $working_test/$s1.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > $working_test/$s1.vcf
java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp $working_test/$s2.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > $working_test/$s2.vcf


#VCF File comparison
grep '#' $working_test/$s1.vcf > $working_test/VCF_header.txt
header=$(grep -c '#' $working_test/$s1.vcf)
header="$((header+1))"

echo 'number of headers'
echo $header
tail -n +$header $working_test/$s1.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > $working_test/$s1_sorted_vcf
tail -n +$header $working_test/$s2.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > $working_test/$s2_sorted_vcf

mkdir $working_test/vcf_output

join -1 1 -2 1 -o 1.3,1.5,1.6 $working_test/$s1_sorted_vcf $working_test/$s2_sorted_vcf > $working_test/vcf_output/common.vcf

join -1 1 -2 1 -v1 $working_test/$s1_sorted_vcf $working_test/$s2_sorted_vcf > $working_test/vcf_output/only_$s1.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' $working_test/vcf_output/only_$s1.vcf > $working_test/vcf_output/only_$s1_orig.vcf

cat $working_test/VCF_header.txt $working_test/vcf_output/only_$s1_orig.vcf > $working_test/vcf_output/only_$s1.vcf

join -1 1 -2 1 -v2 $working_test/$s1_sorted_vcf $working_test/$s2_sorted_vcf > $working_test/vcf_output/only_$s2.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' $working_test/vcf_output/only_$s2.vcf > $working_test/vcf_output/only_$s2_orig.vcf

cat $working_test/VCF_header.txt $working_test/vcf_output/only_$s2_orig.vcf > $working_test/vcf_output/only_$s2.vcf
