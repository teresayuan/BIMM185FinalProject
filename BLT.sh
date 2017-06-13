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
	echo "blt.sh reference sample1forward sample1reverse sample2forward sample2reverse"

	exit 1
fi

echo "Welcome to BLT!"

s=$1
mkdir ${s%%.*}
# cd ${s%%.*}
mkdir ${s%%.*}/fast_qc

cp $1 ${s%%.*}
cp $2 ${s%%.*}
cp $3 ${s%%.*}
cp $4 ${s%%.*}
cp $5 ${s%%.*}

fastqc -o ${s%%.*}/fast_qc $2
fastqc -o ${s%%.*}/fast_qc $3
fastqc -o ${s%%.*}/fast_qc $4
fastqc -o ${s%%.*}/fast_qc $5

for f in $(ls $DIR/working_test/*.fna); do echo $f; bwa index $f; done

bwa index $1

s1=$(basename $2 .fastq)
s2=$(basename $4 .fastq)

bwa mem $1 $2 $3 > ${s%%.*}/$s1.sam
bwa mem $1 $4 $5 > ${s%%.*}/$s2.sam

samtools view -o ${s%%.*}/$s1.bam -b ${s%%.*}/$s1.sam
samtools view -o ${s%%.*}/$s2.bam -b ${s%%.*}/$s2.sam

samtools flagstat ${s%%.*}/$s1.bam > ${s%%.*}/$s1_align_stats.txt
samtools flagstat ${s%%.*}/$s2.bam > ${s%%.*}/$s2_align_stats.txt

samtools sort ${s%%.*}/$s1.bam -o ${s%%.*}/$s1_sorted.bam 
samtools sort ${s%%.*}/$s2.bam -o ${s%%.*}/$s2_sorted.bam 

samtools mpileup -f $1 ${s%%.*}/$s1_sorted.bam > ${s%%.*}/$s1.mpileup
samtools mpileup -f $1 ${s%%.*}/$s2_sorted.bam > ${s%%.*}/$s2.mpileup

java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp ${s%%.*}/$s1.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > ${s%%.*}/$s1.vcf
java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp ${s%%.*}/$s2.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > ${s%%.*}/$s2.vcf


#VCF File comparison
grep '#' ${s%%.*}/$s1.vcf > ${s%%.*}/VCF_header.txt
header=$(grep -c '#' ${s%%.*}/$s1.vcf)
header="$((header+1))"

echo 'number of headers'
echo $header
tail -n +$header ${s%%.*}/$s1.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > ${s%%.*}/$s1_sorted_vcf
tail -n +$header ${s%%.*}/$s2.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > ${s%%.*}/$s2_sorted_vcf

mkdir ${s%%.*}/vcf_output

join -1 1 -2 1 -o 1.3,1.5,1.6 ${s%%.*}/$s1_sorted_vcf ${s%%.*}/$s2_sorted_vcf > ${s%%.*}/vcf_output/common.vcf

join -1 1 -2 1 -v1 ${s%%.*}/$s1_sorted_vcf ${s%%.*}/$s2_sorted_vcf > ${s%%.*}/vcf_output/only_$s1.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' ${s%%.*}/vcf_output/only_$s1.vcf > ${s%%.*}/vcf_output/only_$s1_orig.vcf

cat ${s%%.*}/VCF_header.txt ${s%%.*}/vcf_output/only_$s1_orig.vcf > ${s%%.*}/vcf_output/only_$s1.vcf

join -1 1 -2 1 -v2 ${s%%.*}/$s1_sorted_vcf ${s%%.*}/$s2_sorted_vcf > ${s%%.*}/vcf_output/only_$s2.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' ${s%%.*}/vcf_output/only_$s2.vcf > ${s%%.*}/vcf_output/only_$s2_orig.vcf

cat ${s%%.*}/VCF_header.txt ${s%%.*}/vcf_output/only_$s2_orig.vcf > ${s%%.*}/vcf_output/only_$s2.vcf
