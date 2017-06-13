#!/bin/bash

# 1 --> reference
# 2 --> 

# Set working directory
DIR=/home/bimm185

# Input argument parsing
if [ $1 = "-h" ]
then
    echo "USAGE: BLT reference sample1forward sample1reverse sample2forward sample2reverse"

    exit 1
fi

if [ $# != 5 ] 
then
    echo "wrong number of files"
    echo "USAGE: BLT reference sample1forward sample1reverse sample2forward sample2reverse"

    exit 1
fi

mkdir $DIR/working_test
mkdir $DIR/working_test/fastqc

# for f in $(ls $DIR/data/*.fastq); do echo $f; fastqc -o $DIR/working_test/fastqc $f ; done

for f in $(ls $DIR/data/*.fastq); do echo $f; python $DIR/BIMM185FinalProject/modfastq.py $f ${f/data/working_test}; done



#bwa index -p $f{${f/reference/working_test}%%.fna}; done
cp $DIR/reference/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/.
for f in $(ls $DIR/working_test/*.fna); do echo $f; bwa index $f; done

# UC8
bwa mem $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC8_1.fastq $DIR/working_test/UC8_2.fastq > $DIR/working_test/UC8.sam
bwa mem $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC9_1.fastq $DIR/working_test/UC9_2.fastq > $DIR/working_test/UC9.sam

samtools view -o $DIR/working_test/UC8.bam -b $DIR/working_test/UC8.sam
samtools view -o $DIR/working_test/UC9.bam -b $DIR/working_test/UC9.sam

samtools flagstat $DIR/working_test/UC8.bam > $DIR/working_test/UC8_align_stats.txt
samtools flagstat $DIR/working_test/UC9.bam > $DIR/working_test/UC9_align_stats.txt

samtools sort $DIR/working_test/UC8.bam -o $DIR/working_test/UC8_sorted.bam 
samtools sort $DIR/working_test/UC9.bam -o $DIR/working_test/UC9_sorted.bam 

samtools mpileup -f $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC8_sorted.bam > $DIR/working_test/UC8.mpileup
samtools mpileup -f $DIR/working_test/GCF_000006945.2_ASM694v2_genomic.fna $DIR/working_test/UC9_sorted.bam > $DIR/working_test/UC9.mpileup

java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp $DIR/working_test/UC8.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > $DIR/working_test/UC8.vcf
java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp $DIR/working_test/UC9.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > $DIR/working_test/UC9.vcf



#VCF File comparison
#get header                                                                     

grep '#' $DIR/working/UC8.vcf > $DIR/working_test/VCF_header.txt
header=$(grep -c '#' $DIR/working/UC8.vcf)
header="$((header+1))"

echo 'number of headers'
echo $header
tail -n +$header $DIR/working_test/UC8.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > $DIR/working_test/UC8_sorted_vcf
tail -n +$header $DIR/working_test/UC9.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > $DIR/working_test/UC9_sorted_vcf

mkdir $DIR/working_test/vcf_output

join -1 1 -2 1 -o 1.3,1.5,1.6 $DIR/working_test/UC8_sorted_vcf $DIR/working_test/UC9_sorted_vcf > $DIR/working_test/vcf_output/common.vcf

join -1 1 -2 1 -v1 $DIR/working_test/UC8_sorted_vcf $DIR/working_test/UC9_sorted_vcf > $DIR/working_test/vcf_output/only_UC8.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' $DIR/working_test/vcf_output/only_UC8.vcf > $DIR/working_test/vcf_output/only_UC8_orig.vcf

cat $DIR/working_test/VCF_header.txt $DIR/working_test/vcf_output/only_UC8_orig.vcf > $DIR/working_test/vcf_output/only_UC8.vcf

join -1 1 -2 1 -v2 $DIR/working_test/UC8_sorted_vcf $DIR/working_test/UC9_sorted_vcf > $DIR/working_test/vcf_output/only_UC9.vcf

awk '{$(NF+1)=$1;$1=""}sub(FS,"")' $DIR/working_test/vcf_output/only_UC9.vcf > $DIR/working_test/vcf_output/only_UC9_orig.vcf

cat $DIR/working_test/VCF_header.txt $DIR/working_test/vcf_output/only_UC9_orig.vcf > $DIR/working_test/vcf_output/only_UC9.vcf
