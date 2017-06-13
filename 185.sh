#!/bin/bash

# Input argument parsing
if [ $1 = "-h" ]
then
    echo "USAGE: BLT reference.fasta data_dir sample1root sample2root"

    exit 1
fi

if [ $# != 4 ] 
then
    echo "wrong number of files"
    echo "USAGE: BLT reference.fasta data_dir sample1root sample2root"

    exit 1
fi

# Set working directory and other variables
HOMEDIR=/home/`whoami`
REFERENCE=$1
DATADIR=$2
SAMPLE1=$3
SAMPLE2=$4
HOMEDIR=${HOMEDIR%%\/}
DATADIR=${DATADIR%%\/}

echo "Welcome to BLT!"

mkdir $HOMEDIR/working
mkdir $HOMEDIR/working/fastqc

echo "===========Starting FastQC==========="
for f in $(ls $DATADIR/*.fastq); do echo "Performing fastqc on $f..."; fastqc -o $HOMEDIR/working/fastqc $f ; done
echo "===========Done with FastQC==========="


echo "===========Performing Alignment==========="
# Copy all required files to the working directory
cp $REFERENCE $HOMEDIR/working/
cp $DATADIR/*.fastq $HOMEDIR/working/
REFERENCE=${REFERENCE##*/}

pushd $HOMEDIR/working/

# Index the reference
bwa index $REFERENCE

echo "Align using bwa mem..."
bwa mem $REFERENCE ${SAMPLE1}_1.fastq ${SAMPLE1}_2.fastq > ${SAMPLE1}.sam
bwa mem $REFERENCE ${SAMPLE2}_1.fastq ${SAMPLE2}_2.fastq > ${SAMPLE2}.sam

echo "Convert to bam..."
samtools view -o ${SAMPLE1}.bam -b ${SAMPLE1}.sam
samtools view -o ${SAMPLE2}.bam -b ${SAMPLE2}.sam

echo "Get alignment stats..."
samtools flagstat ${SAMPLE1}.bam > ${SAMPLE1}_align_stats.txt
samtools flagstat ${SAMPLE2}.bam > ${SAMPLE2}_align_stats.txt

echo "Sort the alignment..."
samtools sort ${SAMPLE1}.bam -o ${SAMPLE1}_sorted.bam 
samtools sort ${SAMPLE2}.bam -o ${SAMPLE2}_sorted.bam 

popd

echo "===========Done with Alignment==========="

echo "===========Creating mpileup==========="
pushd $HOMEDIR/working/

samtools mpileup -f $REFERENCE ${SAMPLE1}_sorted.bam > ${SAMPLE1}.mpileup
samtools mpileup -f $REFERENCE ${SAMPLE2}_sorted.bam > ${SAMPLE2}.mpileup

popd
echo "===========Done with mpileup==========="

echo "===========Calling Variants==========="
pushd $HOMEDIR/working/

java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp ${SAMPLE1}.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > ${SAMPLE1}.vcf
java -Xmx4g -jar ~/VarScan.v2.4.3.jar mpileup2snp ${SAMPLE2}.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > ${SAMPLE2}.vcf

popd
echo "===========Done with Variant Calling==========="



echo "===========Analyzing SNPs==========="
#VCF File comparison
#get header                                                                     

pushd $HOMEDIR/working/
mkdir tmp
mkdir vcf_output

# Assume VCF headers are similar, so save an arbitrary one
grep '#' ${SAMPLE1}.vcf > tmp/VCF_header.txt
header=$(grep -c '#' ${SAMPLE1}.vcf)
header="$((header+1))"

tail -n +$header ${SAMPLE1}.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > ${SAMPLE1}_sorted_vcf
tail -n +$header ${SAMPLE2}.vcf | awk '{printf("%s:%s:%s %s %s %s %s %s %s %s %s %s %s\n", $2, $4, $5, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10);}' | sort > ${SAMPLE2}_sorted_vcf

join -1 1 -2 1 -o 1.3,1.5,1.6 ${SAMPLE1}_sorted_vcf ${SAMPLE2}_sorted_vcf > vcf_output/common.vcf
join -1 1 -2 1 -v1 ${SAMPLE1}_sorted_vcf ${SAMPLE2}_sorted_vcf > vcf_output/only_${SAMPLE1}.vcf
awk '{$(NF+1)=$1;$1=""}sub(FS,"")' vcf_output/only_${SAMPLE1}.vcf > tmp/only_${SAMPLE1}_orig.vcf
cat tmp/VCF_header.txt tmp/only_${SAMPLE1}_orig.vcf > vcf_output/only_${SAMPLE1}.vcf
join -1 1 -2 1 -v2 ${SAMPLE1}_sorted_vcf ${SAMPLE2}_sorted_vcf > vcf_output/only_${SAMPLE2}.vcf
awk '{$(NF+1)=$1;$1=""}sub(FS,"")' vcf_output/only_${SAMPLE2}.vcf > tmp/only_${SAMPLE2}_orig.vcf
cat tmp/VCF_header.txt tmp/only_${SAMPLE2}_orig.vcf > vcf_output/only_${SAMPLE2}.vcf

rm -rf tmp

popd
echo "===========Done analyzing SNPs==========="
