#!/bin/bash

#VCF File comparison
#get header                                                                     

DIR=/home/bimm185 

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
