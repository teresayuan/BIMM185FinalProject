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


for f in $(ls $DIR/data/*.fastq); do echo $f; python modfastq.py $f ${f/data/working_test}; done