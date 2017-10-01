# BIMM185 Final Project

BLT.sh is a bash script to automate alignment and variant calling of two 
samples to a reference genome. 

###Usage
To run the pipeline, run the command: 

```sh
BLT.sh </path/to/reference.fasta> </path/to/data_directory> <sample1root> <sample2root>
```

The two samples can be single-end or paired-end reads. If the samples are 
single-end reads, they are expected to be named as:

```sh
<sample1root>.fastq
<sample2root>.fastq
```

Paired-end reads are expected to be named as:

```sh
<sample1root>_1.fastq
<sample1root>_2.fastq
<sample2root>_1.fastq
<sample2root>_2.fastq
```

This script uses FASTQC, BWA, samtools, and VarScan to analyze the two samples.
