# SuperPlus

**SuperPlus** is an stLFR assembly pipeline based on Supernova, with a fastq parsing and sorting module constumized for stLFR data, a re-scaffolding module and an ONT reads gap-closing module. 

### Dependencies:

* gcc (tested on version 4.8.5)
* bwa
* biobambam
* samtools
* minimap2

### Installation: 

To install SuperPlus, run:
```
make
```
### Configuration:

For configuring, put bwa, biobambam, samtools and minimap2 in ./bin/ or system path, or modify the runall.sh and runall_ont.sh to specify the path to find them.

### Usage:

The SuperPlus pipeline requires two stLFR read file and an optional ONT read file, all formatted in fastq. runall.sh and runall_ont.sh are two scripts for running pipeline.

If you do not have ONT reads, run
```
./runall.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]
```
and if you have ONT reads, run
 ```
./runall_ont.sh -i read1.fq read2.fq ont.fq -o output.fasta [-t number_of_threads]
```

There are three steps to the pipeline:

1. Run Supernova to generate initial scaffolds.

2. Map scaffolds to themself to remove redundant scaffolds and map reads to initial scaffolds to do re-scaffolding based on common barcodes.

3. Use ONT reads to fill gaps in scaffolds, if you do not have ONT reads, skip this step.

All the temporary files will be stored in folder ./tmp. By default the memory limit of SuperPlus is 640GB, if you want to change it, modify the line 
```
MEM=640
```
in scripts.

If you have any questions, please contact wupei1@genomics.cn.

