# SuperPlus

**SuperPlus** is an stLFR assembly pipeline based on Supernova, with a fastq parsing and sorting module constumized for stLFR data, a re-scaffolding module and an ONT reads gap-closing module. 

### Dependencies:

* SOAPfilter_v2.2(already in ./split_barcode/ directory)
* gcc (tested on version 4.8.5)
* bwa
* biobambam(bammarkduplicates and bamsort)
* samtools
* minimap2

### Installation: 

To install SuperPlus, run:
```
make
```
### Configuration:

For configuring, put bwa, biobambam, samtools and minimap2 in ./bin/ or system path, or modify the runall.sh and runall_ont.sh to specify the path to find them.

### Quick Guide:

Before SuperPlus pipeline, you need to split barcodes, filter duplicate and adpater for stLFR reads. You can run run_split_unsort_filter.sh under ./split_barcode/. And you need to change program directory, Read1, Read2 and Output.
```
sh ./run_split_unsort_filter.sh
```

The SuperPlus pipeline requires two stLFR read file and an optional ONT read file, all formatted in fastq. runall.sh and runall_ont.sh are two scripts for running pipeline.

* Need to change directory for program.

For monoploid:

If you do not have ONT reads, run
```
./runall.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]
```
and if you have ONT reads, run
```
./runall_ont.sh -i read1.fq read2.fq ont.fq -o output.fasta [-t number_of_threads]
```

For diploid:

If you do not have ONT reads, run
```
./runall_diploid.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]
```
and if you have ONT reads, run
```
./runall_ont_diploid.sh -i read1.fq read2.fq ont.fq -o output.fasta [-t number_of_threads]
```

There are three steps to the pipeline:

1. Run Supernova to generate initial scaffolds.

2. Map scaffolds to themself to remove redundant scaffolds and map reads to initial scaffolds to do re-scaffolding based on common barcodes.

3. Use ONT reads to fill gaps in scaffolds, if you do not have ONT reads, skip this step.

All the temporary files will be stored in folder ./tmp. 

By default the memory limit of SuperPlus is 640GB, if you want to change it, modify the line 
```
MEM=640
```
in scripts.

And for change the number of threads, modify the line 
```
THREADS=40
```
in scripts. By default the parameter *BUCKETS* is 2/3 of number of threads, which works well for most of our case.


### Usage:

You can modify the script runall.sh or runall_ont.sh for costumization.

```
$BINDIR/ParseBarcodedFastqs FASTQS="{$INPUT1,$INPUT2}" OUT_HEAD=$TMPPATH/reads NUM_THREADS=$THREADS NUM_BUCKETS=$BUCKETS

$BINDIR/DF ROOT=$TMPPATH LR=$TMPPATH/reads.fastb PIPELINE=cs ALIGN=False NUM_THREADS=$THREADS MAX_MEM_GB=$MEM

$BINDIR/CP DIR=$TMPPATH/GapToy/1/a.base NUM_THREADS=$THREADS MAX_MEM_GB=$MEM

$BINDIR/MakeFasta DIR=$TMPPATH/GapToy/1/a.base/final FLAVOR=pseudohap OUT_HEAD=$TMPPATH/original
```

Those four lines runs supernova 1.1, and generate an initial assembly at **$TMPPATH/original.fasta.gz**. In most cases you do not need to change parameters in them. 


Then we run external scaffolding. First we run bwa to align stLFR reads to assemblies.

```
zcat $TMPPATH/original.fasta.gz | sed -e "s/ /_/g" >$TMPPATH/original_underscore.fasta

$BWA index $TMPPATH/original_underscore.fasta

$BWA mem -t $THREADS $TMPPATH/original_underscore.fasta $INPUT1 $INPUT2 | $SAMTOOLS view -Sb > $TMPPATH/original.bam
```

Then sort bam and markdup

```
$BAMSORT I=$TMPPATH/original.bam O=$TMPPATH/original_sort.bam blockmb=$( expr 1024 \* $MEM ) sortthreads=$THREADS

$BAMMDP I=$TMPPATH/original_sort.bam O=$TMPPATH/original_sort_mdp.bam markthreads=$THREADS index=1
```

We also align the assembly to itself for removing duplicates in outputs

```
$MINIMAP -x asm10 $TMPPATH/original_underscore.fasta $TMPPATH/original_underscore.fasta > $TMPPATH/overlap.paf
```

Finally do rescaffolding to produce final ouputs.

```
$BINDIR/ReScaffold FASTA_IN=$TMPPATH/original_underscore.fasta BAM_IN=$TMPPATH/original_sort_mdp.bam PAF_IN=$TMPPATH/overlap.paf LWML_IN=$TMPPATH/GapToy/1/a.base/records/lwml FASTA_OUT=$OUTPUT
```

If you have any questions, please contact wupei1@genomics.cn.

