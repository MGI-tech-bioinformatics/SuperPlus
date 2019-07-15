#!/usr/bin/bash

date
echo

TMPPATH=./tmp/
BINDIR=$(pwd)/lib/bin/
BWA=$(pwd)/bin/bwa
MINIMAP=$(pwd)/bin/minimap2
SAMTOOLS=$(pwd)/bin/samtools
BAMMDP=$(pwd)/bin/bammarkduplicates
BAMSORT=$(pwd)/bin/bamsort

#test if bin exist
if [ ! -e "$BWA" ];
then
	if [ hash bwa 2>/dev/null ];
	then
	$BWA=bwa
	else
	echo "Error: bwa not detected!"
	exit
	fi
fi

if [ ! -e "$MINIMAP" ];
then
        if [ hash minimap2 2>/dev/null ];
        then
        $MINIMAP=minimap2
        else
        echo "Error: minimap not detected!"
        exit
        fi
fi

if [ ! -e "$SAMTOOLS" ];
then
        if [ hash samtools 2>/dev/null ];
        then
        $SAMTOOLS=samtools
        else
        echo "Error: samtools not detected!"
        exit
        fi
fi

if [ ! -e "$BAMSORT" ];
then
        if [ hash bamsort 2>/dev/null ];
        then
        $BAMSORT=bamsort
        else
        echo "Error: bamsort not detected!"
        exit
        fi
fi



THREADS=$( expr $( nproc )  / 2 )
BUCKETS=$( expr $THREADS \* 2 / 3 )
MEM=640

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -i)
    INPUT1="$2"
    INPUT2="$3"
    shift # past argument
    shift # past value
    shift
    ;;
    -o)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -t)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    shift # past argument
    ;;
esac
done

BUCKETS=`expr $THREADS \* 2 / 3`

#echo $INPUT1
#echo $INPUT2
#echo $OUTPUT
#echo $TMPPATH

if [ "$TMPPATH" = "" ]; then 
	echo "Error! no temp path!"
	echo "Usage: ./runall_diploid.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]"
	exit
fi

if [ ! -e "$INPUT1" ]; then
        echo "Error! read1 fastq file not exist!"
	echo "Usage: ./runall_diploid.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]"
	exit
fi

if [ ! -e "$INPUT2" ]; then
        echo "Error! read2 fastq file not exist!"
	echo "Usage: ./runall_diploid.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]"
        exit
fi

if [ -e "$OUTPUT" ]; then
        echo "Error! output file already exist!"
	echo "Usage: ./runall_diploid.sh -i read1.fq read2.fq -o output.fasta [-t number_of_threads]"
        exit
fi

echo "Read file: $INPUT1 and $INPUT2, results saved in $OUTPUT, using $THREADS threads. "

$BINDIR/ParseBarcodedFastqs FASTQS="{$INPUT1,$INPUT2}" OUT_HEAD=$TMPPATH/reads NUM_THREADS=$THREADS NUM_BUCKETS=$BUCKETS

$BINDIR/DF ROOT=$TMPPATH LR=$TMPPATH/reads.fastb PIPELINE=cs ALIGN=False NUM_THREADS=$THREADS MAX_MEM_GB=$MEM

$BINDIR/CP DIR=$TMPPATH/GapToy/1/a.base NUM_THREADS=$THREADS MAX_MEM_GB=$MEM

$BINDIR/MakeFasta DIR=$TMPPATH/GapToy/1/a.base/final FLAVOR=pseudohap OUT_HEAD=$TMPPATH/original

$BINDIR/MakeFasta DIR=$TMPPATH/GapToy/1/a.base/final FLAVOR=pseudohap2 OUT_HEAD=$TMPPATH/original

zcat $TMPPATH/original.fasta.gz | sed -e "s/ style=3//g" | sed -e "s/ /_/g" >$TMPPATH/original_underscore.fasta

zcat $TMPPATH/original.1.fasta.gz | sed -e "s/ style=4//g" | sed -e "s/ /_/g" >$TMPPATH/original_underscore.1.fasta

zcat $TMPPATH/original.2.fasta.gz | sed -e "s/ style=4//g" | sed -e "s/ /_/g" >$TMPPATH/original_underscore.2.fasta

$BWA index $TMPPATH/original_underscore.fasta

$BWA mem -t $THREADS $TMPPATH/original_underscore.fasta $INPUT1 $INPUT2 | $SAMTOOLS view -Sb - > $TMPPATH/original.bam

$BAMSORT I=$TMPPATH/original.bam O=$TMPPATH/original_sort.bam blockmb=$( expr 1024 \* $MEM ) sortthreads=$THREADS

$BAMMDP I=$TMPPATH/original_sort.bam O=$TMPPATH/original_sort_mdp.bam markthreads=$THREADS index=1

$MINIMAP -x asm10 $TMPPATH/original_underscore.fasta $TMPPATH/original_underscore.fasta > $TMPPATH/overlap.paf

$BINDIR/ReScaffold_diploid FASTA_IN=$TMPPATH/original_underscore BAM_IN=$TMPPATH/original_sort_mdp.bam PAF_IN=$TMPPATH/overlap.paf LWML_IN=$TMPPATH/GapToy/1/a.base/records/lwml FASTA_OUT=$OUTPUT

date
echo
