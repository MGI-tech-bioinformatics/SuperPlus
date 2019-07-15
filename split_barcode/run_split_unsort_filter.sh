date
Reads1=./read_1.fq
Reads2=./read_2.fq
Output=output

time perl ./split_barcode_PEXXX_42_unsort_reads.pl ./barcode.list ./barcode_RC.list $Reads1 $Reads2 100 $Output 

echo "split barcode have done"

echo "$Output.1.fq.gz" >lane.lst

echo "$Output.2.fq.gz" >>lane.lst

#WARNING: this is the old version adapter filter.
#time ./SOAPfilter_v2.2 -t 50 -q 33 -y -F CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG -R TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC -p -M 2 -f -1 -Q 10 lane.lst stat.txt 1>log 2>err

time ./SOAPfilter_v2.2 -t 50 -q 33 -y -F CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA -R TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG -p -M 2 -f -1 -Q 10 lane.lst stat.txt 1>log 2>err

perl -e '@A;$n=-1; while(<>){$n++;chomp;@t=split; for($i=0;$i<@t;$i++){$A[$n][$i]=$t[$i]; }} for($i=0;$i<@t;$i++){print "$A[0][$i]\t$A[1][$i]\n";}' stat.txt >stat.csv

date
echo Filter reads done!


