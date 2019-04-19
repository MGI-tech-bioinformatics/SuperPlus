use strict;

if(@ARGV != 6)
{
        print "Example: perl Split_SingleTube_reads.pl barcode.list barcode_RC.list read_1.fq.gz read_2.fq.gz 100 split_read \n";
        exit(0);
}

my $read_len = $ARGV[4];
my ($n1, $n2, $n3, $n4, $n5) = (10, 6, 10, 6, 10);
my %barcode_hash;
my %barcode_RC_hash;

open IN,"$ARGV[0]" or die "cann't not open barcode.list";
my $n = 0;
while(<IN>){
  $n ++;
  my @line = split;
  my @barcode = split(//,$line[0]);
  my $barcode_ID = $line[1];
  for(my $num = 0; $num <= 9; $num++){
    my @barcode_mis = @barcode;
    $barcode_mis[$num] = "A";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "G";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "C";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "T";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
  }
}
close IN;
my $barcode_types = $n * $n *$n;
my $barcode_each = $n;

open IN_RC,"$ARGV[1]" or die "cann't not open barcode_RC.list";
while(<IN_RC>){
  my @line = split;
  my @barcode = split(//,$line[0]);
  my $barcode_ID = $line[1];
  for(my $num = 0; $num <= 9; $num++){
    my @barcode_mis = @barcode;
    $barcode_mis[$num] = "A";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_RC_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "G";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_RC_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "C";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_RC_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "T";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_RC_hash{$barcode_mis} = $barcode_ID;
  }
}
close IN_RC;

open IN1,"gzip -dc $ARGV[2] |" or die "cannot open file";
open IN2,"gzip -dc $ARGV[3] |" or die "cannot open file";
$n = 0;
my $reads_num;
my $progress;
my %index_hash;
my %index_hash_reverse;
my $split_barcode_num;
my (@Read1, @Read2, @Read3, @Read4);
my $T;
my $id;
my $reads_num;
my @line;
my @Read_num;
$Read_num[0] = 0;
my $split_reads_num;


while(<IN2>){
  chomp;
  @line = split;
  $n ++;
  if($n % 4 == 1){
    $reads_num ++;
    my @A  = split(/\//,$line[0]);
         $id = $A[0];
         if($reads_num % 1000000 == 1)
         {
              print "reads_1 processed $progress (M) reads ...\n";
              $progress ++;
         }

  }
  if($n % 4 == 2){
    my $b1 = substr($line[0], $read_len, $n1);
    my $b2 = substr($line[0], $read_len+$n1+$n2, $n3);
    my $b3 = substr($line[0], $read_len+$n1+$n2+$n3+$n4, $n5);
    if((exists $barcode_hash{$b1}) && (exists $barcode_hash{$b2}) && (exists $barcode_hash{$b3})){
      my $hash = $barcode_hash{$b1}."_".$barcode_hash{$b2}."_".$barcode_hash{$b3};
      if(!(exists $index_hash{$hash})){
        $split_barcode_num ++;
        $index_hash{$hash} = $split_barcode_num;
        $index_hash_reverse{$split_barcode_num} = $hash;
        $Read_num[$index_hash{$hash}] = 0;
      }
      $split_reads_num ++;
      $Read_num[$index_hash{$hash}] ++;
      $T = <IN1>; chomp($T);
      $Read1[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $id."\#$hash\/1\t$index_hash{$hash}\t1";
      $T = <IN1>; chomp($T);
      $Read2[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
      $T = <IN1>; chomp($T);
      $Read3[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
      $T = <IN1>; chomp($T);
      $Read4[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
    }
    else{
      $Read_num[0] ++;
      $T = <IN1>; chomp($T);
      $Read1[0][$Read_num[0]] = $id."\#0_0_0\/1\t0\t1";
      $T = <IN1>; chomp($T);
      $Read2[0][$Read_num[0]] = $T;
      $T = <IN1>; chomp($T);
      $Read3[0][$Read_num[0]] = $T;
      $T = <IN1>; chomp($T);
      $Read4[0][$Read_num[0]] = $T;
    }

  }
}
close IN1;
close IN2;

open OUT2, ">split_stat_read1.log" or die "Can't write file";
open OUT, "| gzip > $ARGV[5].1.fq.gz" or die "Can't write file";
print OUT2 "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
my $r;
$r = 100 *  $split_barcode_num/$barcode_types;
print OUT2 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 *  $split_reads_num/$reads_num;
print OUT2 "Reads_pair_num  = $reads_num \n";
print OUT2 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
  print OUT2 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
  for(my $j=1;$j<=$Read_num[$i];$j++){
    print OUT "$Read1[$i][$j]\n$Read2[$i][$j]\n$Read3[$i][$j]\n$Read4[$i][$j]\n";
  }
}

for(my $j=1;$j<=$Read_num[0];$j++){
    print OUT "$Read1[0][$j]\n$Read2[0][$j]\n$Read3[0][$j]\n$Read4[0][$j]\n";
  }


close OUT;
close OUT2;

print "all done!\n";







