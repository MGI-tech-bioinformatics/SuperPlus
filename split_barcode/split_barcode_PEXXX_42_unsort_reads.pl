use strict;

if(@ARGV != 6)
{
        print "Example: perl Split_SingleTube_reads.1.pl barcode.list barcode_RC.list read_1.fq.gz read_2.fq.gz 100 split_read \n";
        exit(0);
}

my $read_len = $ARGV[4];
my ($n1, $n2, $n3,$n4,$n5) = (10, 6, 10,6,10);
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

open OUT, "| gzip > $ARGV[5].1.fq.gz" or die "Can't write file";
open OUT3, "| gzip > $ARGV[5].2.fq.gz" or die "Can't write file";

if ($ARGV[2] =~ /gz/){
  open IN1,"gzip -dc $ARGV[2] |" or die "cannot open file";
  open IN2,"gzip -dc $ARGV[3] |" or die "cannot open file";
}else{
  open IN1,"$ARGV[2]" or die "cannot open file";
  open IN2,"$ARGV[3]" or die "cannot open file"; 
}
$n = 0;
my $reads_num;
my $progress;
my %index_hash;
my %index_hash_reverse;
my $split_barcode_num;
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
    my $read = substr($line[0], 0, $read_len);
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
      print OUT "$id\#$hash\/1\t$index_hash{$hash}\t1\n";
      print OUT3 "$id\#$hash\/2\t$index_hash{$hash}\t1\n";
      print OUT3 "$read\n";
#      $Read1[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $id."\#$hash\/1\t$index_hash{$hash}\t1";
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
      $T = <IN2>; $n++; chomp($T);
      print OUT3 "$T\n";
      $T = <IN2>; $n++; chomp($T);
      my $qual = substr($T,0,$read_len);
      print OUT3 "$qual\n";
#      $Read2[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
#      $Read3[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
#      $Read4[$index_hash{$hash}][$Read_num[$index_hash{$hash}]] = $T;
    }
    else{
      $Read_num[0] ++;
      print OUT3 "$id\#0_0_0\/2\t0\t1\n";
      print OUT3 "$read\n";
      $T = <IN2>; $n++;chomp($T);
      print OUT3 "$T\n";
      $T = <IN2>; $n++;chomp($T);
      my $qual = substr($T,0,$read_len);
      print OUT3 "$qual\n";


      $T = <IN1>; chomp($T);
      print OUT "$id\#0_0_0\/1\t0\t1\n";
#      $Read1[0][$Read_num[0]] = $id."\#0_0_0\/1\t0\t1";
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
#      $Read2[0][$Read_num[0]] = $T;
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
#      $Read3[0][$Read_num[0]] = $T;
      $T = <IN1>; chomp($T);
      print OUT "$T\n";
#      $Read4[0][$Read_num[0]] = $T;
    }

  }
}
close IN1;
close IN2;
close OUT3;

open OUT2, ">split_stat_read1.log" or die "Can't write file";
print OUT2 "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
my $r;
$r = 100 *  $split_barcode_num/$barcode_types;
print OUT2 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 *  $split_reads_num/$reads_num;
print OUT2 "Reads_pair_num  = $reads_num \n";
print OUT2 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
  print OUT2 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}


close OUT;
close OUT2;

print "all done!\n";







