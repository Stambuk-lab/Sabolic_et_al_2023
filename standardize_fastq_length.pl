#!/usr/bin/perl

use strict;

print "opening gzcat  $ARGV[0] |\n";
my $s =  $ARGV[0];
my $p =  $ARGV[1];
open (IN, "zcat  $s |") || die "$!";
open (O, ">", $p ) || die "$!";

my $n = $ARGV[2];

my $head;
my $seq;
my $q;
my $t = 0;
while(<IN>){
  chomp;
  $t++;
  if(/\@/){
    if($head){
#     if(length($seq) != length($q)){
        my $first_n = substr($q, 0, $n);     #substr extracts a substring from string i.e. all characters from beggining of row (0) to $n (98) -> 98 reads left -> From Stacks manual:Stacks is optimized for short-read, Illumina-style sequencing. There is no limit to the length the sequences can be, although there is a hard-coded limit of 1024bp in the source code now for efficency reasons, but this limit could be raised if the technology warranted it.""
        $q = $first_n;
        $first_n =  substr($seq, 0, $n);
        $seq = $first_n;
#      }

    print O "$head\n$seq\n+\n$q\n";
    }
    $t = 1;
    ($head,$seq,$q)=($_,"","")
  }elsif($t==2){
    $seq = $_;
  }elsif($t == 4){
    $q = $_
  }
}

#if(length($seq) != length($q)){
        my $first_n = substr($q, 0, $n);
        $q = $first_n;
        $first_n =  substr($seq, 0, $n);
        $seq = $first_n;



print O "$head\n$seq\n+\n$q\n";

close IN;
close O;

print  $p . "\n";
system("gzip $p");
system("mv  $p.gz $s");
