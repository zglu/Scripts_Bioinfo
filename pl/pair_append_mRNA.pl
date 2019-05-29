#!/usr/bin/perl
# usage: perl ~zl3/scripts/pair_append_mRNA.pl [id-note file] [INPUT gff] > OUTPUTFILE
# currently this scripts creates a blank line after each matched line

use strict;
use warnings;

my %hsh=();

open (MYFILE, $ARGV[0]); # Smp_121100.1	;product=Hypothetical protein
open (MYFILE1, $ARGV[1]);

while (<MYFILE>) {
my@arr = split/\t+/; # tab separated
$hsh{$arr[0]} = $arr[1];
}
my $flag;
while(<MYFILE1>)
{
$flag=0;
my $line=$_;
foreach my $key (keys %hsh)
{
   if($line=~/ID=$key;/) # eg. ID=Smp_123450.1;Parent=Smp_123450
   {
    $flag=1;
    $line=~s/$/$hsh{$key}/g; # ID=Smp_123450.1;Parent=Smp_123450 ---> ID=Smp_123450.1;Parent=Smp_123450;product=Hypothetical protein...
    print $line;
   }
}
  if($flag!=1)
  {
  print $line;
  $flag=0;
  }
}
close(MYFILE);
close(MYFILE1);
