#!/usr/bin/perl

# usage: perl ~zl3/scripts/pair_append.pl [id pairs] [gff input] > OUTPUTFILE
# currently this scripts creates a blank line after each matched line

use strict;
use warnings;

my %hsh=();

open (MYFILE, $ARGV[0]); # Smp_121100	;previous_systematic_id=...
open (MYFILE1, $ARGV[1]); # gff3 file

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
   if($line=~/ID=$key$/) # ID=Smp_121100 OR Parent=Smp_121100 for mRNA
   {
    $flag=1;
    $line=~s/$key$/$key$hsh{$key}/g; # ID=Smp_121100 ---> ID=Smp_121100;previous_systematic_id=...
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
