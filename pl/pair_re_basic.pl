#!/usr/bin/perl

# usage: perl ~zl3/scripts/pair_re_basic.pl [Aug-RATT id pairs] [INPUT FILE] > OUTPUTFILE
# replace simply just the whole word

use strict;
use warnings;

my %hsh=();

open (MYFILE, $ARGV[0]);
open (MYFILE1, $ARGV[1]);

while (<MYFILE>) {
my@arr = split/\s+/;
$hsh{$arr[0]} = $arr[1];
}
my $flag;
while(<MYFILE1>)
{
$flag=0;
my $line=$_;
foreach my $key (keys %hsh)
{
   if($line=~/$key/)
   {
    $flag=1;
    $line=~s/$key/$hsh{$key}/g; # replace each string with paired-string
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
