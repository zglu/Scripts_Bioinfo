#!/usr/bin/perl

# usage: perl ~zl3/scripts/pair_replace.pl [Aug-RATT id pairs] [gff input] > OUTPUTFILE
# replace the ids in GFF file with Smp number (replaces g10 to Smp_XXXX, and g10.t1 to Smp_XXXX.1)

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
   if($line=~/$key($|\.t)/) # finding lines with string at line end or with string + ".t"
   {
    $flag=1;
    $line=~s/$key$/$hsh{$key}/g; # replace string at line end with paired-string
    $line=~s/$key\.t/$hsh{$key}\./g; # replace string + ".t" with paired-string + "."
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
