#! /usr/bin/perl -w
use strict;

open(IN,"NM_007294.txt");

my @st;
my @end;
my @t;

foreach(<IN>){
	@t = split("\t",$_);
	@st = split(",",$t[9]);
	@end = split(",",$t[10]);
}

my $i = 0;

while($i < @st){
	print "$t[2]\t$st[$i]\t$end[$i]\t$t[1]\t$t[8]\t$t[3]\n";
	$i++;
}
close IN;
exit;