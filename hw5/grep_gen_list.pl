#! /usr/bin/perl -w
use strict;

open(IN,"gen_list.txt");
my $com = "grep \"";
foreach(<IN>){
	if($_ =~ /(NM_.+)\.[0-9]/){
		$com = "$com$1\\|";
	}
}
chop $com;
chop $com;
$com = "$com\" ../resources/hg19_refGene.txt > bed_list.bed";
#print $com;
system($com);
close IN;

open(BED,"bed_list.bed");



foreach(<BED>){
	my @st;
	my @end;
	my @t;
    @t = split("\t",$_);
    @st = split(",",$t[9]);
    @end = split(",",$t[10]);
    
    my $i = 0;

	while($i < @st){
        print "$t[2]\t$st[$i]\t$end[$i]\t$t[1]\t$t[8]\t$t[3]\n";
        $i++;
	}
}


close BED;
exit;