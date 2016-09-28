#! /usr/bin/perl -w
use strict;

open(IN,"gen_list.txt");
my $com = "grep \"";
foreach(<IN>){
	my @t = split(" ",$_);
	chomp $t[1];
	$com = "$com$t[1]\\|";
	
}
chop $com;
chop $com;
$com = "$com\" ../resources/hg19_refGene.txt > bed_list.txt";
#print $com;
system($com);
close IN;

open(BED,"bed_list.txt");



foreach(<BED>){
	my @st;
	my @end;
	my @t;
    @t = split("\t",$_);
    @st = split(",",$t[9]);
    @end = split(",",$t[10]);
    
    my $i = 0;

	while($i < @st){
	    my $start = $st[$i]-20;
	    my $end = $end[$i]+20;
        print "$t[2]\t$start\t$end\t$t[1]\t$t[8]\t$t[3]\n";
        $i++;
	}
}


close BED;
exit;
