#! /usr/bin/perl -w
use strict;

open(IN,"NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf");

while(<IN>){
    if($_ =~ /^\#/){
	print;
    }else{
	print "chr$_";
    }
}

close IN;
exit;
