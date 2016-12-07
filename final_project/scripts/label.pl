#! /usr/bin/perl -w
use strict;

my $inname = $ARGV[0];
my $outname = ">".$ARGV[1];

my @exomchr;
my @exomstart;
my @exomend;
my @exomname;
my @exomlabel;
my @label = (0,0,0,0,0,0);
my %index;


open(BED,"scripts/temp/exom_list.bed");

my $gene = 0;
my $exom = 1;
foreach(<BED>){
    my @t = split("\t",$_);
    if(!$index{$t[3]}){
        $exom = 1;
        $gene++;
        $index{$t[3]} = $gene;
        push(@exomchr,$t[0]);
        push(@exomstart,$t[1]);
        push(@exomend,$t[2]);
        push(@exomname,$t[3]);
        $label[$gene-1] = $exom;
        my $label = "$label[0]\t$label[1]\t$label[2]\t$label[3]\t$label[4]\t$label[5]";
        push(@exomlabel,$label);
       @label = (0,0,0,0,0,0);
    }else{
        push(@exomchr,$t[0]);
        push(@exomstart,$t[1]);
        push(@exomend,$t[2]);
        push(@exomname,$t[3]);
        $label[$gene-1] = $exom;
        my $label = "$label[0]\t$label[1]\t$label[2]\t$label[3]\t$label[4]\t$label[5]";
        push(@exomlabel,$label);
        @label = (0,0,0,0,0,0);
    }
    $exom++;
}

close BED;

open(SITE,$inname);
open(OUT,$outname);
print OUT "chr\tsite\tdepth\tLMNA\tTNNT2\tSCNSA\tMYBPC3\tMYH6\tMYH7\n";
my $n = 0;
=pod
foreach(@exomchr){
    print $n."\t".$exomchr[$n]."\t".$exomstart[$n]."\t".$exomend[$n]."\t".$exomname[$n]."\t".$exomlabel[$n]."\n";
    $n++;
}

$n = 0;
=cut
while(<SITE>){
    chomp;
    my @t = split("\t",$_);
    RS:
    if($t[0] eq $exomchr[$n]){
        if($t[1] >= $exomstart[$n] and $t[1] < $exomend[$n]){
            print OUT $_;
            print OUT "\t".$exomlabel[$n]."\t".$exomname[$n]."\n";
        }elsif($t[1] < $exomstart[$n]){
            print "next";
        }else{
            $n++;
            goto RS;
        }
    }else{
        $n++;
        goto RS;
    }
}

close SITE;
close OUT;

exit;





