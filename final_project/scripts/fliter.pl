 
my %db;
my %snp;

sub score{
	if($_[0] eq "Mutation confirmed"){
		return 1;
	}elsif($_[0] eq "Ambiguous sequencing result"){
		return 0.6;
	}elsif($_[0] eq "Mutation untrusted"){
		return 0;
	}
}

sub loaddb{
	open(IN,$_[0]);
	while(<IN>){
		my @t = split("\t",$_);
		my $str = $t[0]."-".$t[1]."-".$t[3]."-".$t[4];
		if(!$db{$str}){
			$db{$str} = score($t[10]);
		}else{
			$db{$str} = $snp{$str} + score($t[10]);
		}
	}
	close IN;
}

sub overlap{
	open(IN,$_[0]);
	while(<IN>){
		my @t = split("\t",$_);
		my $str = $t[0]."-".$t[1]."-".$t[3]."-".$t[4];
		if(!$snp{$str}){
			$snp{$str} = score($t[10]);
		}else{
			$snp{$str} = $snp{$str} + score($t[10]);
		}
	}
	close IN;
}

loaddb("SNP_report/C1_afterinter.vcfnew.vcf");
loaddb("SNP_report/C2_afterinter.vcfnew.vcf");

overlap("SNP_report/P1_afterinter.vcfnew.vcf");
overlap("SNP_report/P2_afterinter.vcfnew.vcf");
overlap("SNP_report/P3_afterinter.vcfnew.vcf");
overlap("SNP_report/P4_afterinter.vcfnew.vcf");

open(OUT,">SNP_report/SNP_report.txt");

#print OUT "chr\tsite\toriganl\tmutation\tcontrol_score\tSNP_score\n";

foreach(keys %snp){
	if(!$db{$_}){
		my @t = split("-",$_);
		$, = "\t";
		print OUT @t;
		print OUT "\t".-1;
		print OUT "\t".$snp{$_}."\n";
	}elsif($db{$_} < 2){
		my @t = split("-",$_);
		$, = "\t";
		print OUT @t;
		print OUT "\t".$db{$_};
		print OUT "\t".$snp{$_}."\n";
	};
}

my $com = "sort -k5n -k6nr -o SNP_report/SNP_report.txt SNP_report/SNP_report.txt";

system($com);

exit;
