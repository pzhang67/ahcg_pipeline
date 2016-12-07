#echo "Download bam file of patient1"
#wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
#wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai
#echo "Finished downloading bam file of patient1"

echo "Perform variant call analysis using GATK-HaplotyeCaller"
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o p1_raw_variants.vcf
echo "Finished performing variant call analysis using GATK-HaplotyeCaller"
echo "Perform variants recalibration analysis using GATK-VariantRecalibrator"
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches -o p1_recalibrated_snps_raw_indels.vcf
echo "Finished perform variants recalibration analysis using GATK-VariantRecalibrator"

perl scripts/grep_gen_list.pl
echo "dcm-related 6 genes' genome position bed file created."

bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > patient1_dcm_final.vcf
echo "extract the dcm-related 6 genes' genomic variants"

samtools view -L scripts/temp/exom_list.bed Patient1_RG_MD_IR_BQ.bam -b > patient1.dcm.bam
echo "Depth for each gene"

bedtools genomecov -ibam patient1.dcm.bam -bga > patient1.dcm.bed
echo "Reporting genome coverage"

bedtools intersect -split -a patient1.dcm.bed -b scripts/temp/exom_list.bed -bed > scripts/temp/patient1.final.bed
echo "Overlaping"

python scripts/depth_for_each_site.py
echo "Depth for each site"

perl scripts/label.pl scripts/temp/patient1.final.full.txt patient1.final.full.labeled.txt
cp patient1.final.full.labeled.txt scripts/txt/p1.final.txt
rm patient1.final.full.labeled.txt
echo "Depth labeled"

R CMD BATCH '--args scripts/txt/p1.final.txt' scripts/plots.R
cp boxplot.exons.LMNA.jpg SNP_report/boxplots
cp boxplot.exons.MYBPC3.jpg SNP_report/boxplots
cp boxplot.exons.MYH6.jpg SNP_report/boxplots
cp boxplot.exons.MYH7.jpg SNP_report/boxplots
cp boxplot.exons.SCNSA.jpg SNP_report/boxplots
cp boxplot.exons.TNNT2.jpg SNP_report/boxplots
rm boxplot.exons.LMNA.jpg boxplot.exons.MYBPC3.jpg boxplot.exons.MYH6.jpg boxplot.exons.MYH7.jpg boxplot.exons.SCNSA.jpg boxplot.exons.TNNT2.jpg
cp cutoff.MYH7.csv SNP_report/cutoffs
cp cutoff.MYH6.csv SNP_report/cutoffs
cp cutoff.LMNA.csv SNP_report/cutoffs
cp cutoff.TNNT2.csv SNP_report/cutoffs
cp cutoff.SCNSA.csv SNP_report/cutoffs
cp cutoff.MYBPC3.csv SNP_report/cutoffs
rm cutoff.MYH7.csv cutoff.MYH6.csv cutoff.LMNA.csv cutoff.TNNT2.csv cutoff.SCNSA.csv cutoff.MYBPC3.csv
echo "Cutoff for sequencing result created"

cp P1cutoff.final.csv scripts/cutoff/P1cutoff.final.csv
rm P1cutoff.final.csv
cp patient1_dcm_final.vcf scripts/vcf/P1_afterinter.vcf
python scripts/min.py
cp scripts/vcf/*new.vcf SNP_report/
rm scripts/vcf/*new.vcf
echo "SNP report created"

