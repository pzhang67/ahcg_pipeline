# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878
```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

### Set up github
- Changing a remote's URL

https://help.github.com/articles/changing-a-remote-s-url/
```
git remote set-url https://github.com/<your_account>/ahcg_pipeline.git
```
    Check the remote's URL with : 

```
git remote -v
```
- Set the directory will not pull to the Github
 
 add directory name into hidden file .gitignore (in folder ahcg_pipeline)

### Set up environment
- Pull the pipeline (Also pull the tools if needed)

https://github.com/shashidhar22/ahcg_pipeline

- Fasta index using Samtools

  Install Samtools

  use samtools faidx to generate index:
```
samtools faidx hg19.fa
```

- Genome dict file using picard 

  Genome dict file using picard will produce hg19.dict
```
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict
```
- tips: make sure java version is higher than 1.8.

### Running pipeline

- Find help with:
```
python3 ahcg_pipeline.py -h
```
- Running code example:
```
python3 ahcg_pipeline.py 
-t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar 
-b ./lib/bowtie2-2.2.9/bowtie2 
-p ./lib/picard.jar 
-g ./lib/GenomeAnalysisTK.jar 
-i ./resources/test/test_r1.fastq ./resources/test/test_r2.fastq 
-w ./resources/genome/hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf 
-r ./resources/genome/hg19.fa 
-a ./lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa 
-o ./hw1
```

## Extract sequences for the gene of interest: BRCA1

### Download gene coordinates file for hg19
```
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
```

### Grep NM_007294 exome
```
grep NM_007294 ./resources/hg19_refGene.txt
```

### Write a code to generate bed file for NM_007294
```
./get_bed.pl > onelineBRCA1.bed
```

### Extract sequence of NM_007294
```
bedtools getfasta -s -fi ./resources/genome/hg19.fa -bed ./hw2/onelineBRCA1.bed -fo ./hw2/NM_007294.out
```

## Extracting reads mapping to BRCA1 from NA12878 HiSeq Exome dataset

### Download gene coordinates file for hg19
```
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bam
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.bam
```

### Merge BAM files
- -R option means only merge the read of on chr17
```
samtools merge -R chr17  NA12878_Chr17.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_1_NA12878.bwa.markDuplicates.bam project.NIST_NIST7086_H7AP8ADXX_CGTACTAG_2_NA12878.bwa.markDuplicates.bam
```

### Extract region of interest from whole bam file
```
samtools view NA12878_Chr17.bam -L onelineBRCA1.bed -b -o BRCA1.bam
```

### Convert BAM file to FASTQ file
```
bedtools bamtofastq -i BRCA1.bam -fq <Fastq_Read1> -fq1 <Fastq_Read2>
```

### Run the ahcg_pipeline.py
```
python3 ahcg_pipeline.py 
-t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar 
-b ./lib/bowtie2-2.2.9/bowtie2 
-p ./lib/picard.jar 
-g ./lib/GenomeAnalysisTK.jar 
-i ./resources/test/test_r1.fastq ./resources/test/test_r2.fastq 
-w ./resources/genome/hg19 
-d ./resources/dbsnp/dbsnp_138.hg19.vcf 
-r ./resources/genome/hg19.fa 
-a ./lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa 
-o ./hw3
```

## Compare the variant calling with published data

### Download the variants called on the NA12878 exome samples
```
wget http://vannberg.biology.gatech.edu/data/variants.vcf
```

### Download the verified variants calling
```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz
```

### Download the published exom data
```
wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed
```

### Preprocess vcf format
- use Preprocess.pl to add “chr” before the chromosome number in verified variants calling downloaded from GIAB
```
./Preprocess.pl > NA12878.vcf
```

### Find the intersection between the three files
```
bedtools intersect -wa -header -a variants.vcf 
-b nexterarapidcapture_exome_targetedregions_v1.2.bed 
> exom_in_variant-1.vcf

bedtools intersect -wa -header -a NA12878.vcf  
-b nexterarapidcapture_exome_targetedregions_v1.2.bed 
> exom_in_variant-2.vcf

bedtools intersect -wa -header -a exom_in_variant-1.vcf 
-b exom_in_variant-2.vcf 
> overlap_variant.vcf
``` 
- tips: the file  exom_in_variant.vcf have to have header
- tips: do three-step intersect(showed here) if your system storage is not large enough 

## Process a list of cancer gene

### Copy a list of interest gene in to a file
- generate gen_list.txt by merge the gene in the following to reference
	Reference : http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf
	Reference : http://www.genedx.com/wp-content/uploads/2016/07/Breast-Ovarian-Panel-Fact-Sheet.pdf

### Get exom in bed file
- run grep_gen_list.pl in folder get_exon_list will generate a bed file with interest exom (exon_list.bed) and the original data of the interest genes extracted from hg19_refGene.txt.
- 20 bases were add to both of the start and end.
```
./grep_gen_list.pl > exom_list.bed
```

### Variant call verification
- find the difference between our variant calling and verified variant calling on GIAB in the exom area

```
bedtools intersect -wa -header -a variants.vcf 
-b exom_list.bed 
> exom_variant-1.vcf

bedtools intersect -wa -header -a NA12878.vcf  
-b exom_list.bed 
> exom_variant-2.vcf

bedtools intersect -wa -header -a exom_in_variant-1.vcf 
-b exom_in_variant-2.vcf 
> overlap_variant.vcf
```

- Overlap results were summarized in number\_of\_overlapping_calls.txt .


## Verify Variant Calling
### Download variant calling using the fastq files for NA12878

```
wget http://vannberglab.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf
```

### Do the overlap between our course calling and reported calling

```
bedtools intersect -wa -header -a NA12878_variants.vcf 
-b exom_list.bed 
> exom_variant_course.vcf

bedtools intersect -wa -header -a NA12878.vcf  
-b exom_list.bed 
> exom_variant_repoert.vcf

bedtools intersect -wa -header -a exom_in_variant_course.vcf 
-b exom_in_variant_repoert.vcf 
> overlap.vcf
```

## VariantRecalibrator

### index and dict for gemome.ga
```
samtools faidx ../resources/vrlb/genome.fa

java -jar ../lib/picard.jar CreateSequenceDictionary REFERENCE=../resources/vrlb/genome.fa OUTPUT=../resources/vrlb/genome.dict
```


### Run VariantRecalibrator

```
java -Xmx4g -jar ../lib/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ../resources/vrlb/hg19.fa \
-input NA12878_variants_ahcg.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../resources/vrlb/hapmap_3.3.hg19.sites.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 ../resources/vrlb/1000G_omni2.5.hg19.sites.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ../resources/vrlb/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../resources/vrlb/dbsnp_138.hg19.vcf -an DP \
-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile recalibrate_SNP.recal \
-tranchesFile recalibrate_SNP.tranches

```

### Run GATK again to get vcf file
```
java -jar ../lib/GenomeAnalysisTK.jar -T ApplyRecalibration \ 
-R ../resources/vrlb/hg19.fa \ 
-input NA12878_variants_ahcg.vcf \ 
-mode SNP \ -ts_filter_level 99.0 \ -recalFile recalibrate_SNP.recal \ 
-tranchesFile recalibrate_SNP.tranches \ -o recalibrated_snps_raw_indels.vcf
```

## Calculate Read Depth Based on Alignment File

### Extract BRCA1 gene chromosome coordinates
- Extract BRCA1 gene chromosome coordinates from "exom_list.bed"
```
grep 'NM_007294' exom_list.bed > brca1.bed
```

### Extract brca1 alignments
```
samtools view -L brca1.bed project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > na12878.brca1.bam
```
-L: only output alignments overlapping in the input bed file
-b: output alignments in the bam format

### Computes and summarize coverage for brca1
```
bedtools genomecov -ibam na12878.brca1.bam -bga > na12878.brca1.bga.bed
```
-ibam BAM file as input for coverage.
-bga Reporting genome coverage for all positions in BEDGRAPH format.

- This is BedGraph format: chrom chromStart chromEnd dataValue

### Intersection between two bed files
```
bedtools intersect -split -a na12878.brca1.bga.bed -b brca1.bed -bed > brca1.final.bed
```
-split : only the exon overlaps are reported

## How to annotate the vcf file with pathogenicity

### Download the BRCA variant pathogenic annotation file
```
wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv
wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv
```

### Use the PathgenicAnnotation.py script to match the variants between vcf file and pathogenic annotation file



## Final Project

### 1.Download bam file of patient1
```
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai
```

### 2.Perform variant call analysis using GATK-HaplotyeCaller
```
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o p1_raw_variants.vcf
```

### 3.Perform variants recalibration analysis using GATK-VariantRecalibrator
```
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files
lib/jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input p1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile p1_recalibrate_SNP.recal -tranchesFile p1_recalibrate_SNP.tranches -o p1_recalibrated_snps_raw_indels.vcf
```

### 4.Create dcm-related 6 genes' genome position bed file
```
scripts/grep_gen_list.pl
```

### 5.Extract the dcm-related 6 genes' genomic variants
```
bedtools intersect -wa -header -a p1_recalibrated_snps_raw_indels.vcf -b scripts/temp/exom_list.bed > patient1_dcm_final.vcf
```

### 6.Calculat the reads depth information for DCM genes
#### 6.1.Get depth for each gene
```
samtools view -L scripts/temp/exom_list.bed Patient1_RG_MD_IR_BQ.bam -b > patient1.dcm.bam
```

#### 6.2.Reporting genome coverage
```
bedtools genomecov -ibam patient1.dcm.bam -bga > patient1.dcm.bed
```

#### 6.3.Read coverage for each gene
```
bedtools intersect -split -a patient1.dcm.bed -b scripts/temp/exom_list.bed -bed > scripts/temp/patient1.final.bed
python scripts/depth_for_each_site.py
scripts/label.pl scripts/temp/patient1.final.full.txt patient1.final.full.labeled.txt
```

### 7.Coverage plotting
```
cp patient1.final.full.labeled.txt scripts/txt/p1.final.txt
R CMD BATCH '--args scripts/txt/p1.final.txt' scripts/plots.R
```
- R 3.14 with package ggplot2 needed
- Cutoff for each exom was created in cutoff.final.csv



### 8.Vcf analysis reporting
#### 8.1.SNP varification by depth
```
cp patient1_dcm_final.vcf scripts/vcf/P1_afterinter.vcf
python scripts/min.py
```
- Depth - Cutoff > 5 : Mutation confirmed
- Depth - Cutoff between 0 ~ 5 : Ambiguous sequencing result
- Depth - Cutoff < 0 : Mutation untrusted

#### 8.2.Potential Pathogenic SNP repoert
```
cp scripts/vcf/*new.vcf SNP_report/
perl scripts/fliter.pl
```
- Collum 5 stands for the existence of the SNP in control.
- Collum 5 stands for the existence of the SNP in patients.

(Mutation confirmed += 1; Ambiguous sequencing result += 0.6; Mutation untrusted += 0; Not exist in control = -1)


