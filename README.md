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
    tips: make sure java version is higher than 1.8.

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

## Compare the variant calling with published data

