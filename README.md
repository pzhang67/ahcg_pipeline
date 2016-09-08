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
## How to run it
## Set up github
-Changing a remote's URL

https://help.github.com/articles/changing-a-remote-s-url/
```
git remote set-url https://github.com/<your_account>/ahcg_pipeline.git
```
Check the remote's URL with : 
```
git remote -v
```
-Set the directory will not pull to the Github
 
add directory name into hidden file .gitignore (in folder ahcg_pipeline)

## Set up environment
-Pull the pipeline (Also pull the tools if needed)

https://github.com/shashidhar22/ahcg_pipeline

-Fasta index using Samtools

Install Samtools

use samtools faidx to generate index:
```
samtools faidx hg19.fa
```

-Genome dict file using picard 

Genome dict file using picard will produce hg19.dict
```
java -jar picard.jar CreateSequenceDictionary R=hg19.fa O=reference.dict
```
tips: make sure java version is higher than 1.8.

## Running pipeline

-Find help with:
```
python3 ahcg_pipeline.py -h
```
-Running code example:
```
python3 ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i ./resources/test/test_r1.fastq ./resources/test/test_r2.fastq -w ./resources/genome/hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o ./hw1
```
