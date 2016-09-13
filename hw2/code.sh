# get BRCA1.bam
samtools view -L onelineBRCA1.bed -b -o BRCA1.bam ../resources/NA12878/bam/NA12878_Chr17.bam

# get fasta
bedtools bamtofastq -i BRCA1.bam -fq ../resources/NA12878/NA12878_brca_r1.fastq -fq2 ../resources/NA12878/NA12878_brca_r2.fastq