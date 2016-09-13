samtools faidx ../resources/chr17/chr17.fa
java -jar ../lib/picard.jar CreateSequenceDictionary R=../resources/chr17/chr17.fa O=../resources/chr17/chr17.dict
