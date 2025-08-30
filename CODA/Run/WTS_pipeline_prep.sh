#!/bin/bash

NAME=$1
WTS_DIR=$2
INDEX_DIR=$3
REF_DIR=$4


## data gunzip

R1_gz=$(ls $WTS_DIR/*R1.fastq.gz)
R1=${R1_gz//.gz/}

R2_gz=$(ls $WTS_DIR/*R2.fastq.gz)
R2=${R2_gz//.gz/}

gunzip -d < $R1_gz > $R1
gunzip -d < $R2_gz > $R2


## Alignment (by STAR)

# 1 pass alignment
mkdir ${WTS_DIR}/${NAME}
mkdir ${WTS_DIR}/${NAME}/1pass_${NAME}

cd ${WTS_DIR}/${NAME}/1pass_${NAME}

STAR \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM \
    --runThreadN 4  \
    --genomeDir ${INDEX_DIR} \
    --readFilesIn $R1 $R2

# Generate genome index
mkdir ${WTS_DIR}/${NAME}/1pass_genome_${NAME}
cd ${WTS_DIR}/${name}/1pass_genome_${NAME}

STAR \
    --runMode genomeGenerate \
    --sjdbOverhang 100 \
    --genomeChrBinNbits 12 \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles ${REF_DIR}/hg38.fa \
    --sjdbGTFfile ${REF_DIR}/gencode.v47.annotation.gtf \
    --genomeDir ${WTS_DIR}/${NAME}/1pass_genome_${NAME} \
    --sjdbFileChrStartEnd ${WTS_DIR}/${NAME}/1pass_${NAME}/SJ.out.tab

# 2 pass alignment
mkdir ${WTS_DIR}/${NAME}/2pass_${NAME}
cd ${WTS_DIR}/${NAME}/2pass_${NAME}

STAR \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM \
    --runThreadN 4 \
    --genomeDir ${WTS_DIR}/${NAME}/1pass_genome_${NAME} \
    --readFilesIn $R1 $R2


## Quantification (by RSEM)
cd ${WTS_DIR}/${NAME}/2pass_${NAME}

rsem-calculate-expression \
    --alignments \
    --no-bam-output \
    --paired-end ${WTS_DIR}${NAME}/2pass_${NAME}/Aligned.toTranscriptome.out.bam \
    ${INDEX_DIR}\
    Quantification
