#!/bin/bash

Ref_DIR=$1
OUTPUT_DIR=$2

STAR \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir ${OUTPUT_Dir}/STAR_index \
    --genomeFastaFiles ${Ref_DIR}/hg38.fa \
    --sjdbGTFfile ${Ref_DIR}/gencode.v47.annotation.gtf \
    --sjdbOverhang 100
