name=$1
WTS_Dir=$2


##data gunzip

R1_gz=$(ls $WTS_Dir/*R1.fastq.gz)
R1=${R1_gz//.gz/}

R2_gz=$(ls $WTS_Dir/*R2.fastq.gz)
R2=${R2_gz//.gz/}

gunzip -d < $R1_gz > $R1
gunzip -d < $R2_gz > $R2


##1 pass aligned
mkdir ${WTS_Dir}/${name}
mkdir ${WTS_Dir}/${name}/1pass_${name}

cd ${WTS_Dir}/${name}/1pass_${name}

STAR --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --runThreadN 4  --genomeDir /BiO/LabData/Pipeline/RefDB/FusionCaller/star --readFilesIn $R1 $R2

#1 pass_genome
mkdir ${WTS_Dir}/${name}/1pass_genome_${name}
cd ${WTS_Dir}/${name}/1pass_genome_${name}

STAR --runMode genomeGenerate --sjdbOverhang 100 --genomeChrBinNbits 12 --genomeSAindexNbases 10  --genomeFastaFiles /BiO/LabData/Pipeline/tools/Nextflow/rnaseq/Ref/genome/hg38.fa --sjdbGTFfile /BiO/LabData/Pipeline/tools/Nextflow/rnaseq/Ref/genome/gencode.v46.annotation.gtf --genomeDir ${WTS_Dir}/${name}/1pass_genome_${name} --sjdbFileChrStartEnd ${WTS_Dir}/${name}/1pass_${name}/SJ.out.tab
#

#2pass
mkdir ${WTS_Dir}/${name}/2pass_${name}
cd ${WTS_Dir}/${name}/2pass_${name}

STAR --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --runThreadN 4  --genomeDir ${WTS_Dir}/${name}/1pass_genome_${name} --readFilesIn $R1 $R2



#RSEM

cd ${WTS_Dir}/${name}/2pass_${name}

rsem-calculate-expression --alignments --no-bam-output --paired-end ${WTS_Dir}${name}/2pass_${name}/Aligned.toTranscriptome.out.bam /BiO/LabData/Pipeline/tools/Nextflow/rnaseq/Ref/genome/rsem Quant



