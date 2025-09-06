[TOC]

# RNA-seq 분석

해당 디렉토리는 CODA에서 제공하는 RNA-seq 데이터 분석을 위해 작성되었습니다. 해당 문서에는 FASTQ 전처리, DEG 분석, Pathway 분석을 위한 conda pack 파일, 실행 스크립트에 대한 설명이 포함되어 있습니다. 



## 1. 분석 준비 

CODA에서는 원시데이터(FASTQ, count matrix)를 제공하지 않습니다. 분석 파일은 FASTQ이며 행렬 데이터(count matrix)로 만든 후 DEG 분석을 진행한 결과를 반출할 수 있습니다. CODA에서 분석을 위해 아래 conda pack 파일 두 개와 human reference를 준비합니다. 



해당 파일은 `/tier4/DSC/shpark/geodeg/CODA` 경로에 있습니다. 

```bash
/tier4/DSC/shpark/geodeg/CODA
├── geodeg_prep.tar.gz     # FASTQ 전처리를 위한 conda pack 파일
├── geodeg.tar.gz          # GEO 분석을 위한 conda pack 파일
└── reference
    ├── gencode.v47.annotation.gtf.gz
    ├── hg38.fa.gz
    └── RSEM_index
        └── ...
```

**Reference 디렉토리 구조**

Reference 디렉토리는 아래와 같이 구성되어야 합니다:

```
Reference 
├── gencode.v47.annotation.gtf.gz
├── hg38.fa.gz
├── RSEM_index
│   ├── chrLength.txt
│   ├── chrName.txt
│   ├── chrNameLength.txt
│   ├── chrStart.txt
│   ├── exonGeTrInfo.tab
│   ├── exonInfo.tab
│   ├── geneInfo.tab
│   ├── Genome
│   ├── genomeParameters.txt
│   ├── RSEM_index.chrlist
│   ├── RSEM_index.grp
│   ├── RSEM_index.idx.fa
│   ├── RSEM_index.n2g.idx.fa
│   ├── RSEM_index.seq
│   ├── RSEM_index.ti
│   ├── RSEM_index.transcripts.fa
│   ├── RSEM_indexLog.out
│   ├── SA
│   ├── SAindex
│   ├── sjdbInfo.txt
│   ├── sjdbList.fromGTF.out.tab
│   ├── sjdbList.out.tab
│   └── transcriptInfo.tab
└── RSEM_index.tar.gz
```



## 2. 분석 환경 세팅

아래는 CODA에서 분석하기 위한 환경을 구축하는 방법입니다. 파일은  `geodeg_prep.tar.gz`과 `geodeg_tar.gz` 이며 각각 FASTQ 전처리 및 DEG 분석을 하기위한 conda pack 파일입니다.  



### 2.1 FASTQ 전처리 분석 환경 구축

```bash
# (1) anaconda3 또는 miniconda의 envs 경로 확인 및 이동
$ which conda 
~/miniconda3/bin/conda # [conda env path]
$ cd [conda env path]/envs

# (2) geodeg_prep라는 디렉토리 생성 후, geodeg_prep.tar.gz 파일 압축 해제
$ mkdir -p [conda env path]/geodeg
$ tar -xzf geodeg_prep.tar.gz -C [conda env path]/geodeg_prep

# (3) conda 경로 정리하기
$ [conda env path]/envs/geodeg/bin/conda-unpack

# (4) conda 환경 activate
$ conda activate geodeg_prep

# (5) 활성화된 환경에서 R 실행 확인
(geodeg_prep) $ R
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-conda-linux-gnu (64-bit)
...
> q()
Save workspace image? [y/n/c]: n
```



### 2.2 DEG 분석 환경 구축

```shell
# (1) anaconda3 또는 miniconda의 envs 경로 확인 및 이동
$ which conda 
~/miniconda3/bin/conda # [conda env path]
$ cd [conda env path]/envs

# (2) geodeg라는 디렉토리 생성 후, geodeg.tar.gz 파일 압축 해제
$ mkdir -p [conda env path]/geodeg
$ tar -xzf geodeg.tar.gz -C [conda env path]/geodeg

# (3) conda 경로 정리하기
$ [conda env path]/envs/geodeg/bin/conda-unpack

# (4) conda 환경 activate
$ conda activate geodeg

# (5) 활성화된 환경에서 R 실행 확인
(geodeg) $ R
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-conda-linux-gnu (64-bit)
...
> q()
Save workspace image? [y/n/c]: n
```



## 3. 분석 스크립트

환경 세팅이 다 되었다면 분석을 시작합니다. `Run` 디렉토리의 `shell` 스크립트를 구동합니다. 



### 3.1 FASTQ 전처리 

앞서 준비한 `geodeg_prep` Conda 환경을 활성화한 후, `WTS_pipeline_prep.sh`를 실행하여 파이프라인을 구동합니다. 필요에 따라 `WTX_pipeline_index.sh`를 실행하여 RSEM 또는 STAR 인덱스를 생성합니다. 

```bash
dsc/DSC-human/RNA-KDCA/CODA
├── Run
│   ├── WTS_pipeline_index.sh
│   ├── WTS_pipeline_prep.sh
│   ├── WTS_pipeline_prep_RSEM_only.sh
│   └── ..
└── ...
```



### 3.2 DEG 분석 

앞서 준비한 `geodeg` Conda 환경을 활성화한 후, DEG 분석을 위해 `WTS_pipeline_DEG.sh` 스크립트를 실행합니다. 

```bash
dsc/DSC-human/RNA-KDCA/CODA
├── Run
│   ├── WTS_pipeline_DEG.sh
│   ├── 1_gene_mapping.sh
│   ├── 2_batch_effect.sh
│   ├── 3_gsva.sh
│   ├── 4_deg_analysis.sh
│   └── ...
└── Script
    ├── 1_gene_mapping.r
    ├── 2_batch_effect.r
    ├── 3_gsva.r
    ├── 4_deg_analysis.r
    └── ...
```



### 3.3 시각화 

다음은 DEG 결과를 시각화하기 위한 스크립트입니다. 아래 스크립트에서 Volcano plot, Heatmap, Vendiagram을 생성할 수 있습니다.  

```bash
dsc/DSC-human/RNA-KDCA/CODA
├── Run
│   └── ...
└── Script
    ├── visualization_muscle.R
    ├── visualization_pancreas.R
    └── ...
```



### 3.4 Pathway 분석

다음은 Pathway 분석을 위한 스크립트입니다. DEG 분석 결과 파일이 필요하며, FDR, PValue, LogFC 값을 한 후, DEG에 대한 GO, KEGG, Reactome, WikiPathway, MsigDB 기능 분석을 수행합니다. 분석 결과는 각각의 엑셀 시트로 저장됩니다. 

```bash
dsc/DSC-human/RNA-KDCA/CODA
├── Run
│   └── ...
└── Script
    ├── pathway_analysis.R
    └── ...
```



## 4. 참고

아래는 위 환경과 스크립트가 정상적으로 실행되는지 확인하기 위해 테스트한 환경입니다. 

|      | 항목           | FASTQ 전처리                                                 | DEG 분석                                       |
| ---- | -------------- | ------------------------------------------------------------ | ---------------------------------------------- |
| 1    | 로컬 분석 환경 | shpark@incogpu                                               | user@incogpu                                   |
| 2    | 스크립트       | WTS_pipeline_prep.sh                                         | WTS_pipeline_DEG.sh                            |
| 3    | 패키지 경로    | **STAR**: /home/shpark/anaconda3/envs/geodeg/bin/STAR<br />**rsem**: /home/shpark/anaconda3/envs/geodeg/bin/rsem-calculate-expression | /data/user/anaconda3/envs/geodeg/lib/R/library |
| 4    | 전달 파일      | geodeg_prep.tar.gz (360M)                                    | geodeg.tar.gz (1.2G)                           |
