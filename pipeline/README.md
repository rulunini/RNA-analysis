# GEO RNA-seq Analysis Pipeline

GEO 데이터베이스의 공개 데이터(public data)를 활용한 RNA-seq 차등발현유전자(DEG) 분석 및 경로 분석 파이프라인입니다.

## 데이터셋

공개 데이터베이스 GEO(Gene Expression Omnibus)에서 다운로드한 데이터셋:

- **GSE167186**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167186
- **GSE111016**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111016

## 분석 워크플로우

### 1. 차등발현유전자 분석 (`deg_analysis.r`)
- 공개 데이터 전처리
- DEG 분석 수행
- 통계적 유의성 검정 (FDR < 0.05, |logFC| > 1)

**사용법:**
```bash
Rscript deg_analysis.r -i <input_directory> -o <output_directory>
```

### 2. 경로 분석 (`pathway_analysis.R`)
- GO 분석 (Gene Ontology)
- KEGG 경로 분석
- Reactome 경로 분석
- MSigDB 분석

**주요 기능:**
- clusterProfiler를 이용한 enrichment 분석
- 다중 데이터베이스 통합 분석
- 결과 시각화 및 엑셀 파일 출력

### 3. 시각화 (`visualization.R`)
- Volcano plot 생성
- Heatmap 생성 (Z-score 변환)
- Venn diagram을 통한 데이터셋 간 비교
- ComplexHeatmap을 이용한 고급 히트맵

## 필요 패키지

```r
# 기본 분석
library(dplyr)
library(ggplot2)
library(tidyverse)

# DEG 분석
library(argparse)

# 경로 분석
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(msigdbr)

# 시각화
library(ComplexHeatmap)
library(EnhancedVolcano)
library(VennDiagram)
library(openxlsx)
```

## 출력 파일

- `*_DEG_analysis.csv`: 차등발현유전자 분석 결과
- `*_pathway_analysis.xlsx`: 경로 분석 결과 (다중 시트)
- `*.png`: 각종 시각화 결과 (volcano plot, heatmap, venn diagram)

## 분석 파라미터

- **FDR threshold**: 0.05
- **LogFC threshold**: ±1.0
- **Z-score 변환**: 히트맵 생성 시 적용

## 사용 예시

1. 전체 파이프라인 실행:
```bash
# DEG 분석
Rscript deg_analysis.r -i ./data -o ./results

# 경로 분석 (R 스크립트 내 경로 수정 필요)
Rscript pathway_analysis.R

# 시각화 (R 스크립트 내 경로 수정 필요)
Rscript visualization.R
```

2. 개별 분석 단계별 실행 가능

## 주의사항

- 공개 데이터 사용 시 해당 연구의 라이선스 및 인용 정보 확인
- `pathway_analysis.R`와 `visualization.R`의 경로 설정을 프로젝트에 맞게 수정 필요
- 대용량 데이터 처리 시 메모리 사용량 고려
- R 버전 및 패키지 호환성 확인 권장
