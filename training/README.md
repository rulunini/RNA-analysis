# RNA-seq Analysis Training Materials

RNA-seq 데이터 분석 교육을 위한 Jupyter Notebook 실습 자료입니다. 단계별 학습을 통해 전사체 분석의 전 과정을 익힐 수 있습니다.

## 📚 학습 단계

### 1단계: 기초 분석
- **`1-1_EDA.ipynb`**: 탐색적 데이터 분석 (EDA) 기초
- **`1-2_DEG_GSE167186.ipynb`**: GSE167186 데이터셋 DEG 분석
- **`1-3_DEG_GSE111016.ipynb`**: GSE111016 데이터셋 DEG 분석

### 2단계: 고급 분석 및 시각화
- **`2-1_EDA.ipynb`**: 심화 탐색적 데이터 분석
- **`2-2_Visualization.ipynb`**: 기본 시각화 기법
- **`2-2_DEG_Combined.ipynb`**: 다중 데이터셋 DEG 분석
- **`2-3_Batchcorrection.ipynb`**: 배치 효과 보정
- **`2-4_Visualization.ipynb`**: 고급 시각화
- **`2-5_Visualization2.ipynb`**: 추가 시각화 기법

### 3단계: 통합 분석
- **`3-1_All.ipynb`**: 전체 파이프라인 통합 실습

## 📊 실습 데이터

### 포함된 데이터셋
- **GSE111016**: 
  - `GSE111016_counts_matrix.csv`: Count matrix 데이터
  - `GSE111016_metadata_matrix.csv`: 메타데이터
  - `GSE111016.ipynb`: 전체 분석 워크플로우

- **GSE167186**: 
  - `GSE167186.ipynb`: 전체 분석 워크플로우

## 🎯 학습 목표

### 기초 단계
- RNA-seq 데이터 구조 이해
- 기본적인 EDA 수행
- 단일 데이터셋 DEG 분석

### 중급 단계
- 다양한 시각화 기법 습득
- 배치 효과 이해 및 보정
- 다중 데이터셋 비교 분석

### 고급 단계
- 전체 분석 파이프라인 구축
- 결과 해석 및 생물학적 의미 도출

## 🛠️ 필요 패키지

```python
# 데이터 처리
import pandas as pd
import numpy as np

# 생물정보학
import scanpy as sc
import anndata as ad

# 통계 분석
from scipy import stats
import statsmodels.api as sm

# 시각화
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

# DEG 분석 (R 패키지 - rpy2 필요)
# DESeq2, edgeR, limma
```

## 📖 사용 방법

### 순차적 학습 (권장)
```bash
# 1단계부터 순서대로 실행
jupyter notebook 1-1_EDA.ipynb
jupyter notebook 1-2_DEG_GSE167186.ipynb
jupyter notebook 1-3_DEG_GSE111016.ipynb

# 2단계 진행
jupyter notebook 2-1_EDA.ipynb
# ... 계속
```

### 특정 주제 학습
- **EDA만**: `1-1_EDA.ipynb`, `2-1_EDA.ipynb`
- **시각화만**: `2-2_Visualization.ipynb`, `2-4_Visualization.ipynb`, `2-5_Visualization2.ipynb`
- **DEG 분석만**: `1-2_DEG_GSE167186.ipynb`, `1-3_DEG_GSE111016.ipynb`
- **배치 보정**: `2-3_Batchcorrection.ipynb`

## 📋 주요 학습 내용

### 데이터 전처리
- Count matrix 로딩 및 검증
- 메타데이터 연결
- 품질 관리 및 필터링

### 탐색적 데이터 분석
- 기본 통계량 계산
- 분포 확인 및 이상치 탐지
- 샘플 간 상관관계 분석

### 차등발현 분석
- 정규화 방법 비교
- 통계 검정 수행
- 다중 검정 보정

### 시각화
- PCA, t-SNE, UMAP
- Volcano plot, MA plot
- Heatmap, Box plot
- 네트워크 분석 시각화

### 배치 효과 보정
- 배치 효과 탐지
- ComBat, Limma 보정
- 보정 전후 비교

## 💡 학습 팁

- 각 노트북은 독립적으로 실행 가능하지만 순서대로 학습 권장
- 코드를 직접 수정해보며 매개변수 변화에 따른 결과 관찰
- 실제 연구 데이터를 사용하여 실무 경험 축적
- 결과 해석에 중점을 두고 생물학적 의미 고민
