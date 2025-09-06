# RNA-seq Analysis Repository

RNA-seq 데이터 분석을 위한 파이프라인, 교육 자료, 그리고 실습 코드를 포함한 통합 저장소입니다.

## 디렉토리 구조

### 📁 [pipeline/](./pipeline/)
실제 RNA-seq 분석을 위한 프로덕션 파이프라인
- **DEG**: 차등발현유전자 분석 스크립트
- **Pathway**: GO, KEGG, Reactome 등 기능 분석
- **Visualization**: Volcano plot, Heatmap, PCA 등

### 📁 [training/](./training/)
RNA-seq 분석 교육 및 학습을 위한 Jupyter Notebook 자료
- **기초 분석**: EDA, 데이터 전처리
- **배치 효과**: 배치 보정 기법
- **통합 분석**: 여러 데이터셋 통합 분석
- **DEG**: 단계별 차등발현 분석 실습
- **Visualization**: 다양한 플롯 생성 방법

### 📁 [offline-server/](./offline-server/)

보안 환경에서의 RNA-seq 전처리 환경 설정 및 DEG 파이프라인
