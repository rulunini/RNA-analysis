R.version.string
.libPaths()
library(argparse)

# 인수 파서 설정
parser <- ArgumentParser()
parser$add_argument("-i", "--input_directory", help="The results directory", required=TRUE)
parser$add_argument("-o", "--output_directory", help="The results directory", required=TRUE)

args <- parser$parse_args()

# 인수를 사용하여 경로 설정
file_dir <- args$input_directory
work_dir <- args$output_directory

# 라이브러리 로드 =========

library(dplyr)
library(ggplot2)
library("tidyverse")
library(biomaRt)
library(org.Hs.eg.db)
library('msigdbr')
library(tibble)
library(cowplot)
library(xlsx)
library(edgeR)
set.seed(100)


# 데이터 불러오기 =========
merged_concat_df = readRDS(paste0(file_dir, "results_concat_df.Rds"))
metadata = readRDS(paste0(work_dir, "results_metadata_concat_df.Rds"))
merged_concat_df = read.delim(paste0(file_dir, "results_count_matrix.txt"))

merged_concat_df <- merged_concat_df %>%
  group_by(Gene_sylbol) %>%
  summarise(across(everything(), ~ round(mean(., na.rm = TRUE))))
merged_concat_df <- merged_concat_df %>% column_to_rownames(var = "Gene_sylbol")

# DEG 분석 ========
library(edgeR)
sample_group_info = as.vector(t(metadata)[, 1])
batch <- factor(sample_group_info) 

# 샘플 타입 선정
Sample_type <- as.character(metadata[1, ])
Sample_type[Sample_type %in% c("group: Sarcopenia", "sarcopenia status: yes")] = "Sarcopenia"
Sample_type[Sample_type %in% c("group: Old Healthy", "sarcopenia status: no")] = "Healthy"
sprintf("Sample_type %s", Sample_type) # Sarcopenia, Healthy, Sarcopenia, ....

# (2) DGEList 객체 생성 및 정규화
dge <- DGEList(counts=merged_concat_df,
               gene=rownames(merged_concat_df),
               group=sample_group_info)

dge$samples$group <- factor(Sample_type)
dge$samples$batch <- batch
sprintf("dge group %s", dge$samples$group) # Healthy, Healthy, ....
sprintf("dge batch %s", dge$samples$batch) # sarcopenia status: no, sarcopenia, status: no, ....

# (3) 데이터 정규화
dge <- calcNormFactors(dge)

# (4) 디자인 매트릭스 생성 (배치 + 처리)
design <- model.matrix(~ group, data=dge$samples)
# design <- model.matrix(~ batch + group, data=dge$samples)
print(dge$samples)
print("Made design matrix")
print(design)
print("Printed design")

# (5) 디스퍼전 추정 
dge <- estimateDisp(dge, design)
print("Estimated Dispersion")

# (6) GLM 적합
fit <- glmFit(dge, design)

# (7) 그룹 간 비교를 위한 대비 설정 
#  groupA - groupB
print(colnames(design))
contrast <- makeContrasts(groupSarcopenia, levels=colnames(design))
# contrast <- makeContrasts(groupSarcopenia-groupHealthy, levels=colnames(design))

# (8) LRT 테스트 실행
lrt <- glmLRT(fit, contrast=contrast)

# (9) 결과 확인
topTags(lrt)
result_table <- topTags(lrt, n=Inf)  # 'n=Inf' 를 해서 전체 목록 도출
print("Analyzed DEG")

# (10) 볼케이노 플롯 생성
library(EnhancedVolcano)
p=EnhancedVolcano(result_table$table,
                  lab = rownames(result_table$table),
                  x = 'logFC',  # log fold change
                  y = 'PValue',  # p-value
                  xlim = c(-log2(2), log2(2)),
                  title = 'Volcano Plot',
                  pCutoff = 0.05,  # P-value cutoff
                  FCcutoff = 0.5,  # Fold change cutoff
                  pointSize = 2.0,
                  labSize = 3.0, ylim = c(0, 10))
ggsave(paste0(work_dir, "results_Enhanced_volanoPlot.png"), plot = p, width = 10, height = 9)
