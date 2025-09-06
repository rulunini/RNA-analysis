R.version.string
.libPaths()

# 라이브러리 로드
library(dplyr)
library(ggplot2)
library("tidyverse") 
library(biomaRt) 
library(org.Hs.eg.db) 
library('msigdbr')
library(tibble) 
library(cowplot) 
library(xlsx)
library(argparse)
set.seed(100)

# 인수 파서 설정
parser <- ArgumentParser()
parser$add_argument("-i", "--input_directory", help="The results directory", required=TRUE)
parser$add_argument("-c", "--input_counts", help="The path of counts.csv", required=TRUE)
parser$add_argument("-m", "--input_matrix", help="The path of series matrix.txt", required=TRUE)
parser$add_argument("-o", "--output_directory", help="The results directory", required=TRUE)

args <- parser$parse_args()

# 인수를 사용하여 경로 설정
file_dir <- args$input_directory
counts_file <- args$input_counts
matrix_file <- args$input_matrix
work_dir <- args$output_directory


# 데이터 로드
counts_df = read.delim(paste0(file_dir, counts_file), sep = ',')
counts_sheet = read.delim(paste0(file_dir, matrix_file), sep = ',')
counts_sheet = read.delim(paste0(file_dir, matrix_file), sep ='\t', head = FALSE)




# (1) 데이터 전처리 ===========
library(edgeR)
library(limma)

# Ensemble에서 유전자 매핑
# Ensembl Mart 설정
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 유전자 ID 추출
gene_ids <- counts_df[,1]  # Gene ID
# Ensembl에서 유전자 정보 가져오기
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)
# 매칭된 유전자 이름 추출
matched_gene_list = sapply(counts_df[,1], function(x) {gene_info[gene_info[,1] %in% x,2]})
# 매칭된 유전자 수 계산
matched_gene_list_length = sapply(counts_df[,1], function(x) {nrow(gene_info[gene_info[,1] %in% x,2])})
# 총 매칭된 유전자 수 확인
sum(matched_gene_list_length %>% unlist())

# 데이터 프레임에 유전자 기호 추가 및 NA 값 제거
counts_df$Gene_sylbol = as.character(matched_gene_list)
counts_df = counts_df[!is.na(counts_df$Gene_sylbol),]
counts_df = counts_df[counts_df$Gene_sylbol != "", ]

# 불필요한 열 삭제
counts_df$X = NULL

# 데이터 저장
saveRDS(merged_concat_df, paste0(work_dir, "results_concat_df.Rds"))
saveRDS(metadata, paste0(work_dir, "results_metadata_concat_df.Rds"))
write.table(merged_concat_df, paste0(work_dir, "results_count_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)





# 데이터 불러오기 
merged_concat_df = readRDS(paste0(file_dir, "results_concat_df.Rds"))
metadata = readRDS(paste0(work_dir, "results_metadata_concat_df.Rds"))

merged_concat_df = read.delim(paste0(file_dir, "results_count_matrix.txt"))





# (2) 배치보정 ===========
# 샘플 그룹 정보 설정
sample_group_info = c(rep("counts_df", 40), rep("GSE167186", 76))

# DGEList 객체 생성
y <- DGEList(counts=merged_concat_df,
             gene=rownames(merged_concat_df),
             group=merged_concat_df)
dim(y)

# 유전자 필터링
#keep <- filterByExpr(y, group=sample_group_info)
y <- y[,, keep.lib.sizes=FALSE]

# 정규화
# TMM 방법으로 정규화 요인 계산
y <- calcNormFactors(y, method="TMM")
#y$samples

# 효과적인 라이브러리 크기 계산
eff.lib.size <- y$samples$lib.size*y$samples$norm.factors

# 로그 변환한 CPM 계산
Cpm.count <- edgeR::cpm(y, log=TRUE)
# write.table(Cpm.count, paste0(work_dir, "results_cpm_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)

# 배치효과 제거
logcpm_corrected <- removeBatchEffect(Cpm.count, batch=sample_group_info)
write.table(logcpm_corrected, paste0(work_dir, "results_cpm_matrix_corrected.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)





# (3) PCA 분석 ===========
# PCA 수행
pca_result <- prcomp(t(Cpm.count), scale. = FALSE)

# PCA 결과를 데이터 프레임으로 변환
pca_scores <- as.data.frame(pca_result$x)

# PCA 시각화
library(ggsci)
p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(counts_sheet[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() + 
  labs(title = "PCA of counts_df",
       x = "Principal Component 1",
       y = "Principal Component 2")
ggsave(paste0(work_dir, "results_raw_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)

# 정정된 CPM으로 PCA 수행
# logcpm_corrected는 counts_df 데이터를 기반으로 계산되었다고 가정합니다.
pca_result <- prcomp(t(logcpm_corrected), scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x)

# 정정된 PCA 시각화
p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(counts_sheet[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() +
  labs(title = "PCA of counts_df (Corrected)",
       x = "Principal Component 1",
       y = "Principal Component 2")
ggsave(paste0(work_dir, "results_raw_corrected_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)





# (4) GSVA분석 ===========
# 데이터 수집 민 준비
library(GSVA)
msigdb_data <- msigdbr(species = "Homo sapiens")#, category = "H") # "H"는 hallmark gene sets
pathways <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)

expr_matrix <- as.matrix(logcpm_corrected)

# 결과 탐색 및 요약
ssgsea_results <- gsva(expr_matrix, pathways, method = "ssgsea")

# ssGSEA 결과 확인
head(ssgsea_results)

ssgsea_df = data.frame(t(ssgsea_results), "Annotation"=as.character(metadata[2,]))
table(as.character(metadata[2,]))
#group: Old Healthy      group: Sarcopenia    group: UNCLASSIFIED   group: Young Healthy  sarcopenia status: no sarcopenia status: yes 
#29                     24                      4                     19                     20                     20 
head(ssgsea_df[seq(1,10), seq(1,10)])
write.xlsx(ssgsea_df, paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "ssGSEA", append=FALSE)

# 시각화 및 통계 분석
pvalue = apply(ssgsea_df[, seq(1, ncol(ssgsea_df)-1)], 2,
      function(x) {wilcox.test(x[ssgsea_df$Annotation %in% c("group: Sarcopenia", "sarcopenia status: yes")],
                               x[ssgsea_df$Annotation %in% c("group: Old Healthy")])$p.value})
pvalue = pvalue[order(pvalue)]

p = ggplot(ssgsea_df, aes(x=Annotation, y= HALLMARK_PROTEIN_SECRETION, fill=Annotation)) + geom_boxplot()
p = p + theme_bw()
write.xlsx(ssgsea_df, paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "ssGSEA")





# (5) DEG 분석 ========
# (1) DEG 분석을 위한 데이터 준비
batch <- factor(sample_group_info) 
Sample_type = as.character(metadata[2,])

# 샘플 타입 설정
Sample_type[Sample_type %in% c("group: Sarcopenia", "sarcopenia status: yes")] = "Sarcopenia"
Sample_type[Sample_type %in% c("group: Old Healthy", "sarcopenia status: no")] = "Healthy"
Sample_type[Sample_type %in% c("group: Young Healthy")] = "Young_Healthy"

# (2) DGEList 객체 생성 및 정규화
dge <- DGEList(counts=merged_concat_df,
               gene=rownames(merged_concat_df))
dge$samples$group <- factor(Sample_type)
dge$samples$batch <- batch

# (3) 데이터 정규화
dge <- calcNormFactors(dge)

# (4) 디자인 매트릭스 생성 (배치 + 처리)
design <- model.matrix(~ batch + group, data=dge$samples)

# (5) 디스퍼전 추정
dge <- estimateDisp(dge, design)

# (6) GLM 적합
fit <- glmFit(dge, design)

# (7) 그룹 간 비교를 위한 대비 설정
#  groupA - groupB
contrast <- makeContrasts(groupSarcopenia-groupHealthy, levels=design)

# (8) LRT 테스트 실행
lrt <- glmLRT(fit, contrast=contrast)

# (9) 결과 확인
topTags(lrt)
result_table <- topTags(lrt, n=Inf)  # 'n=Inf' 를 해서 전체 목록 도출

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





# (6) GSEA 분석 ===========
library(msigdbr)

deg_results <- topTags(lrt, n = Inf)$table

# 3. Prepare Gene List for GSEA
gene_list <- deg_results$logFC
names(gene_list) <- rownames(deg_results)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA 분석 수행

# 4. Load MSigDB Gene Sets
# 예: hallmark gene sets 사용
msigdb_data <- msigdbr(species = "Homo sapiens")#, category = "H") # "H"는 hallmark gene sets
pathways <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)

# 5. GSEA Analysis with fgsea
gsea_results <- fgsea(pathways = pathways,
                      stats = gene_list,
                      minSize = 10,
                      maxSize = 500,
                      nperm = 1000)
head(gsea_results)
dim(gsea_results)

# Heatmap 시각화 및 결과 저장
library(ComplexHeatmap)

z_score = apply(logcpm_corrected, 1, function(x) {(x-mean(x))/sd(x)}) %>% data.frame() %>% t() %>% data.frame()

Heatmap(z_score[deg_results$genes[seq(1,50)],], name = "Corrected", 
        column_split = as.character(metadata[2,]))

write.xlsx(data.frame(gsea_results[seq(1,100),]), paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "GSEA", append = TRUE)
openxlsx::write.xlsx(gsea_results, file = paste0(work_dir, "GSEA, ssGSEA results2.xlsx"))


