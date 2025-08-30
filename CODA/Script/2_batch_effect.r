R.version.string
.libPaths()
library(argparse)
start_time <- Sys.time()

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
library(tibble)
library(ggplot2)
set.seed(100)

# 데이터 불러오기 =========
merged_concat_df = readRDS(paste0(file_dir, "results_concat_df.Rds"))
metadata = readRDS(paste0(work_dir, "results_metadata_concat_df.Rds"))
#merged_concat_df = read.delim(paste0(file_dir, "results_count_matrix.txt"))

# 배치보정 ===========
library(edgeR)
library(limma)

# 샘플 그룹 정보 설정
sample_group_info =  as.vector(t(metadata)[, 1])
print(sprintf("sample_group_info: %s", sample_group_info))
# sample_group_info = c(rep("counts_df", 40), rep("GSE167186", 76))

merged_concat_df <- merged_concat_df %>%
  group_by(Gene_sylbol) %>%
  summarise(across(everything(), ~ round(mean(., na.rm = TRUE))))
merged_concat_df <- merged_concat_df %>% column_to_rownames(var = "Gene_sylbol") 

# DGEList 객체 생성
y <- DGEList(counts=merged_concat_df,
             gene=rownames(merged_concat_df),
             group=sample_group_info)
dim(y)
print(sprintf("Created DEGList object"))

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
write.table(Cpm.count, paste0(work_dir, "results_cpm_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
print("Calculated CPM")

# 배치효과 제거
logcpm_corrected <- removeBatchEffect(Cpm.count, batch=sample_group_info)
write.table(logcpm_corrected, paste0(work_dir, "results_cpm_matrix_corrected.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
print("Batch effect corrected")

# PCA  ==========
library(ggsci)

# PCA 시각화 (1)
pca_result <- prcomp(t(Cpm.count), scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x)

p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[1,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() + 
  labs(title = "PCA of counts_df",
       x = "Principal Component 1",
       y = "Principal Component 2")
ggsave(paste0(work_dir, "results_raw_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)
print("Visualization of PCA (1)")

# PCA 시각화 (2): batch effect corrected table
# logcpm_corrected는 counts_df 데이터를 기반으로 계산되었다고 가정합니다.
pca_result <- prcomp(t(logcpm_corrected), scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x)

p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[1,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() +
  labs(title = "PCA of counts_df (Corrected)",
       x = "Principal Component 1",
       y = "Principal Component 2")
ggsave(paste0(work_dir, "results_raw_corrected_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)
print("Visualilzation of PCA (2)")

end_time <- Sys.time()
print(end_time - start_time)
