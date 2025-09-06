.libPaths("/BiO/users/jjg/tools/anaconda3_R_4.3.2_LibPaths/")
library(tibble)
library(dplyr)
library(cowplot)
library(ggplot2)
library('xlsx')
library("tidyverse")
library(org.Hs.eg.db)

set.seed(100)


work_dir = "/BiO/users/jjg/ISG/Aging_GSE/240503_Use_CPM_matrix/"


GSE111016 = read.delim("/BiO/users/jjg/ISG/Aging_GSE/GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv", sep = ',')
GSE111016_sheet = read.delim("/BiO/users/jjg/ISG/Aging_GSE/GSE111016_series_matrix.txt", sep ='\t', head = FALSE)
GSE111016_sheet = GSE111016_sheet[c(1, 11), seq(2, ncol(GSE111016_sheet))]

GSE167186 = read.delim("/BiO/users/jjg/ISG/Aging_GSE/GSE167186_counts.csv", sep = ',')
GSE167186_sheet = read.delim("/BiO/users/jjg/ISG/Aging_GSE/GSE167186-GPL20301_series_matrix.txt", sep ='\t', head = FALSE)
GSE167186_sheet = GSE167186_sheet[c(1, 10), seq(2, ncol(GSE167186_sheet))]
#GSE167186_sheet = GSE167186_sheet[,order(as.numeric(gsub("X_", "", GSE167186_sheet[1,])))]

identical(as.character(GSE167186_sheet[1,]), colnames(GSE167186)[seq(2, ncol(GSE167186))])
library(edgeR)
library(limma)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- GSE111016[,1]  # 예시 Gene ID
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)


matched_gene_list = sapply(GSE111016[,1], function(x) {gene_info[gene_info[,1] %in% x,2]})

matched_gene_list_length = sapply(GSE111016[,1], function(x) {nrow(gene_info[gene_info[,1] %in% x,2])})
sum(matched_gene_list_length %>% unlist())

GSE111016$Gene_sylbol = as.character(matched_gene_list)
GSE111016 = GSE111016[!is.na(GSE111016$Gene_sylbol),]
GSE111016 = GSE111016[GSE111016$Gene_sylbol != "", ]

GSE111016$X = NULL
GSE167186$Gene_sylbol = GSE167186$Symbol
GSE167186$Symbol = NULL

merged_df = merge(GSE111016, GSE167186, by = "Gene_sylbol")

merged_concat_df = sapply(unique(merged_df$Gene_sylbol), function(x) {colSums(merged_df[ merged_df$Gene_sylbol %in% c(x), colnames(merged_df)!="Gene_sylbol"] )}) %>% data.frame()

merged_concat_df = t(merged_concat_df)
rownames(merged_concat_df) = unique(merged_df$Gene_sylbol)
colnames(merged_concat_df) = colnames(merged_df)[colnames(merged_df)!="Gene_sylbol"]

#saveRDS(merged_concat_df, paste0(work_dir, "merged_concat_df.Rds"))
#metadata = cbind(GSE111016_sheet, GSE167186_sheet)
#saveRDS(metadata, paste0(work_dir, "metadata_concat_df.Rds"))


merged_concat_df = readRDS(paste0(work_dir, "merged_concat_df.Rds"))
metadata = readRDS(paste0(work_dir, "metadata_concat_df.Rds"))

write.table(merged_concat_df, paste0(work_dir, "merged_count_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)



sample_group_info = c(rep("GSE111016", 40), rep("GSE167186", 76))

y <- DGEList(counts=merged_concat_df,
             gene=rownames(merged_concat_df),
             group=sample_group_info)
dim(y)
#keep <- filterByExpr(y, group=sample_group_info)
y <- y[,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
#y$samples
eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
Cpm.count <- edgeR::cpm(y, log=TRUE)
write.table(Cpm.count, paste0(work_dir, "merged_cpm_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
logcpm_corrected <- removeBatchEffect(Cpm.count, batch=sample_group_info)
write.table(logcpm_corrected, paste0(work_dir, "merged_cpm_matrix_corrected.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)

pca_result <- prcomp(t(Cpm.count), scale. = FALSE)


pca_scores <- as.data.frame(pca_result$x)

p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() + 
  labs(title = "PCA of GSE111016, GSE167186",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave(paste0(work_dir, "raw_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)


pca_result <- prcomp(t(logcpm_corrected), scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x)

p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() +
  labs(title = "PCA of GSE111016, GSE167186",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave(paste0(work_dir, "raw_correctedCPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)



# 카운트 데이터와 그룹 정보, 배치 정보

batch <- factor(sample_group_info)  # 배치 정보
Sample_type = as.character(metadata[2,])
Sample_type[Sample_type %in% c("group: Sarcopenia", "sarcopenia status: yes")] = "Sarcopenia"
Sample_type[Sample_type %in% c("group: Old Healthy", "sarcopenia status: no")] = "Healthy"
Sample_type[Sample_type %in% c("group: Young Healthy")] = "Young_Healthy"

# DGEList 객체 생성 및 샘플 정보 입력

dge <- DGEList(counts=merged_concat_df,
             gene=rownames(merged_concat_df))
dge$samples$group <- factor(Sample_type)
dge$samples$batch <- batch

# 데이터 정규화
dge <- calcNormFactors(dge)

# 디자인 매트릭스 생성 (배치 + 처리)
design <- model.matrix(~ batch + group, data=dge$samples)

# 디스퍼전 추정
dge <- estimateDisp(dge, design)

# GLM 적합
fit <- glmFit(dge, design)

# 이 경우 groupC는 기준이므로 계수 비교는 groupA - groupB
contrast <- makeContrasts(groupSarcopenia-groupHealthy, levels=design)

# LRT 테스트 실행
lrt <- glmLRT(fit, contrast=contrast)
# 결과 확인
topTags(lrt)



result_table <- topTags(lrt, n=Inf)  # 'n=Inf'는 모든 유전자 결과를 가져옴

# 볼케이노 플롯 생성
p=EnhancedVolcano(result_table$table,
                lab = rownames(result_table$table),
                x = 'logFC',  # log fold change
                y = 'PValue',  # p-value
                xlim = c(-log2(2), log2(2)),  # x축 범위 설정
                title = 'Volcano Plot',
                pCutoff = 0.05,  # P-value cutoff
                FCcutoff = 0.5,  # Fold change cutoff
                pointSize = 2.0,
                labSize = 3.0, ylim = c(0, 10))
ggsave(paste0(work_dir, "Enhanced_volanoPlot.png"), plot = p, width = 10, height = 9)




