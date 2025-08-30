.libPaths("C:/Users/ASUS/Desktop/Rlib/")
#library(tibble)
library(dplyr)
#library(cowplot)
#library(ggplot2)
library('xlsx')
library(biomaRt)
library("tidyverse")
library(org.Hs.eg.db)
library(ggplot2)
library(msigdbr)

set.seed(100)


file_dir = "C:\\Users\\ASUS\\Dropbox\\ISG\\질병연구청\\240517\\"
work_dir = paste0(file_dir, "output\\")

GSE111016 = read.delim(paste0(file_dir, "GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv"), sep = ',')
GSE111016_sheet = read.delim(paste0(file_dir,"GSE111016_series_matrix.txt"), sep ='\t', head = FALSE)
GSE111016_sheet = GSE111016_sheet[c(1, 11), seq(2, ncol(GSE111016_sheet))]

GSE167186 = read.delim(paste0(file_dir, "GSE167186_counts.csv"), sep = ',')
GSE167186_sheet = read.delim(paste0(file_dir,"GSE167186-GPL20301_series_matrix.txt"), sep ='\t', head = FALSE)
GSE167186_sheet = GSE167186_sheet[c(1, 10), seq(2, ncol(GSE167186_sheet))]
#GSE167186_sheet = GSE167186_sheet[,order(as.numeric(gsub("X_", "", GSE167186_sheet[1,])))]

identical(as.character(GSE167186_sheet[1,]), colnames(GSE167186)[seq(2, ncol(GSE167186))])

library(edgeR)
library(limma)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- GSE111016[,1]  # Gene ID
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
metadata = cbind(GSE111016_sheet, GSE167186_sheet)
#saveRDS(metadata, paste0(work_dir, "metadata_concat_df.Rds"))


#merged_concat_df = readRDS(paste0(file_dir, "merged_concat_df.Rds"))
#metadata = readRDS(paste0(work_dir, "metadata_concat_df.Rds"))

#write.table(merged_concat_df, paste0(work_dir, "merged_count_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
merged_concat_df = read.delim(paste0(file_dir, "merged_count_matrix.txt"))


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
#write.table(Cpm.count, paste0(work_dir, "merged_cpm_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
logcpm_corrected <- removeBatchEffect(Cpm.count, batch=sample_group_info)
#write.table(logcpm_corrected, paste0(work_dir, "merged_cpm_matrix_corrected.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)

pca_result <- prcomp(t(Cpm.count), scale. = FALSE)


pca_scores <- as.data.frame(pca_result$x)

library(ggsci)
p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() + 
  labs(title = "PCA of GSE111016, GSE167186",
       x = "Principal Component 1",
       y = "Principal Component 2")
p

#ggsave(paste0(work_dir, "raw_CPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)


pca_result <- prcomp(t(logcpm_corrected), scale. = FALSE)
pca_scores <- as.data.frame(pca_result$x)

p = ggplot(pca_scores, aes(PC1, PC2, color = as.character(metadata[2,]))) +
  geom_point() +
  theme_bw() + scale_color_nejm() +
  labs(title = "PCA of GSE111016, GSE167186",
       x = "Principal Component 1",
       y = "Principal Component 2")

p

#ggsave(paste0(work_dir, "raw_correctedCPM_pre_batch_correction.png"), width = 8, height = 7, plot = p)
library(GSVA)
msigdb_data <- msigdbr(species = "Homo sapiens")#, category = "H") # "H"는 hallmark gene sets
pathways <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)

expr_matrix <- as.matrix(logcpm_corrected)

ssgsea_results <- gsva(expr_matrix, pathways, method = "ssgsea")

# ssGSEA 결과 확인
head(ssgsea_results)

ssgsea_df = data.frame(t(ssgsea_results), "Annotation"=as.character(metadata[2,]))
table(as.character(metadata[2,]))
#group: Old Healthy      group: Sarcopenia    group: UNCLASSIFIED   group: Young Healthy  sarcopenia status: no sarcopenia status: yes 
#29                     24                      4                     19                     20                     20 
head(ssgsea_df[seq(1,10), seq(1,10)])
#write.xlsx(ssgsea_df, paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "ssGSEA", append=FALSE)


pvalue = apply(ssgsea_df[, seq(1, ncol(ssgsea_df)-1)], 2,
      function(x) {wilcox.test(x[ssgsea_df$Annotation %in% c("group: Sarcopenia", "sarcopenia status: yes")],
                               x[ssgsea_df$Annotation %in% c("group: Old Healthy")])$p.value})
pvalue = pvalue[order(pvalue)]

p = ggplot(ssgsea_df, aes(x=Annotation, y= HALLMARK_PROTEIN_SECRETION, fill=Annotation)) + geom_boxplot()
p = p + theme_bw()
p

write.xlsx(ssgsea_df, paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "ssGSEA")



batch <- factor(sample_group_info) 
Sample_type = as.character(metadata[2,])
Sample_type[Sample_type %in% c("group: Sarcopenia", "sarcopenia status: yes")] = "Sarcopenia"
Sample_type[Sample_type %in% c("group: Old Healthy", "sarcopenia status: no")] = "Healthy"
Sample_type[Sample_type %in% c("group: Young Healthy")] = "Young_Healthy"

# DGEList 

dge <- DGEList(counts=merged_concat_df,
             gene=rownames(merged_concat_df))
dge$samples$group <- factor(Sample_type)
dge$samples$batch <- batch

dge <- calcNormFactors(dge)

design <- model.matrix(~ batch + group, data=dge$samples)
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)

#  groupA - groupB
contrast <- makeContrasts(groupSarcopenia-groupHealthy, levels=design)

lrt <- glmLRT(fit, contrast=contrast)

topTags(lrt)
result_table <- topTags(lrt, n=Inf)  # 'n=Inf' 를 해서 전체 목록 도출

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
p
#ggsave(paste0(work_dir, "Enhanced_volanoPlot.png"), plot = p, width = 10, height = 9)

library(msigdbr)

deg_results <- topTags(lrt, n = Inf)$table

# 3. Prepare Gene List for GSEA
gene_list <- deg_results$logFC
names(gene_list) <- rownames(deg_results)
gene_list <- sort(gene_list, decreasing = TRUE)

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



library(ComplexHeatmap)
z_score = apply(logcpm_corrected, 1, function(x) {(x-mean(x))/sd(x)}) %>% data.frame() %>% t() %>% data.frame()

Heatmap(z_score[deg_results$genes[seq(1,50)],], name = "Corrected", 
        column_split = as.character(metadata[2,]))

write.xlsx(data.frame(gsea_results[seq(1,100),]), paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "GSEA", append = TRUE)

#openxlsx::write.xlsx(gsea_results, file = paste0(work_dir, "GSEA, ssGSEA results2.xlsx"))

