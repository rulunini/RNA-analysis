R.version.string
.libPaths()
library(argparse)

# 인수 파서 설정
parser <- ArgumentParser()
parser$add_argument("-i", "--input_directory", help="The results directory", required=TRUE)
parser$add_argument("-c", "--input_counts_cpm", help="The path of counts.csv", required=TRUE)
parser$add_argument("-m", "--input_matrix", help="The path of series matrix.txt", required=TRUE)
parser$add_argument("-o", "--output_directory", help="The results directory", required=TRUE)

args <- parser$parse_args()

# 인수를 사용하여 경로 설정
file_dir <- args$input_directory
counts_file <- args$input_counts_cpm
matrix_file <- args$input_matrix
work_dir <- args$output_directory

# 라이브러리 로드 =========
library(ggplot2)
library('msigdbr')
library(xlsx)
library(GSVA)
set.seed(100)

# 데이터 불러오기 =========
logcpm_corrected = read.delim(paste0(file_dir, "results_cpm_matrix_corrected.txt"))
metadata = readRDS(paste0(work_dir, "results_metadata_concat_df.Rds"))

# GSVA분석 ===========
msigdb_data <- msigdbr(species = "Homo sapiens")#, category = "H") # "H"는 hallmark gene sets
pathways <- split(msigdb_data$gene_symbol, msigdb_data$gs_name)

expr_matrix <- as.matrix(logcpm_corrected)

# 결과 탐색 및 요약
ssgsea_results <- gsva(expr_matrix, pathways, method = "ssgsea")
head(ssgsea_results)

ssgsea_df = data.frame(t(ssgsea_results), "Annotation"=as.character(metadata[1,]))
table(as.character(metadata[1,]))
#group: Old Healthy      group: Sarcopenia    group: UNCLASSIFIED   group: Young Healthy  sarcopenia status: no sarcopenia status: yes 
#29                     24                      4                     19                     20                     20 
head(ssgsea_df[seq(1,10), seq(1,10)])
#write.xlsx(ssgsea_df, paste0(work_dir, "GSEA, ssGSEA results.xlsx"), sheetName = "ssGSEA", append=FALSE)
write.csv(ssgsea_df, paste0(work_dir, "results_ssGSEA.csv"), row.names = FALSE)

# 시각화 및 통계 분석
#pvalue = apply(ssgsea_df[, seq(1, ncol(ssgsea_df)-1)], 2,
#      function(x) {wilcox.test(x[ssgsea_df$Annotation %in% c("group: Sarcopenia", "sarcopenia status: yes")],
#                               x[ssgsea_df$Annotation %in% c("group: Old Healthy")])$p.value})
#pvalue = pvalue[order(pvalue)]

#p = ggplot(ssgsea_df, aes(x=Annotation, y= HALLMARK_PROTEIN_SECRETION, fill=Annotation)) + geom_boxplot()
#p = p + theme_bw()
