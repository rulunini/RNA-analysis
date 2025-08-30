R.version.string
.libPaths()
start_time <- Sys.time()

library(argparse)

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

# 라이브러리 로드 =========

library(dplyr)
library(biomaRt)
library(tibble)
set.seed(100)

# (0) 데이터 로드 ===========

counts_df = read.delim(paste0(file_dir, counts_file), sep = ',')
counts_sheet = read.delim(paste0(file_dir, matrix_file), sep = ',')


# (1) 데이터 전처리  ===========
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

# 유전자 매핑 완료
merged_concat_df <- counts_df
metadata <- counts_sheet
print("Mapped Ensembl to Genesmbol")

# 데이터 저장
saveRDS(merged_concat_df, paste0(work_dir, "results_concat_df.Rds"))
saveRDS(metadata, paste0(work_dir, "results_metadata_concat_df.Rds"))
write.table(merged_concat_df, paste0(work_dir, "results_count_matrix.txt"), sep = '\t', col.names = TRUE, row.names=TRUE, quote = FALSE)
print('Saved Rds and txt')

end_time <- Sys.time()
print(end_time - start_time)
