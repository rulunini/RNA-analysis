library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(openxlsx)
library(ReactomePA)
library(msigdbr)
library(scales)  

path <- "/Users/sohee/Desktop/02-GEO/Results/CODA/Results_Preprocessed/"
output <- "/Users/sohee/Desktop/02-GEO/Results/CODA/Results_png/"

# DEG 리스트 준비
DEG <- read.csv(paste0(path, "pancreas_results_DEG_analysis.csv"))

# LogFC와 FDR을 변수로 받기
fdr_threshold <- 0.01 
logFC_threshold_up <- 2 
logFC_threshold_down <- -2 

# LogFC ≥ 2 (Up) & FDR < 0.05 유전자 추출
gene_symbols_up <- DEG$Symbol[DEG$logFC > logFC_threshold_up & DEG$PValue < fdr_threshold]
gene_entrez_up <- bitr(gene_symbols_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list_up <- gene_entrez_up$ENTREZID

# LogFC ≤ -2 (Down) & FDR < 0.05 유전자 추출
gene_symbols_down <- DEG$Symbol[DEG$logFC < logFC_threshold_down & DEG$PValue < fdr_threshold]
gene_entrez_down <- bitr(gene_symbols_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list_down <- gene_entrez_down$ENTREZID

# GO & KEGG 분석 함수 정의
perform_enrichment <- function(gene_list, ont) {
  if (length(gene_list) == 0) return(NULL)  
  enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH")
}

perform_kegg <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  enrichKEGG(gene = gene_list, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
}

perform_reactome <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  enrichPathway(gene = gene_list, organism = "human", pvalueCutoff = 0.05)
}

perform_msigdb <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  msigdb <- msigdbr(species = "Homo sapiens", category = "C2")  # MSigDB C2 (Canonical pathways) 사용
  enricher(gene_list, TERM2GENE = msigdb[, c("gs_name", "entrez_gene")], pvalueCutoff = 0.05)
}

perform_wikipathways <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  enrichWP(gene = gene_list, organism = "Homo sapiens", pvalueCutoff = 0.05)
}

# Up-regulated 분석 (먼저 실행!)
go_bp_up <- perform_enrichment(gene_list_up, "BP")
go_mf_up <- perform_enrichment(gene_list_up, "MF")
go_cc_up <- perform_enrichment(gene_list_up, "CC")
kegg_up <- perform_kegg(gene_list_up)
reactome_up <- perform_reactome(gene_list_up)
msigdb_up <- perform_msigdb(gene_list_up)
wikipathways_up <- perform_wikipathways(gene_list_up)

# Down-regulated 분석
go_bp_down <- perform_enrichment(gene_list_down, "BP")
go_mf_down <- perform_enrichment(gene_list_down, "MF")
go_cc_down <- perform_enrichment(gene_list_down, "CC")
kegg_down <- perform_kegg(gene_list_down)
reactome_down <- perform_reactome(gene_list_down)
msigdb_down <- perform_msigdb(gene_list_down)
wikipathways_down <- perform_wikipathways(gene_list_down)

# Excel 파일 생성 및 저장 ----------------
wb <- createWorkbook()

# 데이터를 엑셀에 추가하는 함수
add_to_excel <- function(wb, sheet_name, data) {
  if (!is.null(data)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(data))
  }
}

# Up-regulated 결과 저장 (먼저 저장!)
add_to_excel(wb, paste0("up_GO_BP"), go_bp_up)
add_to_excel(wb, paste0("up_GO_MF"), go_mf_up)
add_to_excel(wb, paste0("up_GO_CC"), go_cc_up)
add_to_excel(wb, paste0("up_KEGG"), kegg_up)
add_to_excel(wb, paste0("up_reactome"), reactome_up)
add_to_excel(wb, paste0("up_msigdb"), msigdb_up)
add_to_excel(wb, paste0("up_wikipathways"), wikipathways_up)

# Down-regulated 결과 저장 (나중에 저장!)
add_to_excel(wb, paste0("down_GO_BP"), go_bp_down)
add_to_excel(wb, paste0("down_GO_MF"), go_mf_down)
add_to_excel(wb, paste0("down_GO_CC"), go_cc_down)
add_to_excel(wb, paste0("down_KEGG"), kegg_down)
add_to_excel(wb, paste0("down_reactome"), reactome_down)
add_to_excel(wb, paste0("down_msigdb"), msigdb_down)
add_to_excel(wb, paste0("down_wikipathways"), wikipathways_down)

# Excel 저장
saveWorkbook(wb, paste0(path, "muscle_pathway_results_PValue_LogFC", logFC_threshold_up, ".xlsx"), overwrite = TRUE)

# Visualization ----------------
plot_go_enrichment <- function(go_data, title, showCategory = 10) {
  dotplot(go_data, x = "GeneRatio", showCategory = showCategory) + 
    scale_y_discrete(labels = wrap_format(100)) +  # 50자로 제한
    ggtitle(title)
}

p1 <- plot_go_enrichment(go_bp_up, "GO:BP of Up-Regulated Genes")
p2 <- plot_go_enrichment(go_mf_up, "GO:MF of Up-Regulated Genes")
p3 <- plot_go_enrichment(go_cc_up, "GO:CC of Up-Regulated Genes")
p4 <- plot_go_enrichment(kegg_up, "KEGG Pathway of Up-Regulated Genes")

p5 <- plot_go_enrichment(go_bp_down, "GO:BP of Down-Regulated Genes")
p6 <- plot_go_enrichment(go_mf_down, "GO:MF of Down-Regulated Genes")
p7 <- plot_go_enrichment(go_cc_down, "GO:CC of Down-Regulated Genes")
p8 <- plot_go_enrichment(kegg_down, "KEGG Pathway of Down-Regulated Genes")

# dotplot(kegg, showCategory=10) + ggtitle("KEGG Pathway Enrichment")

# Up-regulated 시각화 저장
ggsave(paste0(output, "up_go_bp.png"),
       plot = p1, dpi = 300, width = 8, height = 4)
ggsave(paste0(output, "up_go_mf.png"),
       plot = p2, dpi = 300, width = 12, height = 4)
ggsave(paste0(output, "up_go_cc.png"),
       plot = p3, dpi = 300, width = 7, height = 4)
ggsave(paste0(output, "up_kegg.png"),
       plot = p4, dpi = 300, width = 9, height = 4)

# down-regulated 시각화 저장
ggsave(paste0(output, "down_go_bp.png"),
       plot = p5, dpi = 300, width = 9, height = 4)
ggsave(paste0(output, "down_go_mf.png"),
       plot = p6, dpi = 300, width = 6, height = 4)
ggsave(paste0(output, "down_go_cc.png"),
       plot = p7, dpi = 300, width = 7, height = 4)
ggsave(paste0(output, "down_kegg.png"),
       plot = p8, dpi = 300, width = 6, height = 4)

