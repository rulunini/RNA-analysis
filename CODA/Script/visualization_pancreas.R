library(ComplexHeatmap)
library(EnhancedVolcano)
library(dplyr)
library(VennDiagram)
library(grid) 

path <- "/Users/sohee/Desktop/02-GEO/Results/CODA/Results_Preprocessed/"
DEG <- read.csv(paste0(path, "pancreas_results_DEG_analysis.csv"))
output <- "/Users/sohee/Desktop/02-GEO/Results/CODA/Results_png/"

# Z-score 변환
logcpm_corrected <- read.csv(paste0(path, "pancreas_results_cpm.csv"), 
                             row.names = 1)
z_score <- apply(logcpm_corrected, 1, function(x) {
  (x - mean(x)) / sd(x)
}) %>% t() %>% as.data.frame()

# Metadata 생성
samples <- c("normal1", "normal2", "normal3", "normal4", 
             "patient1", "patient2", "patient3", "patient4", 
             "patient5", "patient6", "patient7", "patient8")
conditions <- ifelse(grepl("^normal", samples), "normal", "patient")
metadata <- data.frame(Sample = samples, 
                       condition = conditions, 
                       row.names = samples)
print(metadata)

# Volcano 생성 ------------------
p = EnhancedVolcano(DEG,
                    lab = DEG$Symbol,
                    x = 'logFC', 
                    y = 'PValue', 
                    FCcutoff = 1, 
                    pCutoff = 0.01,  
                    xlim = c(-10, 10),
                    # ylim = c(0, 10),
                    xlab = bquote(log[2] ~ "FC"),
                    ylab = bquote(-log[10] ~ "P-Value"),
                    pointSize = 2.0,
                    labSize = 3.0, 
                    legendPosition = 'bottom',
                    subtitle = NULL,
                    title = NULL)
p 

p2 = p + 
  annotate(
    "text", x = -7, y = 11, label = "n=2", size = 4) +
  annotate(
    "text", x = 7, y = 11, label = "n=4", size = 4) +
  annotate(
    "segment", x = -12, xend = -2, y = 11-0.5, yend = 11-0.5, 
    arrow = arrow(length = unit(0.3, "cm"), type = "closed", end = "both")) + 
  annotate(
    "segment", x = 12, xend = 2, y = 11-0.5, yend = 11-0.5, 
    arrow = arrow(length = unit(0.3, "cm"), type = "closed", end = "both")) 
p2 

ggsave(paste0(output, "results_volcanoplot.png"),
       plot = p2, dpi = 1200, width = 6.2, height = 6.4, units = 'in')


# VennDiagram 생성 ------------------
down_count <- nrow(DEG[DEG$logFC < -1 & DEG$PValue < 0.01, ])
up_count <- nrow(DEG[DEG$logFC > 1 & DEG$PValue < 0.01, ])
common_count <- nrow(DEG[DEG$PValue < 0.01, ]) - down_count - up_count  
grid.newpage()

venn.plot <- draw.pairwise.venn(
  area1 = down_count + common_count,  # Down-regulated 총 개수
  area2 = up_count + common_count,  # Up-regulated 총 개수
  cross.area = common_count,  # 겹치는 개수
  category = c("Down-regulated", "Up-regulated"),
  fill = c("lightblue", "lightpink"),
  alpha = 0.5,
  cat.pos = c(180, 0),
  cat.dist = 0.05,
  lwd = 2,
  cex = 1.5,  # 숫자 크기 조정
  cat.cex = 1,  # 라벨 크기 조정
  fontfamily = "Arial",  # 전체 글꼴 설정
  cat.fontfamily = "Arial"  # 카테고리 라벨 글꼴 설정
)
grid.draw(venn.plot)

# Heatmap 생성 ------------------
filtered_deg <- DEG[abs(DEG$logFC) > 1 & DEG$PValue < 0.01, ]
heatmap_df <- na.omit(z_score[filtered_deg$Symbol, ])
column_split <- as.character(metadata$condition)  
h = Heatmap(heatmap_df, 
            name = "Expression (Z-score)", 
            column_split = column_split,
            show_row_names = FALSE, 
            show_column_names = TRUE, 
            cluster_rows = TRUE,
            cluster_columns = FALSE) 
h

png(paste0(output, "results_heatmap.png"), 
    width = 6 * 1200, height = 3 * 1200, res = 1200)
draw(h)
dev.off()
