GSVA ÇØº¸±â
library(CellChat)
library(Seurat)

load('/home/sohee/analysis/data/HNSCC/mydata/seurat/v3/HNSCC_subcelltype.RData')

# set database #
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 


# (1) Create a CellChat object ----------


my.HPV <- c('HPV-','HPV+')
Idents(HNSCC) <- 'HPV'
subset <- subset(HNSCC, subset = HPV == my.HPV[2])
my.ch <- c('pos.pt','pos.lnmt')
for (i in 1:2){
  my.site <- c('PT','LNMT')
  
  Idents(subset) <- 'site'
  subset <- subset(subset, subset = site == my.site[1])
  assay <- GetAssayData(subset, assay = 'RNA', slot = 'data')
  
  Idents(subset) <- 'sub.celltype'
  labels <- Idents(subset)
  meta <- data.frame(labels, row.names = names(labels))
  cellchat <- createCellChat(object = assay, meta = meta, group.by = 'labels')
  cellchat <- setIdent(cellchat, ident.use = "labels") 
  cellchat@idents
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # (2) Cell to cell interaction ----------
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  # cellchat <- netAnalysis_computeCentrality(cellchat)

    assign(paste0(my.ch)[i], cellchat)
}


# (3) Merging cellchat object-------------------------------------------------- 
neg.pt, neg.lnmt
object.list <- list(PT = neg.pt, LNTM = neg.lnmt)
ch.neg.site <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
ch.neg.site <- setIdent(ch.neg.site, ident.use = 'labels')


pos.pt, pos.lnmt
object.list <- list(PT = pos.pt, LNTM = pos.lnmt)
ch.pos.site <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
ch.pos.site <- setIdent(ch.pos.site, ident.use = 'labels')



#  (4) rankNet ---------------------------------------------------------------------
cellchat.site <- ch.pos.site

gg1 <- rankNet(cellchat.site, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.site, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

gg3 <- rankNet(cellchat.site, mode = "comparison", 
               stacked = F, do.stat = TRUE,  cutoff.pvalue = 0.05, tol = 0.5,
               return.data = T, show.raw = T,
               title = 'CellChat signaling pathway (ligand-receptor pair) from primary and lymphnode'
               # color.use = c()
)
rankNET_pos_site_result <- gg3$signaling.contribution