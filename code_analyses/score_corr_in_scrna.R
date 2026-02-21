load("/ix1/wchen/liutianhao/work_rdata/work_10_05.RData")
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")

endo_marker_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/endo_marker_LS.rds")
DE_sc_endo_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_endo_inter_SSC.rds")
DE_sc_fib_inter_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_LS.rds")
DE_sc_fib_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_SSC.rds")

LS44 <- readRDS("/ix1/wchen/liutianhao/data/scleroderma_sc_bulk/single_cell_all/LS44.rds")



library(Seurat)
library(Matrix)
LS44$sample_celltype <- paste(LS44$library_id, LS44$celltype1, sep = "||")
seu<-LS44[,which(LS44$onset=='Peds')]
pb <- AggregateExpression(seu, group.by = "sample_celltype", slot = "counts", return.seurat = FALSE)$RNA
#filter the number of cells in each celltype-sample category (threshold can be adjusted)
pb<-pb[,names(table(seu$sample_celltype))[which(table(seu$sample_celltype)>20)]]
#meta data for each celltype-sample category
meta_pb <- data.frame(sample_celltype = colnames(pb))
tmp <- do.call(rbind, strsplit(meta_pb$sample_celltype, "\\|\\|"))
meta_pb$library_id <- tmp[,1]
meta_pb$celltype  <- tmp[,2]

scrna_score<-read.csv('/ix1/wchen/liutianhao/data/scleroderma_sc_bulk/clinical_data/scRNA-44samples-meta.csv')
scrna_score$Sample.ID[which(scrna_score$Sample.ID=='SC222AC')]<-'SC222'
clinical_info<-c('Sample.Type','Active..Inactive','Subtype',colnames(scrna_score)[20:24])
LS44@meta.data[,clinical_info]<-scrna_score[match(LS44$library_id,scrna_score$Sample.ID),clinical_info]
#borrow meta_pb and pb from pseudobulk analyses
meta_pb[,clinical_info]<-LS44@meta.data[match(meta_pb$sample_celltype,paste(LS44$library_id, seu$celltype1, sep = "||")),clinical_info]
meta_pb_case<-meta_pb[which(meta_pb$Subtype!='Healthy'),]
pb_case<-pb[,meta_pb_case$sample_celltype]
#may need to scale the gene expression raw counts

ct<-'Fibroblasts'
keep_cols <- meta_pb$celltype == ct
counts_ct <- pb[, keep_cols, drop = FALSE]
meta_ct   <- meta_pb[keep_cols, , drop = FALSE]
sample_meta <- unique(seu@meta.data[, c("library_id", "health")])
meta_ct <- merge(meta_ct, sample_meta, by = "library_id", all.x = TRUE)
rownames(meta_ct) <- meta_ct$sample_celltype

counts_ct <- counts_ct[, rownames(meta_ct)]
meta_ct$health<-factor(meta_ct$health,levels = c('Healthy','LS'))

dge <- DGEList(counts_ct)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE)
colMeans(logCPM[c('CCL18','CD74','B2M','CXCL9','C3'),])->meta_ct$score
log1p(counts_ct/colSums(counts_ct)*10000)->counts_ct
colMeans(counts_ct[rownames(marker_all)[which(marker_all$avg_log2FC>0)][1:20],])->meta_ct$score
meta_ct_case<-meta_ct[which(meta_ct$health=='LS'),]
summary(aov(score~Active..Inactive,meta_ct))->test

cor.test(meta_ct$score,as.numeric(meta_ct$PGA.A),method = 'spearman')

AddModuleScore()

log1p(counts_ct/colSums(counts_ct)*10000)-PGA.Alog1p(counts_ct/colSums(counts_ct)*10000)->counts_ct

cor(colMeans(counts_ct[c('CCL18','CD74','B2M','CXCL9','C3'),]),meta_ct$LoSAI..mLoSSI)

