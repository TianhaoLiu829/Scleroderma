library(harmony)
library(Seurat)

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_fibthelium" = "fibthelium","Lymphatic_fibthelium" = "fibthelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "fibthelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_fibthelium"="fibthelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]


PercentageFeatureSet(reference,pattern = 'MT-')->reference$percent.mt
reference[,which(reference$cell_type_high=='Fibroblasts'&reference$percent.mt<20&reference$nCount_RNA<20000&reference$nCount_RNA>200)]->fib_sc
I1<-Read10X_h5('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/all_data/flex/result/outs/per_sample_outs/NRCOS356/count/sample_filtered_feature_bc_matrix.h5')
fib_sc<-fib_sc[rownames(I1),]
sc_pro3(fib_sc,T)->fib_sc_harmony
DimPlot(fib_sc_harmony,label = T)
FindAllMarkers(fib_sc_harmony,only.pos = T)->marker_all
fib_sc_harmony<-fib_sc_harmony[,which(!fib_sc_harmony$seurat_clusters%in%c(8,15,16))]

AddModuleScore(fib_sc_harmony,features = list(c('CD74','C3','B2M','LYZ')),name = 'inf_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'inf_fib_marker1')
AddModuleScore(fib_sc_harmony,features = list(c("COL1A1","COL3A1","COL1A2","SPARC","CXCL14")),name = 'other_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'other_fib_marker1')

fib_sc_harmony<-fib_sc_harmony[,which(fib_sc_harmony@reductions$umap@cell.embeddings[,1]>0-10)]

for (resolution in seq(0.3,2,0.1)) {
  resolution<-1.7
  FindClusters(fib_sc_harmony,resolution = resolution,seed=100)->fib_sc_harmony
  DimPlot(fib_sc_harmony[,which(fib_sc_harmony@reductions$umap@cell.embeddings[,2]<7.5&fib_sc_harmony@reductions$umap@cell.embeddings[,1]<7)],label = T)->p
  aggregate(fib_sc_harmony@meta.data$inf_fib_marker1,by=list(fib_sc_harmony$seurat_clusters),mean)->inf_fib
  aggregate(fib_sc_harmony@meta.data$other_fib_marker1,by=list(fib_sc_harmony$seurat_clusters),mean)->other_fib
  cbind(inf_fib$x,other_fib$x)->a
  colnames(a)<-c('inf_fib','other_fib')
  rownames(a)<-inf_fib$Group.1
  pheatmap::pheatmap(a,show_rownames = T,scale = 'column',show_colnames = T,border_color=NA)->p
}
test_all<-as.data.frame(matrix(ncol = 2))
as.data.frame(prop.table(table(fib_sc_harmony$orig.ident,fib_sc_harmony$immune_proximal),margin = 1))->df
df$PGAA<-meta_ct_case$PGA.A[match(df,meta_ct_case$library_id)]
for (cluster in 0:23) {
  fib_sc_harmony$immune_proximal<-ifelse(fib_sc_harmony$seurat_clusters%in%c(0,3,8,15,14,13,16),'inf_fib','other_fib')
  fib_sc_harmony_base<-fib_sc_harmony[,grep('1',fib_sc_harmony$orig.ident)]
  prop.table(table(fib_sc_harmony_base$immune_proximal,fib_sc_harmony_base$condition),margin = 2)->df
  prop.table(table(fib_sc_harmony_base$immune_proximal,fib_sc_harmony_base$orig.ident),margin = 2)->df
  df<-as.data.frame(t(df))
  df$condition<-reference$condition[match(df$Var1,reference$orig.ident)]
  df<-df[which(df$Var2=='inf_fib'),]
  df<-df[which(!df$Var1%in%names(table(fib_sc_harmony$orig.ident))[which(table(fib_sc_harmony$orig.ident)<50)]),]
  t.test(Freq~condition,df)->test
  ggplot(df, aes(x = condition, y = Freq, fill = condition)) +
    geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
    geom_jitter(width = 0.1, size = 1, color = "black") + # raw points
    theme_classic() +
    labs(title = "Frequency of inflammed fibthelium")->p
  set.seed(150)
  ggsave('p1.png',p,width = 4,height = 4)
  c(test$p.value,test$estimate)->test_all[cluster,]
}
DimPlot(fib_sc_harmony[,which(fib_sc_harmony@reductions$umap@cell.embeddings[,1]<7&fib_sc_harmony@reductions$umap@cell.embeddings[,2]<7.5)],group.by = 'RNA_snn_res.1.7',label = T)->p
ggsave('p.png',p,width = 6,height = 5,dpi = 300)

FindMarkers(fib_sc_harmony,ident.1 = 'inf_fib',ident.2 = 'other_fib',only.pos = T,group.by = 'immune_proximal')->marker_proximal_fib
out<-DE_wilcox(fib_sc_harmony@assays$RNA@data,fib_sc_harmony$immune_proximal,ident1 = 'inf_fib',ident2 = 'other_fib')
marker_proximal_fib$avg_log2FC<-out$log2FC[match(rownames(marker_proximal_fib),rownames(out))]
res_inf_fib <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/res_inf_fib.rds")
gene_module_fib<-intersect(res_inf_fib$gene[which(res_inf_fib$log2FC>0-5&res_inf_fib$log2FC<0&res_inf_fib$padj<0.05)],rownames(marker_proximal_fib)[which(marker_proximal_fib$avg_log2FC>0&marker_proximal_fib$avg_log2FC<10&marker_proximal_fib$p_val_adj<0.05)])
gene_module_fib<-c(gene_module_fib,'CXCL9')

FindMarkers(endo_sc_harmony,ident.1 = 'inf_endo',ident.2 = 'other_endo',only.pos = T,group.by = 'immune_proximal')->marker_proximal_endo
out<-DE_wilcox(endo_sc_harmony@assays$RNA@data,endo_sc_harmony$immune_proximal,ident1 = 'inf_endo',ident2 = 'other_endo')
marker_proximal_endo$avg_log2FC<-out$log2FC[match(rownames(marker_proximal_endo),rownames(out))]
res_inf_endo <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/res_inf_endo.rds")
gene_module_endo<-intersect(res_inf_endo$gene[which(res_inf_endo$log2FC>0-5&res_inf_endo$log2FC<0&res_inf_endo$padj<0.05)],rownames(marker_proximal_endo)[which(marker_proximal_endo$avg_log2FC>0&marker_proximal_endo$avg_log2FC<10&marker_proximal_endo$p_val_adj<0.05)])
gene_module_endo[2]<-'CXCR4'

LS44 <- readRDS("/ix1/wchen/liutianhao/data/scleroderma_sc_bulk/single_cell_all/LS44.rds")
LS44$sample_celltype <- paste(LS44$library_id, LS44$celltype1, sep = "||")
seu<-LS44[,which(LS44$onset=='Peds')]
pb <- AggregateExpression(seu, group.by = "sample_celltype", slot = "data", return.seurat = FALSE)$RNA
#filter the number of cells in each celltype-sample category (threshold can be adjusted)
pb<-pb[,names(table(seu$sample_celltype))[which(table(seu$sample_celltype)>20)]]
#meta data for each celltype-sample category
meta_pb <- data.frame(sample_celltype = colnames(pb))
meta_pb[,c('library_id','celltype')] <- as.matrix(do.call(rbind, strsplit(meta_pb$sample_celltype, "\\|\\|")))

scrna_score<-read.csv('/ix1/wchen/liutianhao/data/scleroderma_sc_bulk/clinical_data/scRNA-44samples-meta-v2.csv')
scrna_score$Sample.ID[which(scrna_score$Sample.ID=='SC222AC')]<-'SC222'
clinical_info<-c('Sample.Type','Active..Inactive','Subtype','Onset',colnames(scrna_score)[16:20])
#borrow meta_pb and pb from pseudobulk analyses
meta_pb[,clinical_info]<-scrna_score[match(meta_pb$library_id,scrna_score$Sample.ID),clinical_info] 
meta_pb_case<-meta_pb[which(meta_pb$Subtype!='Healthy'),]
pb_case<-pb[,meta_pb_case$sample_celltype]


ct<-'Fibroblasts'
keep_cols <- meta_pb_case$celltype == ct
counts_ct_case <- pb_case[, keep_cols, drop = FALSE]
meta_ct_case   <- meta_pb_case[keep_cols, , drop = FALSE]
counts_ct_case<-log1p(counts_ct_case/colSums(counts_ct_case)*10000)
t(scale(t(counts_ct_case)))->counts_ct_case
colMeans(counts_ct_case[intersect(gene_module_fib,rownames(counts_ct_case)),])->meta_ct_case$score
#dge <- DGEList(counts = counts_ct_case)
#dge <- calcNormFactors(dge)
#cpm(dge, log = TRUE, prior.count = 1)->cpm
#colMeans(cpm[intersect(rownames(cpm),gene_module_fib),])->meta_ct_case$score

ct<-'Endothelial Cells'
keep_cols <- meta_pb_case$celltype == ct
counts_ct_case <- pb_case[, keep_cols, drop = FALSE]
meta_ct_case   <- meta_pb_case[keep_cols, , drop = FALSE]
counts_ct_case<-log1p(counts_ct_case/colSums(counts_ct_case)*10000)
t(scale(t(counts_ct_case)))->counts_ct_case
colMeans(counts_ct_case[intersect(gene_module_endo,rownames(counts_ct_case)),])->meta_ct_case$score



merged_filter_exclude_epi$sample_celltype <- paste(merged_filter_exclude_epi$orig.ident, merged_filter_exclude_epi$cell_type_l1_all, sep = "||")
seu<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$condition=='LS')]
pb <- AggregateExpression(seu, group.by = "sample_celltype", slot = "data", return.seurat = FALSE)$RNA
#filter the number of cells in each celltype-sample category (threshold can be adjusted)
#pb<-pb[,names(table(seu$sample_celltype))[which(table(seu$sample_celltype)>20)]]
#meta data for each celltype-sample category
meta_pb <- data.frame(sample_celltype = colnames(pb))
meta_pb[,c('library_id','celltype')] <- as.matrix(do.call(rbind, strsplit(meta_pb$sample_celltype, "\\|\\|")))

meta_pb[,c('PGAA','PGAD')]<-cbind(PGAA[meta_pb$library_id],PGAD[meta_pb$library_id])
meta_pb_case<-meta_pb
pb_case<-pb

ct<-'Fibroblasts'
keep_cols <- meta_pb_case$celltype == ct
counts_ct_case <- pb_case[, keep_cols, drop = FALSE]
meta_ct_case   <- meta_pb_case[keep_cols, , drop = FALSE]
counts_ct_case<-log1p(counts_ct_case/colSums(counts_ct_case)*10000)
t(scale(t(counts_ct_case)))->counts_ct_case
colMeans(counts_ct_case[intersect(gene_module_fib,rownames(counts_ct_case)),])->meta_ct_case$score




target<-'fib'
setwd('~/part_endo/')
library(ggpubr)
df_exp<-data.frame(score=meta_ct_case$score,PGAA=meta_ct_case$PGA.A)
cor.test(df_exp[which(df_exp$PGAA!=0),1],df_exp[which(df_exp$PGAA!=0),2],method = 'spearman')->test
ggplot(df_exp, aes(x = score, y = PGAA)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  stat_cor(method = "spearman",
           label.x = min(df_exp$score, na.rm = TRUE),
           label.y = max(df_exp$losdi, na.rm = TRUE),hjust = 0,vjust = -17.5,
           size = 5) +
  labs(x = "X value", y = "Y value")->p
ggsave(paste0('p_PGAA_',target,'.png'),p,width = 6,height = 5,dpi = 300)

df_exp<-data.frame(score=meta_ct_case$score,PGAD=meta_ct_case$PGA.D)
cor.test(df_exp[,1],df_exp[,2],method = 'spearman')->test
ggplot(df_exp[which(meta_ct_case$Active..Inactive=='A'),], aes(x = score, y = PGAD)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  stat_cor(method = "spearman",
           label.x = min(df_exp$score, na.rm = TRUE),
           label.y = max(df_exp$losdi, na.rm = TRUE),hjust = 0,vjust = -17,
           size = 5) +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave(paste0('p_PGAD_',target,'.png'),p,width = 6,height = 5,dpi = 300)

df<-data.frame(PGAA=c(test_PGAA_fib$estimate,test_PGAA_endo$estimate),PGAD=c(test_PGAD_fib$estimate,test_PGAD_endo$estimate))
pheatmap::pheatmap(df,cluster_rows = F,cluster_cols = F,breaks = seq(-abs(max(as.matrix(df))),abs(max(as.matrix(df))),length.out=101))->p
pheatmap::pheatmap(df,cluster_rows = F,cluster_cols = F)->p
ggsave('p.png',p,width = 4,height = 3,dpi = 300)


df_exp<-data.frame(score=meta_ct_case$score,lossi=meta_ct_case$LoSAI..mLoSSI)
ggplot(df_exp[which(meta_ct_case$Active..Inactive=='A'),], aes(x = score, y = lossi)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  stat_cor(method = "spearman",
           label.x = min(df_exp$score, na.rm = TRUE),
           label.y = max(df_exp$losdi, na.rm = TRUE),hjust = 0,vjust = -17,
           size = 5) +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave(paste0('p_mLoSSI_',target,'.png'),p,width = 6,height = 5,dpi = 300)

df_exp<-data.frame(score=meta_ct_case$score,losdi=meta_ct_case$LoSDI)
ggplot(df_exp[which(meta_ct_case$Active..Inactive=='A'),], aes(x = score, y = losdi)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  stat_cor(method = "spearman",
           label.x = min(df_exp$score, na.rm = TRUE),
           label.y = max(df_exp$losdi, na.rm = TRUE),hjust = 0,vjust = 2,
           size = 5) +
  labs(x = "X value", y = "Y value")->p
ggsave(paste0('p_losdi_',target,'.png'),p,width = 6,height = 5,dpi = 300)

df_exp<-data.frame(score=meta_ct_case$score,lossi=meta_ct_case$PGA.S)
ggplot(df_exp[which(meta_ct_case$Active..Inactive=='A'),], aes(x = score, y = lossi)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  stat_cor(method = "spearman",
           label.x = min(df_exp$score, na.rm = TRUE),
           label.y = max(df_exp$losdi, na.rm = TRUE),hjust = 0,vjust = -17,
           size = 5) +
  labs(x = "X value", y = "Y value")->p
ggsave(paste0('p_PGAS_',target,'.png'),p,width = 6,height = 5,dpi = 300)



design <- model.matrix(~health, data = meta_ct)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
res <- glmQLFTest(fit, coef = 'healthLS')
topTags(res)
de <- as.data.frame(res$table)
de$gene <- rownames(de)
de$FDR<-p.adjust(de$PValue,method='fdr')
de->res_list[[ct]]




Idents(fib_sc_harmony)<-fib_sc_harmony$RNA_snn_res.1.6
FindAllMarkers(fib_sc_harmony,only.pos = T)->marker_all_fib

#confirm the spatial genes to plot (only genes validated with scRNA)
Idents(fib_sc_harmony)<-fib_sc_harmony$RNA_snn_res.1.6
FindAllMarkers(fib_sc_harmony,only.pos = T)->marker_all_fib

res_fib_pos<-res_fib[which(res_fib$log2FC>0.4),]
marker_all_fib2<-marker_all_fib[which(marker_all_fib$avg_log2FC>0.4),]
marker_all_fib2$p_val_adj<-p.adjust(marker_all_fib2$p_val,method = 'bonferroni')
intersect(rownames(res_fib_pos[which(p.adjust(res_fib_pos$pvalue,method = 'bonferroni')<0.05),]),marker_all_fib2$gene[which(marker_all_fib2$p_val_adj<0.05)])->inf_features

res_fib_neg<-res_fib[which(res_fib$log2FC<0-0.4),]
intersect(rownames(res_fib_neg[which(p.adjust(res_fib_neg$pvalue,method = 'bonferroni')<0.05),]),marker_all_fib2$gene[which(marker_all_fib2$p_val_adj<0.05)])->other_features

features<-c('CD74','C3','CCL19','B2M','CXCL9',"COL1A1","COL3A1","COL1A2","SPARC","DSP")
dot_plot(merged_filter_exclude_epi_fib,features,'inf_fib')+coord_flip()+theme(axis.text.x = element_text(size=13))->p
ggsave('p_inf_fib.png',p,width = 6,height = 3,dpi = 300)

#trash code
AddModuleScore(fib_sc_harmony,features = list(intersect(rownames(DE_sc_fib_inter_LS),inf_features)[-3]),name = 'inf_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'inf_fib_marker1')
AddModuleScore(fib_sc_harmony,features = list(other_features[c(1,2,3,4,6)]),name = 'other_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'other_fib_marker1')

AddModuleScore(fib_sc_harmony,features = list(inf_features[c(1,2,4,5,7)]),name = 'inf_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'inf_fib_marker1')
AddModuleScore(fib_sc_harmony,features = list(other_features[c(1,2,3,4,6)]),name = 'other_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'other_fib_marker1')


AddModuleScore(fib_sc_harmony,features = list(c("CCL19","C3","CD74","B2M","CXCL9")),name = 'inf_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'inf_fib_marker1')
AddModuleScore(fib_sc_harmony,features = list(other_features[c(1,2,3,4,6)]),name = 'other_fib_marker')->fib_sc_harmony
FeaturePlot(fib_sc_harmony,'other_fib_marker1')




for (variable in vector) {
  fib_sc_harmony_base<-fib_sc_harmony[,grep('1',fib_sc_harmony$orig.ident)]
  fib_sc_harmony_base$immune_proximal[which(colnames(fib_sc_harmony_base)%in%colnames(fib_sc_harmony)[which(fib_sc_harmony$seurat_clusters%in%c(1,9,13,15))])]<-'inf_fib'
  fib_sc_harmony_base$immune_proximal[which(colnames(fib_sc_harmony_base)%in%colnames(fib_sc_harmony)[which(fib_sc_harmony$seurat_clusters%in%c(11,0,8,7,5,6))])]<-'other_fib'
  prop.table(table(fib_sc_harmony_base$immune_proximal,fib_sc_harmony_base$condition),margin = 2)->df
  prop.table(table(fib_sc_harmony_base$immune_proximal,fib_sc_harmony_base$orig.ident),margin = 2)->df
  df<-as.data.frame(t(df))
  df$condition<-reference$condition[match(df$Var1,reference$orig.ident)]
  df<-df[which(df$Var2=='inf_fib'),]
  t.test(Freq~condition,df)
  ggplot(df, aes(x = condition, y = Freq, fill = condition)) +
    geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
    geom_jitter(width = 0.1, size = 1.5, color = "black") + # raw points
    theme_classic() +
    labs(title = "Frequency of inflammed fibthelium")->p
  ggsave('p.png',p,width = 4,height = 4)
  pheatmap::pheatmap(df)
}
# FeaturePlot(fib_sc_harmony,'other_fib_marker1')
marker_inf_fib[which(marker_inf_fib$p_val_adj<0.05),]->marker_inf_fib_sig

VlnPlot(fib_sc_harmony,group.by = 'seurat_clusters',features = 'inf_fib_marker1',pt.size = 0)
fib_sc_harmony<-fib_sc_harmony[,which(!fib_sc_harmony$seurat_clusters%in%c(16,12,17))]
fib_sc_harmony<-fib_sc_harmony[,which(fib_sc_harmony@reductions$umap@cell.embeddings[,1]>0-8&fib_sc_harmony@reductions$umap@cell.embeddings[,2]>0-8)]
set.seed(100)
FindClusters(fib_sc_harmony,resolution = 0.65)->fib_sc_harmony
aggregate(fib_sc_harmony_base@meta.data$inf_fib_marker1,by=list(fib_sc_harmony$seurat_clusters),mean)->inf_fib
aggregate(fib_sc_harmony_base@meta.data$other_fib_marker1,by=list(fib_sc_harmony$seurat_clusters),mean)->other_fib

cbind(inf_fib$x,other_fib$x)->a
colnames(a)<-c('inf_fib','other_fib')
rownames(a)<-inf_fib$Group.1
pheatmap::pheatmap(a,scale = 'column',show_rownames = T,show_colnames = T,border_color=NA)->p
ggsave('p.png',p,width=3,height = 5,dpi = 300)
DimPlot(fib_sc_harmony,label = T)->p
ggsave('p1.png',p,width=5,height = 5,dpi = 300)
aggregate()



sc_pro2(fib_sc_harmony,F)->fib_sc_harmony2
AddModuleScore(fib_sc_harmony2,features = list(rownames(marker_inf_fib_sig)[which(marker_inf_fib_sig$avg_log2FC>0)]),name = 'inf_fib_marker')->fib_sc_harmony2
AddModuleScore(fib_sc_harmony2,features = list(rownames(marker_inf_fib_sig)[which(marker_inf_fib_sig$avg_log2FC<0)]),name = 'other_fib_marker')->fib_sc_harmony2
VlnPlot(fib_sc_harmony,group.by = 'seurat_clusters',features = 'inf_fib_marker1',pt.size = 0)
aggregate(fib_sc_harmony2@meta.data$inf_fib_marker1,by=list(fib_sc_harmony2$seurat_clusters),mean)->inf_fib
aggregate(fib_sc_harmony2@meta.data$other_fib_marker1,by=list(fib_sc_harmony2$seurat_clusters),mean)->other_fib


fib_sc_harmony@assays$RNA@data[rownames(marker_inf_fib_sig),which(!fib_sc_harmony$seurat_clusters%in%c(16,12,17))]->mat
fib_sc_harmony@assays$RNA@data[features[1:5],which(!fib_sc_harmony$seurat_clusters%in%c(16,12,17))]->mat
mat_scaled <- t(scale(mat)) # transpose so rows = samples

d <- dist(mat_scaled, method = "euclidean")
hc <- hclust(d, method = "ward.D2")# Plot dendrogram
clusters <- cutree(hc, k = 2)
pheatmap::pheatmap(mat_scaled[order(clusters),rownames(marker_inf_fib_sig)[order(marker_inf_fib_sig$avg_log2FC)]],cluster_row=F,cluster_col=F,scale='column')
apply(mat_scaled[,rownames(marker_inf_fib_sig)[which(marker_inf_fib_sig$avg_log2FC<0)]], 2, function(x) {
  aggregate(x,by=list(clusters),mean)
  t.test(x[which(clusters==1)],x[which(clusters==2)])
}) 


set.seed(123)
km <- kmeans(mat_scaled, centers = 2, nstart = 50)
clusters <- km$cluster
pheatmap::pheatmap(mat_scaled[order(clusters),rownames(marker_inf_fib_sig)[order(marker_inf_fib_sig$avg_log2FC)]],cluster_row=F,cluster_col=F,scale='column',annotation_row=data.frame(row.names=rownames(mat_scaled),cluster=clusters))
pheatmap::pheatmap(mat_scaled[order(clusters),features[1:5]],cluster_row=F,cluster_col=F,scale='column',annotation_row=data.frame(row.names=rownames(mat_scaled),cluster=clusters))
apply(mat_scaled[,features], 2, function(x) {
  aggregate(x,by=list(clusters),mean)
  t.test(x[which(clusters==1)],x[which(clusters==2)])
}) 


sc_pro<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object,features=rownames(object))
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda =1,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}

sc_pro<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object,vars.to.regress='nCount_RNA',features=rownames(object))
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda =1,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 15)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}


sc_pro2<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 1,lambda =1.5,sigma = 0.06)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 15)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}


sc_pro3<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda =1,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}

