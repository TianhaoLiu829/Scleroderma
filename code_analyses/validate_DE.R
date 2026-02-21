merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
object_all<-merged_filter_exclude_epi
object_all<-object_all[,grep('1',object_all$orig.ident)]

#fibroblasts
object_fib<-object_all[,which(object_all$cell_type_l1_all=='Fibroblasts')]
mat_fib<-as.matrix(object_fib@assays$RNA@data)

DE_sc_fib_inter_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_LS.rds")
avg1 <- rowMeans(mat_fib[rownames(DE_sc_fib_inter_LS), which(object_fib$condition=='LS'), drop = FALSE])
avg2 <- rowMeans(mat_fib[rownames(DE_sc_fib_inter_LS), which(object_fib$condition=='SSC'), drop = FALSE])
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
all(log2fc>0)

DE_sc_fib_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_SSC.rds")
avg1 <- rowMeans(mat_fib[rownames(DE_sc_fib_inter_SSC), which(object_fib$condition=='LS'), drop = FALSE])
avg2 <- rowMeans(mat_fib[rownames(DE_sc_fib_inter_SSC), which(object_fib$condition=='SSC'), drop = FALSE])
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
all(log2fc<0)


#endothelial cells
object_endo<-object_all[,which(object_all$cell_type_l1_all=='Endothelium')]
mat_endo<-as.matrix(object_endo@assays$RNA@data)

DE_sc_endo_inter_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_endo_inter_LS.rds")
avg1 <- rowMeans(mat_endo[rownames(DE_sc_endo_inter_LS), which(object_endo$condition=='LS'), drop = FALSE])
avg2 <- rowMeans(mat_endo[rownames(DE_sc_endo_inter_LS), which(object_endo$condition=='SSC'), drop = FALSE])
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
all(log2fc>0)

DE_sc_endo_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_endo_inter_SSC.rds")
avg1 <- rowMeans(mat_endo[rownames(DE_sc_endo_inter_SSC), which(object_endo$condition=='LS'), drop = FALSE])
avg2 <- rowMeans(mat_endo[rownames(DE_sc_endo_inter_SSC), which(object_endo$condition=='SSC'), drop = FALSE])
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
all(log2fc<0)




#solve the problem of upregulated genes in LS in endo
#just manually calculate the fold change, the findmarker embedded method for calculation of P values are the same of manual in "wilcox". 
#But we want to use MAST then
FindMarkers(object_endo,group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'SSC',ident.2 = 'LS')->endo_marker
avg1 <- rowMeans(mat_endo[rownames(endo_marker), which(object_endo$condition=='SSC'), drop = FALSE])
avg2 <- rowMeans(mat_endo[rownames(endo_marker), which(object_endo$condition=='LS'), drop = FALSE])
endo_marker$avg_log2FC_manual<- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
endo_marker_LS<-endo_marker[which(endo_marker$avg_log2FC_manual<0&endo_marker$p_val_adj<0.05),]

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Smooth_muscle_cell", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
LS_endo_sc<-FindMarkers(reference[,which(reference$cell_type_l1=='Endothelium')],group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'SSC',ident.2 = 'LS',logfc.threshold = 0,min.pct = 0,features = rownames(endo_marker_LS))

cbind(endo_marker_LS,LS_endo_sc[rownames(endo_marker_LS),])->endo_marker_LS
colnames(endo_marker_LS)[1:5]<-paste0(colnames(endo_marker_LS)[1:5],'_st')
colnames(endo_marker_LS)[7:11]<-paste0(colnames(endo_marker_LS)[7:11],'_sc')
endo_marker_LS<-endo_marker_LS[which(sign(endo_marker_LS$avg_log2FC_st)==sign(endo_marker_LS$avg_log2FC_sc)),]
saveRDS(endo_marker_LS,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/endo_marker_LS.rds')
dot_plot(object_endo,rownames(endo_marker_LS)[order(0-endo_marker_LS$p_val_adj_st)],"condition")+theme(axis.text.y = element_text(size=12))->p
ggsave('p.png',p,width = 3.8,height = 5,dpi = 300)


#immune proximal and distal
merged_filter_exclude_epi_fib <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_fib.rds")
merged_filter_exclude_epi_endo <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_endo.rds")
immune_fib_marker<-FindMarkers(merged_filter_exclude_epi_fib,group.by = 'cell_type_l1_all',ident.1 = 'inf_fib',ident.2 = 'Fibroblasts',latent.vars = 'orig.ident',test.use = 'MAST')
immune_endo_marker<-FindMarkers(merged_filter_exclude_epi_endo,group.by = 'cell_type_l1_all',ident.1 = 'inf_endo',ident.2 = 'Endothelium',latent.vars = 'orig.ident',test.use = 'MAST')

object_endo_immune<-merged_filter_exclude_epi_endo[,which(merged_filter_exclude_epi_endo$cell_type_l1_all%in%c('Endothelium','inf_endo'))]
object_endo_immune<-object_endo_immune[,grep('1',object_endo_immune$orig.ident)]
mat_endo<-as.matrix(object_endo_immune@assays$RNA@data)
res_immune_endo_marker<-DE_wilcox(mat_endo,object_endo_immune$cell_type_l1_all,'inf_endo','Endothelium')
cbind(immune_endo_marker,res_immune_endo_marker[rownames(immune_endo_marker),])

DE_wilcox<-function(mat,group,ident1,ident2){
  g1 <- which(group == ident1)
  g2 <- which(group == ident2)
  
  # Preallocate results
  n_gene <- nrow(mat)
  pvals <- numeric(n_gene)
  avg1  <- numeric(n_gene)
  avg2  <- numeric(n_gene)
  
  # Compute group means (vectorized, fast)
  avg1 <- rowMeans(mat[, g1, drop = FALSE])
  avg2 <- rowMeans(mat[, g2, drop = FALSE])
  
  pct_1 <- rowMeans(mat[, g1, drop = FALSE]>0)
  pct_2 <- rowMeans(mat[, g2, drop = FALSE]>0)
  pct_all <- rowMeans(mat>0)
  
  # Compute Wilcoxon p-values (use apply over rows)
  pvals <- apply(mat, 1, function(x) {
    wilcox.test(x[g1], x[g2])$p.value
  })
  
  # Fold change
  log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
  
  # Combine results into data.frame
  res <- data.frame(
    gene   = rownames(mat),
    avg_group1 = avg1,
    avg_group2 = avg2,
    pct_1 = pct_1,
    pct_2 = pct_2,
    pct_all = pct_all,
    log2FC = log2fc,
    pvalue = pvals,
    padj   = p.adjust(pvals, method = "fdr")
  )
  return(res)
}
