library(Seurat)
reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Smooth_muscle_cell", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]

colSums(reference@assays$RNA@counts[grep('MT-',rownames(reference@assays$RNA@counts)),])/reference

PercentageFeatureSet(reference,pattern = 'MT-')->reference$percent.mt
reference<-reference[,which(reference$percent.mt<20)]
summary(reference$nCount_RNA)
reference<-reference[,which(reference$nCount_RNA>50)]
endo_sc<-reference[,which(reference$cell_type_high=='Vascular_endothelium')]
sc_pro(endo_sc,FALSE)->endo_sc
DimPlot(endo_sc,group.by = 'orig.ident')
library(harmony)
sc_pro(endo_sc,T)->endo_sc
DimPlot(endo_sc,group.by = 'orig.ident')
DimPlot(endo_sc,group.by = 'seurat_clusters')

endo_sc[,which(!endo_sc$seurat_clusters%in%c(2,9,11,12))]->object
sc_pro(object,T)->object
set.seed(100)
FindClusters(object,resolution = 0.8)->object
as.data.frame(prop.table(table(object$seurat_clusters,object$orig.ident),margin = 2))->a
a$condition<-reference$condition[match(a$Var2,reference$orig.ident)]
a_base<-a[grep('1',a$Var2),]
pvals<-data.frame()
for (cluster in 0:12) {
  a_base[which(a_base$Var1==cluster),]->a_test
  wilcox.test(a_test$Freq[which(a_test$condition=='SSC')],a_test$Freq[which(a_test$condition=='LS')])->test
  test$p.value->pvals[as.character(cluster),'p']
  mean(a_test$Freq[which(a_test$condition=='SSC')])-mean(a_test$Freq[which(a_test$condition=='LS')])->pvals[as.character(cluster),'t']
}

genes <- c(
  "ACKR1", "CLU", "PECAM1", "G0S2", "KLF4", "ANXA1", "CD36", "CA4", 
  "FABP4", "FOS", "DUSP1", "SELE", "SEMA3G", "ICAM2", "JAG1", "TIMP1", 
  "C4orf48", "CCL21", "LGALS5", "RGS5", "APOE", "KRT1", "KRT14", "SFN"
)
genes <- c("PCLAF", "PGF", "PROX1", "SEMA3G", "CD36", "ACKR1")
DotPlot(object,group.by = 'seurat_clusters',features = genes)

object[,grep('1',object$orig.ident)]->object_base
endo_marker_sc<-FindMarkers(object_base,ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',features = SSC_endo_marker,logfc.threshold = 0,min.pct = 0,test.use = 'MAST',latent.vars = 'orig.ident')
endo_marker_sc<-endo_marker_sc[which(endo_marker_sc$p_val_adj<0.05),]
table(sign(endo_marker_sc$avg_log2FC))

object@meta.data[,grep('SSC_marker',coobject_base_markerobject@meta.data[,grep('SSC_marker',colnames(object@meta.data))]<-NULL
AddModuleScore(object,features = list(rownames(DE_sc_endo_inter_SSC)),name = 'SSC_marker_inter')->object
VlnPlot(object,features = 'SSC_marker_inter1',group.by = 'seurat_clusters')

object_base<-object[,grep('1',object$orig.ident)]
as.data.frame(prop.table(table(object_base$seurat_clusters,object_base$condition),margin = 2))->a_condition
cbind(a_condition[1:13,],a_condition[14:26,])->combine
FindAllMarkers(object_base,only.pos = T)->object_base_marker
object_base_marker[which(object_base_marker$cluster==2),]->cluster2_marker
View(cluster2_marker[rownames(DE_sc_endo_inter_SSC),])
object_base_marker[which(object_base_marker$cluster==3),]->cluster3_marker
View(cluster3_marker[rownames(DE_sc_endo_inter_SSC),])
object_base_marker[which(object_base_marker$cluster==7),]->cluster7_marker
View(cluster7_marker[rownames(DE_sc_endo_inter_SSC),])

object_filter<-object[,which(!object$seurat_clusters%in%c(2,10,11,12))]
sc_pro<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda = 4,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}
sc_pro(object_filter,TRUE)->object_filter
object_filter<-AddModuleScore(object_filter,features = list(rownames(DE_sc_endo_inter_SSC)),name = 'SSC_marker_inter')
FeaturePlot(object_filter,'SSC_marker_inter1')
DimPlot(object_filter,label = T)
DimPlot(object_filter,label = T,group.by = 'condition')
DimPlot(object_filter,label = T,group.by = 'orig.ident')

as.data.frame(prop.table(table(object_filter$seurat_clusters[which(object_filter$orig.ident!='G2')],object_filter$orig.ident[which(object_filter$orig.ident!='G2')]),margin = 2))->a
a$condition<-reference$condition[match(a$Var2,reference$orig.ident)]
a_base<-a[grep('1',a$Var2),]
pvals<-data.frame()
for (cluster in 0:11) {
  a_base[which(a_base$Var1==cluster),]->a_test
  wilcox.test(a_test$Freq[which(a_test$condition=='SSC')],a_test$Freq[which(a_test$condition=='LS')])->test
  test$p.value->pvals[as.character(cluster),'p']
  mean(a_test$Freq[which(a_test$condition=='SSC')])-mean(a_test$Freq[which(a_test$condition=='LS')])->pvals[as.character(cluster),'t']
}
object_filter_base<-object_filter[,grep('1',object_filter$orig.ident)]
FindAllMarkers(object_filter_base,only.pos = T)->marker_object_filter_base
FindAllMarkers(object_filter_base,only.pos = T,features = rownames(DE_sc_endo_inter_SSC))->marker_object_filter_base_select
FindMarkers(object_filter_base,ident.1 = '2',group.by = 'seurat_clusters',features = rownames(DE_sc_endo_inter_SSC))
aggregate(marker_object_filter_base_select$avg_log2FC,by=list(marker_object_filter_base_select$cluster),sum)
aggregate(object_filter_base$SSC_marker_inter1,by=list(object_filter_base$seurat_clusters),mean)
VlnPlot(object_filter_base,group.by = 'seurat_clusters',features = 'SSC_marker_inter1')

AddModuleScore(object_filter_base,features = list(SSC_endo_marker),name = 'SSC_marker')->object_filter_base
aggregate(object_filter_base$SSC_marker1,by=list(object_filter_base$seurat_clusters),mean)
FindMarkers(object_filter_base,group.by = 'condition',ident.1 = 'SSC',ident.2 = 'LS',latent.vars = 'orig.ident',test.use = 'MAST',features = SSC_endo_marker,min.pct = 0,logfc.threshold = 0)->endo_marker

#endo fib co-cluster
reference[,grep('1',reference$orig.ident)]->reference_base
sc_pro(reference_base[,which(reference_base$cell_type_high%in%c('Vascular_endothelium','Fibroblasts'))],T)->endo_fib_base
AddModuleScore(endo_fib_base,features = list(rownames(DE_sc_endo_select_SSC)),name = 'SSC_marker_select')->endo_fib_base
FeaturePlot(endo_fib_base,'SSC_marker_select1')

merged_filter_exclude_epi_base<-merged_filter_exclude_epi[,grep('1',merged_filter_exclude_epi$orig.ident)]
sc_pro(merged_filter_exclude_epi_base[,which(merged_filter_exclude_epi_base$cell_type_l1_all%in%c('Endothelium','Fibroblasts'))],F)->endo_fib_st_base
DimPlot(endo_fib_st_base,group.by='condition')
DimPlot(endo_fib_st_base,group.by='cell_type_l1_all')
AddModuleScore(endo_fib_st_base,features = list(rownames(DE_sc_endo_select_SSC)),name = 'SSC_marker_select')->endo_fib_st_base
FeaturePlot(endo_fib_st_base,'SSC_marker_select1')
