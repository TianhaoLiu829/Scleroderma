library(harmony)
library(Seurat)

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_fibthelium" = "fibthelium","Lymphatic_fibthelium" = "fibthelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "fibthelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_fibthelium"="fibthelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
PercentageFeatureSet(reference,pattern = 'MT-')->reference$percent.mt

reference[,which(reference$cell_type_high=='Vascular_endothelium'&reference$percent.mt<20&reference$nCount_RNA<20000&reference$nCount_RNA>200)]->endo_sc
I1<-Read10X_h5('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/all_data/flex/result/outs/per_sample_outs/NRCOS356/count/sample_filtered_feature_bc_matrix.h5')
sc_pro3(endo_sc[rownames(I1),],T)->endo_sc_harmony
DimPlot(endo_sc_harmony,group.by = 'seurat_clusters',label = T)
DimPlot(endo_sc_harmony,group.by = 'orig.ident',label = T)
#FindMarkers(endo_sc_harmony,only.pos = T,logfc.threshold = 0.1,ident.1 = c(6,9))->marker_endo_69
#FindMarkers(endo_sc_harmony,only.pos = T,logfc.threshold = 0.1,ident.1 = c(11))->marker_endo_11

endo_sc_harmony<-endo_sc_harmony[,which(!endo_sc_harmony$seurat_clusters%in%c(6,9))]
sc_pro4(endo_sc_harmony,T)->endo_sc_harmony
DimPlot(endo_sc_harmony,group.by = 'seurat_clusters',label = T)
#DimPlot(endo_sc_harmony,group.by = 'orig.ident',label = T)
#DimPlot(endo_sc_harmony,group.by = 'condition',label = T)

#FindAllMarkers(endo_sc_harmony,only.pos = T,logfc.threshold = 0.1)->marker_all_endo

#res_endo_pos<-res_endo[which(res_endo$log2FC>0.4),]
#marker_all_endo2<-marker_all_endo[which(marker_all_endo$avg_log2FC>0.4),]
#marker_all_endo2$p_val_adj<-p.adjust(marker_all_endo2$p_val,method = 'bonferroni')
#intersect(rownames(res_endo_pos[which(res_endo_pos$pvalue<0.05),]),marker_all_endo2$gene[which(marker_all_endo2$p_val<0.05)])->inf_features

#res_endo_neg<-res_endo[which(res_endo$log2FC<0-0.4),]
#intersect(rownames(res_endo_neg[which(res_endo_neg$pvalue<0.05),]),marker_all_endo2$gene[which(marker_all_endo2$p_val<0.05)])->other_features



#features <- c("CD74","CCL19","LYZ","CXCR4","FOS","FABP4","CSRP1","CAVIN1","COL4A1","COL4A2")
#features <- c("CD74","CXCR4","SERPINE1","FOS","C7","FABP4","CD36","CAVIN1","COL4A1","COL4A2")
features <- c("CD74","CXCR4","SERPINE1","FOS","SELE","FABP4","CD36","CAVIN1","COL4A1","COL4A2")
#features <- c("CD74","CXCR4","SERPINE1","FOS","SELE","FABP4","CD36","CAVIN1","FABP5","ABLIM3")
dot_plot(merged_filter_exclude_epi_endo,features,'inf_endo')+coord_flip()+theme(axis.text.x = element_text(size=13))->p
ggsave('p_inf_endo.png',p,width = 6,height = 3,dpi = 300)


AddModuleScore(endo_sc_harmony,features = list(features[1:5]),name = 'inf_endo_marker')->endo_sc_harmony
FeaturePlot(endo_sc_harmony,'inf_endo_marker1')
AddModuleScore(endo_sc_harmony,features = list(features[6:10]),name = 'other_endo_marker')->endo_sc_harmony
FeaturePlot(endo_sc_harmony,'other_endo_marker1')

#AddModuleScore(endo_sc_harmony,features = list(inf_features),name = 'inf_endo_marker')->endo_sc_harmony
#FeaturePlot(endo_sc_harmony,'inf_endo_marker1')
#AddModuleScore(endo_sc_harmony,features = list(other_features),name = 'other_endo_marker')->endo_sc_harmony
#eaturePlot(endo_sc_harmony,'other_endo_marker1')


for (resolution in seq(0.3,2,0.1)) {
  resolution<-2
  FindClusters(endo_sc_harmony,resolution = resolution,seed=100)->endo_sc_harmony
  DimPlot(endo_sc_harmony,label = T)
  aggregate(endo_sc_harmony@meta.data$inf_endo_marker1,by=list(endo_sc_harmony$seurat_clusters),mean)->inf_endo
  aggregate(endo_sc_harmony@meta.data$other_endo_marker1,by=list(endo_sc_harmony$seurat_clusters),mean)->other_endo
  cbind(inf_endo$x,other_endo$x)->a
  colnames(a)<-c('inf_endo','other_endo')
  rownames(a)<-inf_endo$Group.1
  pheatmap::pheatmap(a,show_rownames = T,scale = 'column',show_colnames = T,border_color=NA)->p
  ggsave('p.png',p,width = 3.5,height = 6,dpi = 300)
  print(p)
}

endo_sc_harmony$immune_proximal<-ifelse(endo_sc_harmony$RNA_snn_res.2%in%c(0,2,3,10,20),'inf_endo','other_endo')

endo_sc_harmony_base<-endo_sc_harmony[,grep('1',endo_sc_harmony$orig.ident)]
prop.table(table(endo_sc_harmony_base$immune_proximal,endo_sc_harmony_base$condition),margin = 2)->df
prop.table(table(endo_sc_harmony_base$immune_proximal,endo_sc_harmony_base$orig.ident),margin = 2)->df
df<-as.data.frame(t(df))
df$condition<-reference$condition[match(df$Var1,reference$orig.ident)]
df<-df[which(df$Var2=='inf_endo'),]
t.test(Freq~condition,df[-8,])


ggplot(df[which(df$Var1!='H1'),], aes(x = condition, y = Freq, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
  geom_jitter(width = 0.1, size = 1.5, color = "black") + # raw points
  theme_classic() +
  geom_point(data = df[which(df$Var1=='H1'),],aes(x = condition, y = Freq),color='black',size = 3) +
  labs(title = "Frequency of inflammed Endothelium")->p
p
t.test(Freq~condition,df[which(df$Var1!='H1'),])
ggsave('p.png',p,width = 4,height = 4)

DimPlot(endo_sc_harmony,group.by = 'RNA_snn_res.2',label = T)->p
ggsave('p.png',p,width = 6.5,height = 5,dpi = 300)


# 
FindAllMarkers(endo_sc_harmony,only.pos = T)->marker_all

marker_inf_endo[which(marker_inf_endo$p_val_adj<0.05),]->marker_inf_endo_sig

VlnPlot(endo_sc_harmony,group.by = 'seurat_clusters',features = 'inf_endo_marker1',pt.size = 0)
set.seed(100)
ggsave('p.png',p,width=3,height = 5,dpi = 300)
DimPlot(endo_sc_harmony,label = T)->p
ggsave('p1.png',p,width=5,height = 5,dpi = 300)
aggregate()

endo_sc_harmony_base<-endo_sc_harmony[,grep('1',endo_sc_harmony$orig.ident)]
endo_sc_harmony_base$immune_proximal<-'other'
endo_sc_harmony_base$immune_proximal[which(colnames(endo_sc_harmony_base)%in%colnames(endo_sc_harmony)[which(endo_sc_harmony$seurat_clusters%in%c(1,9,13,15))])]<-'inf_endo'
endo_sc_harmony_base$immune_proximal[which(colnames(endo_sc_harmony_base)%in%colnames(endo_sc_harmony)[which(endo_sc_harmony$seurat_clusters%in%c(11,0,8,7,5,6))])]<-'other_endo'

sc_pro2(endo_sc_harmony,F)->endo_sc_harmony2
AddModuleScore(endo_sc_harmony2,features = list(rownames(marker_inf_endo_sig)[which(marker_inf_endo_sig$avg_log2FC>0)]),name = 'inf_endo_marker')->endo_sc_harmony2
AddModuleScore(endo_sc_harmony2,features = list(rownames(marker_inf_endo_sig)[which(marker_inf_endo_sig$avg_log2FC<0)]),name = 'other_endo_marker')->endo_sc_harmony2
VlnPlot(endo_sc_harmony,group.by = 'seurat_clusters',features = 'inf_endo_marker1',pt.size = 0)
aggregate(endo_sc_harmony2@meta.data$inf_endo_marker1,by=list(endo_sc_harmony2$seurat_clusters),mean)->inf_endo
aggregate(endo_sc_harmony2@meta.data$other_endo_marker1,by=list(endo_sc_harmony2$seurat_clusters),mean)->other_endo


endo_sc_harmony@assays$RNA@data[rownames(marker_inf_endo_sig),which(!endo_sc_harmony$seurat_clusters%in%c(16,12,17))]->mat
endo_sc_harmony@assays$RNA@data[features[1:5],which(!endo_sc_harmony$seurat_clusters%in%c(16,12,17))]->mat
mat_scaled <- t(scale(mat)) # transpose so rows = samples

d <- dist(mat_scaled, method = "euclidean")
hc <- hclust(d, method = "ward.D2")# Plot dendrogram
clusters <- cutree(hc, k = 2)
pheatmap::pheatmap(mat_scaled[order(clusters),rownames(marker_inf_endo_sig)[order(marker_inf_endo_sig$avg_log2FC)]],cluster_row=F,cluster_col=F,scale='column')
apply(mat_scaled[,rownames(marker_inf_endo_sig)[which(marker_inf_endo_sig$avg_log2FC<0)]], 2, function(x) {
  aggregate(x,by=list(clusters),mean)
  t.test(x[which(clusters==1)],x[which(clusters==2)])
}) 


set.seed(123)
km <- kmeans(mat_scaled, centers = 2, nstart = 50)
clusters <- km$cluster
pheatmap::pheatmap(mat_scaled[order(clusters),rownames(marker_inf_endo_sig)[order(marker_inf_endo_sig$avg_log2FC)]],cluster_row=F,cluster_col=F,scale='column',annotation_row=data.frame(row.names=rownames(mat_scaled),cluster=clusters))
pheatmap::pheatmap(mat_scaled[order(clusters),features[1:5]],cluster_row=F,cluster_col=F,scale='column',annotation_row=data.frame(row.names=rownames(mat_scaled),cluster=clusters))
apply(mat_scaled[,features], 2, function(x) {
  aggregate(x,by=list(clusters),mean)
  t.test(x[which(clusters==1)],x[which(clusters==2)])
}) 

sc_pro4<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:20,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:20)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 1,lambda =1.5,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:20,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:20)
  }
  return(object)
}


sc_pro5<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2000)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object,npcs=30)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:20,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 1,lambda =1.5,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:20,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:20,min.dist=0.35)
  }
  return(object)
}
