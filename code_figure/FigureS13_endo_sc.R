load("/ix1/wchen/liutianhao/work_rdata/work_11_27.RData")
library(Seurat)
library(harmony)
library(ggplot2)
#seurat4.4.0
DimPlot(endo_sc_harmony,label = T)
FindMarkers(endo_sc_harmony,ident.1 = c(4,1,5,9),ident.2 = c(0,6,7,12,8,10,11),only.pos = F)->marker_endo_sc_harmony
endo_sc_harmony_select<-endo_sc_harmony[,which(endo_sc_harmony$seurat_clusters%in%c(4,1,5,9))]
endo_sc_harmony_select$seurat_clusters_old<-endo_sc_harmony$seurat_clusters[match(colnames(endo_sc_harmony_select),colnames(endo_sc_harmony))]
#set lambda as 1.5, other parameters are default
sc_pro(endo_sc_harmony_select,T,lambda = 4)->endo_sc_harmony_select
DimPlot(endo_sc_harmony_select,group.by = 'seurat_clusters',label = T)
DimPlot(endo_sc_harmony_select[,which(endo_sc_harmony_select$seurat_clusters!=7)],group.by = 'seurat_clusters')
sc_pro(endo_sc_harmony_select[,which(endo_sc_harmony_select$seurat_clusters!=7)],T,lambda = 1.5)->endo_sc_harmony_select
DimPlot(endo_sc_harmony_select,group.by = 'seurat_clusters_old',pt.size = 0.3)
DimPlot(endo_sc_harmony_select,group.by = 'seurat_clusters',pt.size = 0.3,label = T)

endo_sc_harmony_select$seurat_clusters_old[which(endo_sc_harmony_select$seurat_clusters_old==9)]<-5
FindAllMarkers(endo_sc_harmony_select,group.by = 'seurat_clusters_old',only.pos = T)->marker_endo_sc_harmony_select
features_sc<-marker_endo_sc_harmony_select%>%group_by(cluster)%>%slice_head(n=5)
FindAllMarkers(seu_ssc,group.by = 'seurat_clusters',only.pos = T,features = intersect(unique(marker_endo_sc_harmony_select$gene),rownames(seu_ssc)))->marker_endo_st_harmony_select
features_st<-marker_endo_st_harmony_select%>%group_by(cluster)%>%slice_head(n=5)

#marker_endo_sc_harmony_select[which(marker_endo_sc_harmony_select$gene=='ACTA2'),]
#intersect(marker_endo_cap$gene[which(marker_endo_cap$cluster==1)],marker_endo_sc_harmony_select$gene[which(marker_endo_sc_harmony_select$cluster==4)])
#intersect(marker_endo_cap$gene[which(marker_endo_cap$cluster==2)],marker_endo_sc_harmony_select$gene[which(marker_endo_sc_harmony_select$cluster==5)])
features_sc$gene[11:15]<-c('COL4A1','TIMP3','IGFBP7','FN1','SPARC')
features_sc$gene[7:8]<-c('ECSCR','MRTFB')
endo_sc_harmony_select$seurat_clusters_old<-factor(endo_sc_harmony_select$seurat_clusters_old,levels = c(4,1,5))
DotPlot(endo_sc_harmony_select,features = c(features_sc$gene[6:10],'S100A6',features_sc$gene[1:5],features_sc$gene[11:15]),group.by = 'seurat_clusters_old')+coord_flip()->p
ggsave('p.png',p,width = 6,height = 6,dpi = 300)
#check the markers in spatial data
DotPlot(seu_ssc,features = c(features_sc$gene[6:10],'S100A6',features_sc$gene[1:5],features_sc$gene[11:15]),group.by = 'seurat_clusters')+coord_flip()->p


#endo_sc_harmony_select3<-endo_sc_harmony_select
#FindClusters(endo_sc_harmony_select3,resolution = 2)->endo_sc_harmony_select3
#chisq.test(table(endo_sc_harmony_select3$seurat_clusters,endo_sc_harmony_select3$seurat_clusters_old))->test
#colnames(as.matrix( test$stdres))[apply(as.matrix( test$stdres), 1,function(x) which.max(x))]->array
#names(array)<-as.character(0:12)

#1-4(LS, CD74, S100A6) {red}, 0-1(intermediate AQP1) {green}, 2-5(col4A1) {blue}
#chisq.test(table(seu_ssc_all$seurat_clusters,seu_ssc_all$condition))$stdres
#chisq.test(table(endo_sc_harmony_select3$seurat_clusters_old,endo_sc_harmony_select3$condition))$stdres
#chisq.test(table(endo_sc_harmony_select3$seurat_clusters,endo_sc_harmony_select3$condition))$stdres->mtx
#ifelse(endo_sc_harmony_select3$seurat_clusters%in%rownames(mtx)[which(mtx[,1]>2)],'SSC_enriched','other')->endo_sc_harmony_select3$enrich
#endo_sc_harmony_select3$seurat_clusters_new<-factor(endo_sc_harmony_select3$seurat_clusters_new,levels = c(4,1,5))
#DimPlot(endo_sc_harmony_select3,group.by = 'seurat_clusters_new')->p
#ggsave('p.png',p,width = 6.5,height = 5.5,dpi = 300)

#endo_sc_harmony_select$seurat_clusters_old<-endo_sc_harmony$seurat_clusters[match(colnames(endo_sc_harmony_select),colnames(endo_sc_harmony))]

endo_sc_harmony_select<-AddModuleScore(endo_sc_harmony_select,features = list(rownames(DE_sc_endo_select_LS)),name = 'LS_endo_marker')
endo_sc_harmony_select<-AddModuleScore(endo_sc_harmony_select,features = list(rownames(DE_sc_endo_select_SSC)),name = 'SSC_endo_marker')
FeaturePlot(endo_sc_harmony_select,'LS_endo_marker1')
FeaturePlot(endo_sc_harmony_select,'SSC_endo_marker1')



#enrichment test
endo_sc_harmony_select2<-endo_sc_harmony_select[,which(endo_sc_harmony_select$condition!='HC')]
chisq.test(endo_sc_harmony_select2$condition,endo_sc_harmony_select2$seurat_clusters_old)->test
test$stdres

#enrichment between immune proximal cells and three states of vascular endothelial cells
merged_filter_exclude_epi_endo <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_endo.rds")
merged_filter_exclude_epi_endo<-merged_filter_exclude_epi_endo[,which(merged_filter_exclude_epi_endo$cell_type_l1_all%in%c('Endothelium','inf_endo'))]
seu_ssc$inf_cluster<-merged_filter_exclude_epi_endo$cell_type_l1_all[match(colnames(seu_ssc),colnames(merged_filter_exclude_epi_endo))]
chisq.test(seu_ssc$inf_cluster,seu_ssc$seurat_clusters)->test
2 * pnorm(-abs(test$stdres))
pheatmap::pheatmap(test$stdres,cluster_rows = F,cluster_cols = F)

chisq.test(seu_ssc$condition,seu_ssc$seurat_clusters)->test
2 * pnorm(-abs(test$stdres))
pheatmap::pheatmap(test$stdres,cluster_rows = F,cluster_cols = F)

run_monocle_sc<-function(condition) {
  cds <- as.cell_data_set(endo_sc_harmony_select[,which(endo_sc_harmony_select$condition==condition)])
  reducedDims(cds)$UMAP <- Embeddings(endo_sc_harmony_select, "umap")[colnames(cds), ]
  reducedDims(cds)$PCA  <- Embeddings(endo_sc_harmony_select, "harmony")[colnames(cds), 1:30]
  cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30, resolution = 0.008, random_seed = 1234)
  cds <- learn_graph(cds, use_partition = FALSE,learn_graph_control = list(prune_graph = TRUE,minimal_branch_len = 10))
  aggregate(endo_sc_harmony_select$SSC_endo_marker1,by=list(endo_sc_harmony_select$seurat_clusters_old),mean)
  DimPlot(endo_sc_harmony_select,group.by = 'seurat_clusters_old')
  root_cells<-colnames(cds)[which(cds@colData$seurat_clusters_old==c(1,4))]
  cds <- order_cells(cds, root_cells = root_cells)  # your chosen root
  plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = T,
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    cell_size = 1
  ) + labs(x = "Component 1", y = "Component 2")->p
  return(cds)
}
run_monocle_sc('SSC')->cds_ssc
run_monocle_sc('LS')->cds_ls
endo_sc_harmony_select$pseudotime<-1
cds_ssc@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->endo_sc_harmony_select$pseudotime[which(endo_sc_harmony_select$condition=='SSC')]
cds_ls@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->endo_sc_harmony_select$pseudotime[which(endo_sc_harmony_select$condition=='LS')]
FeaturePlot(endo_sc_harmony_select[,which(endo_sc_harmony_select$condition=='SSC')],'pseudotime')
FeaturePlot(endo_sc_harmony_select[,which(endo_sc_harmony_select$condition=='LS')],'pseudotime')



sc_pro<-function(object,harmony,lambda){
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
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda = lambda,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}
