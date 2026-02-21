load("/ix1/wchen/liutianhao/work_10_05.RData")
DE_sc_fib_select_SSC <- readRDS("~/DE_sc_fib_select_SSC.rds")
DE_sc_fib_select_LS <- readRDS("~/DE_sc_fib_select_LS.rds")
DE_sc_endo_select_LS <- readRDS("~/DE_sc_endo_select_LS.rds")
DE_sc_endo_select_SSC <- readRDS("~/DE_sc_endo_select_SSC.rds")


library(igraph, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(Seurat, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(SeuratWrappers, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(monocle3)
library(dplyr)
library(harmony)
library(ggplot2)

endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_endo_select_LS)),name = 'LS_endo_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_endo_select_SSC)),name = 'SSC_endo_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_fib_select_LS)),name = 'LS_fib_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_fib_select_SSC)),name = 'SSC_fib_marker')

seu<-endo_harmony[,which(endo_harmony$seurat_clusters%in%c(0,1,2))]
seu_ssc<-seu[,which(seu$orig.ident%in%c('E_1','H_1','K_1'))]
sc_pro2(seu_ssc,T)->seu_ssc
FeaturePlot(seu_ssc,c('SSC_fib_marker1','SSC_endo_marker1'))
DimPlot(seu_ssc,group.by = 'seurat_clusters',pt.size = 0.8)

cds <- as.cell_data_set(seu_ssc)
reducedDims(cds)$UMAP <- Embeddings(seu_ssc, "umap")[colnames(cds), ]
reducedDims(cds)$PCA  <- Embeddings(seu_ssc, "harmony")[colnames(cds), 1:30]
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30, resolution = 0.008, random_seed = 1234)
cds <- learn_graph(cds, use_partition = FALSE,learn_graph_control = list(prune_graph = TRUE,
                                                                         minimal_branch_len = 10))
aggregate(seu_ssc$SSC_endo_marker1,by=list(seu_ssc$seurat_clusters),mean)
aggregate(seu_ssc$SSC_fib_marker1,by=list(seu_ssc$seurat_clusters),mean)
root_cells<-colnames(seu_ssc)[which(seu_ssc$seurat_clusters%in%c(2))]
cds <- order_cells(cds, root_cells = root_cells)  # your chosen root
cds->cds_pre
plot_cells(
  cds_pre,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T,
  label_groups_by_cluster = FALSE,          # leaves still shown (optional)
  cell_size = 1,
  label_roots = F,
  label_leaves =F,
  trajectory_graph_color = "#08F007",       # Make lines black
  trajectory_graph_segment_size = 1.5   # Thicker lines
)->p
ggsave('p1.png',p,width = 5,height = 4.5,dpi = 300)
cds->cds_pre

seu_ssc_post<-seu[,which(seu$orig.ident%in%c('E_2','H_2','K_2'))]
sc_pro2(seu_ssc_post,T)->seu_ssc_post
FeaturePlot(seu_ssc_post,c('SSC_fib_marker1','SSC_endo_marker1'))
DimPlot(seu_ssc_post,group.by = 'seurat_clusters',pt.size = 0.8)

cds <- as.cell_data_set(seu_ssc_post)
reducedDims(cds)$UMAP <- Embeddings(seu_ssc_post, "umap")[colnames(cds), ]
reducedDims(cds)$PCA  <- Embeddings(seu_ssc_post, "harmony")[colnames(cds), 1:30]
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30, resolution = 0.008, random_seed = 1234)
cds <- learn_graph(cds, use_partition = FALSE,learn_graph_control = list(prune_graph = TRUE,
                                                                         minimal_branch_len = 10))
aggregate(seu_ssc_post$SSC_endo_marker1,by=list(seu_ssc_post$seurat_clusters),mean)
aggregate(seu_ssc_post$SSC_fib_marker1,by=list(seu_ssc_post$seurat_clusters),mean)
root_cells<-colnames(seu_ssc_post)[which(seu_ssc_post$seurat_clusters%in%c(2))]
cds <- order_cells(cds, root_cells = root_cells)  # your chosen root
cds->cds_post
plot_cells(
  cds_post,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T,
  label_groups_by_cluster = FALSE,
  label_leaves = F,
  label_branch_points = F,
  label_root = F,
  cell_size = 1,
  trajectory_graph_color = "cyan",       # Make lines black
  trajectory_graph_segment_size = 1   # Thicker lines
  
) + labs(x = "Component 1", y = "Component 2")->p2
ggsave('p2.png',p2,width = 6,height = 5.5,dpi = 300)
cds->cds_post


cds_pre@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->seu_ssc$pseudotime
cds_post@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->seu_ssc_post$pseudotime
seu_ssc<-AddModuleScore(seu_ssc,features = list(rownames(DE_sc_endo_select_LS)),name = 'LS_endo_marker')
seu_ssc<-AddModuleScore(seu_ssc,features = list(rownames(DE_sc_endo_select_SSC)),name = 'SSC_endo_marker')
seu_ssc<-AddModuleScore(seu_ssc,features = list(rownames(DE_sc_fib_select_LS)),name = 'LS_fib_marker')
seu_ssc<-AddModuleScore(seu_ssc,features = list(rownames(DE_sc_fib_select_SSC)),name = 'SSC_fib_marker')
cor.test(log10(seu_ssc@meta.data$min_dis_to_fib/4),seu_ssc$pseudotime)
cor.test(seu_ssc_post@meta.data$min_dis_to_fib,seu_ssc_post$pseudotime)
cor.test(seu_ssc@meta.data$SSC_endo_marker1,seu_ssc$pseudotime)
cor.test(seu_ssc_post@meta.data$SSC_endo_marker1,seu_ssc_post$pseudotime)
cor.test(seu_ssc@meta.data$SSC_fib_marker1,seu_ssc$pseudotime)
cor.test(seu_ssc_post@meta.data$SSC_fib_marker1,seu_ssc_post$pseudotime)
pheatmap::pheatmap(data.frame(endo=c(0.7196631,0.2335714),fib=c(0.6665984,0.1112042)),cluster_rows = FALSE,cluster_cols = FALSE,border_color = NA)





sc_pro2<-function(object,harmony){
  NormalizeData(object)->object
  VariableFeatures(object)<-union(rownames(DE_sc_endo_select_SSC),rownames(DE_sc_fib_select_SSC))
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 1.2, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda = 4,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 1, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}

