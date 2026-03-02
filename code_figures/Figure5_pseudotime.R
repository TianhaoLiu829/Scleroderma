load("/ix1/wchen/liutianhao/work_rdata/work_10_05.RData")
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")

library(igraph, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(Seurat, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(SeuratWrappers, lib.loc = "/software/rhel9/manual/install/r/4.5.0/lib64/R/library")
library(monocle3)
library(dplyr)
library(harmony)
library(ggplot2)

endo_marker_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/endo_marker_LS.rds")
DE_sc_endo_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_endo_inter_SSC.rds")
DE_sc_fib_inter_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_LS.rds")
DE_sc_fib_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_SSC.rds")

endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(endo_marker_LS)),name = 'LS_endo_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_endo_inter_SSC)),name = 'SSC_endo_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_fib_inter_LS)),name = 'LS_fib_marker')
endo_harmony<-AddModuleScore(endo_harmony,features = list(rownames(DE_sc_fib_inter_SSC)),name = 'SSC_fib_marker')


DimPlot(endo_harmony,group.by = 'orig.ident')
DimPlot(endo_harmony,group.by = 'seurat_clusters',label = T)
FeaturePlot(endo_harmony,'SSC_fib_marker1')
FeaturePlot(endo_harmony,'SSC_endo_marker1')

#cor.test(endo_harmony$min_dis_to_fib[which(endo_harmony$orig.ident%in%c('E_2','H_2','K_2'))],endo_harmony$SSC_endo_marker1[which(endo_harmony$orig.ident%in%c('E_2','H_2','K_2'))])
seu<-endo_harmony[,which(endo_harmony$seurat_clusters%in%c(0,1,2))]
seu<-seu[,grep('1',seu$orig.ident)]
seu_ssc<-seu
sc_pro(seu_ssc,F)->seu_ssc
FeaturePlot(seu_ssc,c('SSC_fib_marker1','SSC_endo_marker1'))
DimPlot(seu_ssc,group.by = 'seurat_clusters',pt.size = 0.8)->p
ggsave('p.png',p,width = 6,height = 3.5,dpi = 300)

#Figure 5D
seu_ssc$cluster2<-paste0(ifelse(seu_ssc$seurat_clusters==2,'Cluster2','other'),'_',seu_ssc$condition)
seu_ssc$cluster2<-factor(seu_ssc$cluster2,levels = rev(c('Cluster2_SSC','Cluster2_LS','other_SSC','other_LS')))
dot_plot(seu_ssc,features = rownames(DE_sc_endo_inter_SSC)[c(1:15,16,17,21,23,27,38)],"cluster2")+coord_flip()+theme(axis.text.x = element_text(size = 11))->p
ggsave('p.png',p,width = 8,height = 3.6,dpi = 300)

#Figure 5E
DimPlot(seu_ssc,group.by = 'condition')
set.seed(100)
FindClusters(seu_ssc,resolution=0.5)->seu_ssc
seu_ssc$seurat_clusters<-factor(seu_ssc$seurat_clusters,levels = c(1,0,2))
VlnPlot(seu_ssc[,which(seu_ssc$condition=='SSC')],'SSC_endo_marker1',group.by = 'seurat_clusters',pt.size = 0)


#pseudo time analyses
#combine in pseudotime run
run_monocle<-function(cds,min_cluster){
  reducedDims(cds)$UMAP <- Embeddings(seu_ssc, "umap")[colnames(cds), ]
  reducedDims(cds)$PCA  <- Embeddings(seu_ssc, "harmony")[colnames(cds), 1:30]
  cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30, resolution = 0.008, random_seed = 1234)
  cds <- learn_graph(cds, use_partition = FALSE,learn_graph_control = list(prune_graph = TRUE,minimal_branch_len = 10))
  #aggregate(seu_ssc$SSC_fib_marker1,by=list(seu_ssc$seurat_clusters),mean)
  DimPlot(seu_ssc,group.by = 'seurat_clusters')
  root_cells<-colnames(cds)[which(cds@colData$seurat_clusters==min_cluster)]
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
  return(list(p,cds))
}

aggregate(seu_ssc$SSC_endo_marker1,by=list(seu_ssc$seurat_clusters),mean)->a
a$Group.1[which.min(a$x)]->min_cluster

cds <- as.cell_data_set(seu_ssc[,which(seu_ssc$condition=='SSC')])
run_monocle(cds,min_cluster)->object
object[[2]]->cds_ssc
object[[1]]
cds <- as.cell_data_set(seu_ssc[,which(seu_ssc$condition=='LS')])
run_monocle(cds,min_cluster)->object
object[[2]]->cds_ls
object[[1]]

seu_ssc$pseudotime<-1
cds_ls@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->seu_ssc$pseudotime[which(seu_ssc$condition=='LS')]
cds_ssc@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->seu_ssc$pseudotime[which(seu_ssc$condition=='SSC')]
seu_ssc->seu_ssc_all


cor.test(seu_ssc_all$SSC_endo_marker1[which(seu_ssc$condition=='LS')],seu_ssc$pseudotime[which(seu_ssc$condition=='LS')])
cor.test(seu_ssc_all$SSC_endo_marker1[which(seu_ssc$condition=='SSC')],seu_ssc$pseudotime[which(seu_ssc$condition=='SSC')])
cor.test(seu_ssc_all$SSC_fib_marker1[which(seu_ssc$condition=='LS')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='LS')])
cor.test(seu_ssc_all$SSC_fib_marker1[which(seu_ssc$condition=='SSC')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='SSC')])

cor.test(seu_ssc_all$LS_fib_marker1[which(seu_ssc$condition=='LS')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='LS')])
cor.test(seu_ssc_all$LS_fib_marker1[which(seu_ssc$condition=='SSC')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='SSC')])
cor.test(seu_ssc_all$LS_endo_marker1[which(seu_ssc$condition=='LS')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='LS')])
cor.test(seu_ssc_all$LS_endo_marker1[which(seu_ssc$condition=='SSC')],seu_ssc_all$pseudotime[which(seu_ssc_all$condition=='SSC')])

#Figure 5F
#test the individual object
seu<-endo_harmony[,which(endo_harmony$seurat_clusters%in%c(0,1,2))]
seu<-seu[,grep('1',seu$orig.ident)]

seperate_monocle<-function(condition,features) {
  seu_condition<-seu[,which(seu$condition==condition)]
  sc_pro2(seu_condition,T,features)->seu_condition
  cds <- as.cell_data_set(seu_condition)
  reducedDims(cds)$UMAP <- Embeddings(seu_condition, "umap")[colnames(cds), ]
  reducedDims(cds)$PCA  <- Embeddings(seu_condition, "harmony")[colnames(cds), 1:30]
  cds <- cluster_cells(cds, reduction_method = "UMAP", k = 30, resolution = 0.008, random_seed = 1234)
  cds <- learn_graph(cds, use_partition = FALSE,learn_graph_control = list(prune_graph = TRUE,minimal_branch_len = 10))
  aggregate(seu_condition$SSC_endo_marker1,by=list(seu_condition$seurat_clusters),mean)->a
  a$Group.1[which.min(a$x)]->min_cluster
  root_cells<-colnames(cds)[which(cds@colData$seurat_clusters==min_cluster)]
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
  seu_condition$pseudotime<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  return(list(p,cds,seu_condition))
}

DE_sc_fib_select_SSC <- readRDS("~/DE_sc_fib_select_SSC.rds")
DE_sc_endo_select_SSC <- readRDS("~/DE_sc_endo_select_SSC.rds")
#DE_sc_fib_select_LS <- readRDS("~/DE_sc_fib_select_LS.rds")
#DE_sc_endo_select_LS <- readRDS("~/DE_sc_endo_select_LS.rds")

union(rownames(DE_sc_endo_select_SSC),rownames(DE_sc_fib_select_SSC))->features_ssc
seperate_monocle('SSC',features_ssc)->object_ssc
object_ssc[[3]]->seu_ssc_sep
object_ssc[[1]]
DimPlot(seu_ssc_sep,group.by = 'seurat_clusters')
seu_ssc_sep$seurat_clusters<-factor(seu_ssc_sep$seurat_clusters,levels = c(2,0,1))
DimPlot(seu_ssc_sep,group.by = 'seurat_clusters')->p
ggsave('p.png',p,width = 6,height = 5,dpi = 300)
seu_ssc_sep$cluster_old<-seu_ssc$seurat_clusters[match(colnames(seu_ssc_sep),colnames(seu_ssc))]
DimPlot(seu_ssc_sep,group.by = 'cluster_old')->p
ggsave('p.png',p,width = 6,height = 5,dpi = 300)


#union(rownames(DE_sc_endo_select_LS),rownames(DE_sc_fib_select_LS))->features_ls
seperate_monocle('LS',features_ssc)->object_ls
object_ls[[3]]->seu_ls_sep
object_ls[[1]]
seu_ls_sep$cluster_old<-seu_ssc$seurat_clusters[match(colnames(seu_ls_sep),colnames(seu_ssc))]
DimPlot(seu_ls_sep,group.by = 'cluster_old')->p
ggsave('p.png',p,width = 6,height = 5,dpi = 300)


cor.test(seu_ssc_sep$pseudotime,seu_ssc_sep$SSC_endo_marker1,method = 'spearman')->a_ssc_sep_endo
cor.test(seu_ls_sep$pseudotime,seu_ls_sep$SSC_endo_marker1,method = 'spearman')->a_ls_sep_endo
cor.test(seu_ssc_sep$pseudotime,seu_ssc_sep$SSC_fib_marker1,method = 'spearman')->a_ssc_sep_fib
cor.test(seu_ls_sep$pseudotime,seu_ls_sep$SSC_fib_marker1,method = 'spearman')->a_ls_sep_fib

a_ssc_sep_endo$p.value
a_ls_sep_endo$p.value
a_ssc_sep_fib$p.value
a_ls_sep_fib$p.value

a_ssc_sep_endo$estimate
a_ls_sep_endo$estimate
a_ssc_sep_fib$estimate
a_ls_sep_fib$estimate


library(ggpubr)
ggplot(seu_ssc_sep@meta.data, aes(y = SSC_endo_marker1, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  #stat_cor(method = "pearson", label.x = 0.5, label.y = max(seu_ssc_sep@meta.data$SSC_endo_marker1)*0.9) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p
ggplot(seu_ssc_sep@meta.data, aes(y = SSC_fib_marker1, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  #stat_cor(method = "pearson", label.x = 0.5, label.y = max(seu_ssc_sep@meta.data$SSC_fib_marker1)*0.9) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p

ggplot(seu_ls_sep@meta.data, aes(y = SSC_endo_marker1, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  #stat_cor(method = "pearson", label.x = 0.5, label.y = max(seu_ls_sep@meta.data$SSC_endo_marker1)*0.9) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p
ggplot(seu_ls_sep@meta.data, aes(y = SSC_fib_marker1, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  #stat_cor(method = "pearson", label.x = 0.5, label.y = max(seu_ls_sep@meta.data$SSC_fib_marker1)*0.9) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p

#plot pseudotime~cluster (Figure S13D)
DimPlot(seu_ssc_sep)
seu_ssc_sep$seurat_clusters<-factor(seu_ssc_sep$seurat_clusters,levels = c(2,0,1))
#VlnPlot(seu_ssc_sep,features = 'pseudotime',group.by = 'cluster_old',pt.size = 0)+coord_flip()->p
VlnPlot(seu_ssc_sep,features = 'pseudotime',group.by = 'seurat_clusters',pt.size = 0)+coord_flip()->p
ggsave('p.png',p,width = 5,height = 5,dpi = 300)
#pheatmap::pheatmap(data.frame(endo=c(0.7196631,0.2335714),fib=c(0.6665984,0.1112042)),cluster_rows = FALSE,cluster_cols = FALSE,border_color = NA)


#distance to other immune cells 
test_dist<-function(cell_type){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  lapply(dist_matrix_all, function(x) {
    as.data.frame(apply(x[,colnames(x)%in%names],1,min))
  })->all
  do.call(rbind,all)->all
  return(all[match(colnames(merged_filter_exclude_epi),gsub('.*[.]','',rownames(all))),1])
}
'endo_fib'->merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(seu_ssc)[which(seu_ssc$seurat_clusters==2)])]
test_dist('endo_fib')->merged_filter_exclude_epi$min_dis_to_endo_fib
seu_ssc_all$min_dis_to_endo_fib<-merged_filter_exclude_epi$min_dis_to_endo_fib[match(colnames(seu_ssc_all),colnames(merged_filter_exclude_epi))]
seu_ssc_sep$min_dis_to_endo_fib<-merged_filter_exclude_epi$min_dis_to_endo_fib[match(colnames(seu_ssc_sep),colnames(merged_filter_exclude_epi))]
seu_ls_sep$min_dis_to_endo_fib<-merged_filter_exclude_epi$min_dis_to_endo_fib[match(colnames(seu_ls_sep),colnames(merged_filter_exclude_epi))]

#distance~pseudotime imputed in all LS and SSc
seu_ssc_all2<-seu_ssc_all[,which(!seu_ssc_all$orig.ident%in%c('I_1','J_1'))]
cor.test(seu_ssc_all2$min_dis_to_endo_fib,seu_ssc_all2$pseudotime,method = 'spearman')
seu_ssc_all3<-seu_ssc_all2[,which(!seu_ssc_all2$seurat_clusters%in%c(2))]
cor.test(seu_ssc_all3$min_dis_to_endo_fib,seu_ssc_all3$pseudotime,method = 'spearman')

#distance~pseudotime imputed in LS and SSc separately
seu_ls_sep2<-seu_ls_sep[,which(!seu_ls_sep$orig.ident%in%c('I_1','J_1'))]
cor.test(seu_ls_sep2$min_dis_to_endo_fib,seu_ls_sep2$pseudotime,method = 'spearman')
seu_ls_sep3<-seu_ls_sep2[,which(seu_ls_sep2$min_dis_to_endo_fib!=0)]
cor.test(seu_ls_sep3$min_dis_to_endo_fib,seu_ls_sep3$pseudotime)

cor(seu_ssc_sep$min_dis_to_endo_fib,seu_ssc_sep$pseudotime)
seu_ssc_sep2<-seu_ssc_sep[,which(seu_ssc_sep$min_dis_to_endo_fib!=0)]
cor.test(seu_ssc_sep2$min_dis_to_endo_fib,seu_ssc_sep2$pseudotime)

#compute the density
df_endo<-as.data.frame(seu_ssc_sep@meta.data)
dens <- density(seu_ssc_sep$pseudotime[which(seu_ssc_sep$min_dis_to_endo_fib!=0)], bw = "nrd0",from = min(seu_ssc_sep$pseudotime), to   = max(seu_ssc_sep$pseudotime))
density_df <- data.frame(
  pseudotime = dens$x,
  density = dens$y
)
density_df$density_scaled <- density_df$density / max(density_df$density)

#plot out the figure 5 correlations between distance to fibrotic endothelial cells and pseudo times

library(dplyr)
df<-seu_ssc_sep@meta.data
df<-df[which(!is.infinite(df$min_dis_to_endo_fib)),]
#df<-df[which(!is.na(df$min_dis_to_endo_fib)&df$min_dis_to_endo_fib<400),]
df$mature_numeric <- ifelse(df$min_dis_to_endo_fib == 0, 1, 0)
df$min_dis_to_endo_fib<-df$min_dis_to_endo_fib/4
library(mgcv)
gam_fit <- gam(mature_numeric ~ s(pseudotime, k = 50), 
               data = df, 
               family = binomial)
pt_grid <- seq(min(df$pseudotime), max(df$pseudotime), length.out = 2000)
pred_df <- data.frame(
  pseudotime = pt_grid,
  frac_mature = predict(gam_fit, newdata = data.frame(pseudotime = pt_grid), type = "response")
)
strip_height <- 0.05 * max(df$min_dis_to_endo_fib)
tile_width <- (max(df$pseudotime, na.rm = TRUE) -
                 min(df$pseudotime, na.rm = TRUE)) / 2000


ggplot(df, aes(y = min_dis_to_endo_fib, x = pseudotime)) +
  geom_point(alpha = 0.5, size = 2, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  theme_classic(base_size = 14) +
  geom_tile(
    data = pred_df,
    aes(
      x = pseudotime,
      y= -strip_height,
      fill = frac_mature
    ),
    height = strip_height,
    width = (max(df$pseudotime) - min(df$pseudotime)) / 2000   # matches the grid spacing
  ) +
  scale_y_reverse(
    name = "Left axis (High → Low)",
    sec.axis = sec_axis(~ max(df$min_dis_to_endo_fib) + min(df$min_dis_to_endo_fib) - ., name = "Right axis (Low → High)")
  ) +
  #stat_cor(method = "pearson", label.x = 0.5)+
  scale_fill_gradient(limits = c(0, 1),low = "white", high = "red", name = "% mature (smooth)")->p
ggsave('p.png',p,width = 8.5,height = 6,dpi = 300)
cor.test(df$min_dis_to_endo_fib,df$pseudotime)->test


#pre-post treatment
seu<-endo_harmony[,which(endo_harmony$seurat_clusters%in%c(0,1,2))]
seu_ssc_post<-seu[,which(seu$orig.ident%in%c('E_2','H_2','K_2'))]

sc_pro2(seu_ssc_post,T,features = union(rownames(DE_sc_endo_select_SSC),rownames(DE_sc_fib_select_SSC)))->seu_ssc_post
#sc_pro(seu_ssc_post,T)->seu_ssc_post
FindClusters(seu_ssc_post,resolution = 1)->seu_ssc_post

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

cds_post@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]->seu_ssc_post$pseudotime

cor.test(seu_ssc_sep@meta.data$SSC_endo_marker1,seu_ssc_sep$pseudotime,method = 'spearman')->a_ssc_sep_endo_pre
cor.test(seu_ssc_post@meta.data$SSC_endo_marker1[which(seu_ssc_post$seurat_clusters!=2)],seu_ssc_post$pseudotime[which(seu_ssc_post$seurat_clusters!=2)],method = 'spearman')->a_ssc_sep_endo_post
cor.test(seu_ssc_sep@meta.data$SSC_fib_marker1,seu_ssc_sep$pseudotime,method = 'spearman')->a_ssc_sep_fib_pre
cor.test(seu_ssc_post@meta.data$SSC_fib_marker1[which(seu_ssc_post$seurat_clusters!=2)],seu_ssc_post$pseudotime[which(seu_ssc_post$seurat_clusters!=2)],method = 'spearman')->a_ssc_sep_fib_post
df<-data.frame(endo=c(a_ssc_sep_endo_pre$estimate,a_ssc_sep_endo_post$estimate),fib=c(a_ssc_sep_fib_pre$estimate,a_ssc_sep_fib_post$estimate))
pheatmap::pheatmap(df,cluster_rows = FALSE,cluster_cols = FALSE,border_color = NA,breaks = seq(-max(abs(df))-0.02,max(abs(df))+0.02,length.out=101))->p
ggsave('p_pre_post_gene.png',p,width = 5,height = 4,dpi = 300)
df_p<-data.frame(endo=c(a_ssc_sep_endo_pre$p.value,a_ssc_sep_endo_post$p.value),fib=c(a_ssc_sep_fib_pre$p.value,a_ssc_sep_fib_post$p.value))

#calculate the distance to endo_fib post treatment
seu_ssc_post<- FindClusters(seu_ssc_post,resolution = 1, verbose = FALSE,random.seed = 42)
aggregate(seu_ssc_post$SSC_endo_marker1,by=list(seu_ssc_post$seurat_clusters),mean)
aggregate(seu_ssc_post$SSC_fib_marker1,by=list(seu_ssc_post$seurat_clusters),mean)
DimPlot(seu_ssc_post)
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
'endo_fib'->merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi )%in%colnames(seu_ssc_post)[which(seu_ssc_post$seurat_clusters==0)])]
test_dist('endo_fib')->merged_filter_exclude_epi$min_dis_to_endo_fib
seu_ssc_post$min_dis_to_endo_fib<-merged_filter_exclude_epi$min_dis_to_endo_fib[match(colnames(seu_ssc_post),colnames(merged_filter_exclude_epi))]
cor.test(seu_ssc_post@meta.data$min_dis_to_endo_fib[which(seu_ssc_post$seurat_clusters!=2)],seu_ssc_post$pseudotime[which(seu_ssc_post$seurat_clusters!=2)],method = 'spearman')->a_ssc_sep_dis_post
cor.test(seu_ssc_sep@meta.data$min_dis_to_endo_fib,seu_ssc_sep$pseudotime,method = 'spearman')->a_ssc_sep_dis_pre
library(RColorBrewer)
pheatmap::pheatmap(data.frame(dis=c(a_ssc_sep_dis_pre$estimate,a_ssc_sep_dis_post$estimate)),cluster_rows = FALSE,cluster_cols = FALSE,border_color = NA,color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),breaks = seq(-abs(max(a_ssc_sep_dis_pre$estimate))+0.03,abs(max(a_ssc_sep_dis_pre$estimate))-0.03,length.out=101))->p
ggsave('p_pre_post.png',p,width = 2.7,height = 4,dpi = 300)
df_p_dis<-data.frame(dis=c(a_ssc_sep_dis_pre$p.value,a_ssc_sep_dis_post$p.value))



sc_pro2<-function(object,harmony,features){
  NormalizeData(object)->object
  VariableFeatures(object)<-features
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

