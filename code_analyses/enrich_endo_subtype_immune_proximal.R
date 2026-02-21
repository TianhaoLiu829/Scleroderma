merged_filter_exclude_epi_endo_fib <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_endo_fib.rds")
seu_ssc_all$min_dis_to_B_Plasma<-merged_filter_exclude_epi$min_dis_to_B_Plasma[match(colnames(seu_ssc_all),colnames(merged_filter_exclude_epi))]
seu_ssc_all$min_dis_to_immune<-apply(seu_ssc_all@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')],1,min)
summary(aov(min_dis_to_immune~seurat_clusters,seu_ssc_all@meta.data))
summary(aov(min_dis_to_Macrophage~seurat_clusters,seu_ssc_all@meta.data))
summary(aov(min_dis_to_fib~seurat_clusters,seu_ssc_all@meta.data))

seu_ssc_sep$min_dis_to_B_Plasma<-merged_filter_exclude_epi$min_dis_to_B_Plasma[match(colnames(seu_ssc_sep),colnames(merged_filter_exclude_epi))]
seu_ssc_sep$min_dis_to_immune<-apply(seu_ssc_sep@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')],1,min)
summary(aov(min_dis_to_immune~seurat_clusters,seu_ssc_sep@meta.data))
summary(aov(min_dis_to_Macrophage~seurat_clusters,seu_ssc_sep@meta.data))
summary(aov(min_dis_to_fib~seurat_clusters,seu_ssc_sep@meta.data))

chisq.test(table(seu_ssc_all$seurat_clusters,seu_ssc_all$cell_type_l1_all))
seu_ssc_all$cell_type_l1_all<-merged_filter_exclude_epi_endo_fib$cell_type_l1_all[match(colnames(seu_ssc_all),colnames(merged_filter_exclude_epi_endo_fib))]
chisq.test(table(seu_ssc_all$seurat_clusters,seu_ssc_all$cell_type_l1_all))

chisq.test(table(seu_ssc_all$seurat_clusters[which(seu_ssc_all$condition=='SSC')],seu_ssc_all$cell_type_l1_all[which(seu_ssc_all$condition=='SSC')]))
chisq.test(table(seu_ssc_all$seurat_clusters[which(seu_ssc_all$condition=='LS')],seu_ssc_all$cell_type_l1_all[which(seu_ssc_all$condition=='LS')]))

chisq.test(table(seu_ssc_all$seurat_clusters,seu_ssc_all$cell_type_l1_all))->test
pvals <- 2 * pnorm(abs(test$stdres), lower.tail = FALSE)
pheatmap::pheatmap(as.matrix(test$stdres)[c(2,1,3),c(2,1)],cluster_rows = F,cluster_cols = F)->p
ggsave('p.png',p,width = 3,height = 4,dpi = 300)

seu_ssc_all$seurat_clusters<-factor(seu_ssc_all$seurat_clusters,levels = c(1,0,2))
chisq.test(table(seu_ssc_all$seurat_clusters,seu_ssc_all$condition))->test
pheatmap::pheatmap(test$stdres,cluster_rows = F,cluster_cols = F)->p
ggsave('p.png',p,width = 4,height = 5.3,dpi = 300)
chisq.test(table(seu_ssc_all$cell_type_l1_all,seu_ssc_all$condition))->test
seu_ssc_all$inf_endo<-colnames(seu_ssc_all)%in%colnames(merged_filter_exclude_epi_endo_fib)[which(merged_filter_exclude_epi_endo_fib$cell_type_l1_all=='inf_endo')]
chisq.test(table(seu_ssc_all$inf_endo,seu_ssc_all$condition))->test




