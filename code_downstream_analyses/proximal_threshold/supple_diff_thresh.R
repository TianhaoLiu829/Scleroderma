merged_filter_exclude_epi2 <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi2->merged_filter_exclude_epi
merged_filter_exclude_epi_fib<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts')]
fib_marker<-c('CD74','C3','CCL19','B2M','CXCL9','COL1A1','COL3A1','COL1A2','SPARC','DSP')
res_fib_list<-list()
for (thresh in 80:240) {
  merged_filter_exclude_epi_fib$inf_fib<-FALSE
  merged_filter_exclude_epi_fib$inf_fib[which(apply(merged_filter_exclude_epi_fib@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')], 1, function(x) {sum(x<thresh)>0}))]<-T
  merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_fib)[which(merged_filter_exclude_epi_fib$inf_fib)])]<-'inf_fib'
  res_fib_list[[as.character(thresh)]]<-DE_wilcox(as.matrix(merged_filter_exclude_epi_fib@assays$RNA@data[fib_marker,]),merged_filter_exclude_epi_fib$inf_fib,TRUE,FALSE)
}
res_fib_list2<-lapply(fib_marker, function(gene){
  as.data.frame(do.call(rbind,lapply(res_fib_list, function(x) {x[gene,]})))->df
  rownames(df)->df$thresh
  df
})
names(res_fib_list2)<-fib_marker

lapply(res_fib_list2, function(df) {
  maxY<-max(abs(df$log2FC))
  as.numeric(df$thresh)->df$thresh
  ggplot(df, aes(x = thresh/4, y = log2FC, color = -log10(pvalue))) +
    geom_point(size = 1.5) +
    
    # draw the x-axis at y=0
    geom_hline(yintercept = 0, color='black',linetype='dashed', linewidth = 0.6) +
    scale_color_viridis_c() +
    # symmetric y-axis
    scale_y_continuous(limits = c(-maxY, maxY)) +
    geom_point(data = df['132', ], aes(x = thresh/4, y = log2FC),
               color = "orange", size = 1.4) +
    theme_classic()->p
  ggsave(paste0('~/test_diff_thresh/fib/p_',unique(df$gene),'.png'),p,width = 5,height = 4,dpi = 300)
})



merged_filter_exclude_epi2->merged_filter_exclude_epi
merged_filter_exclude_epi_endo<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')]
endo_marker<-c("CD74","CXCR4","SERPINE1","FOS","SELE","FABP4","CD36","CAVIN1","COL4A1","COL4A2")
res_endo_list<-list()
for (thresh in 80:240) {
  merged_filter_exclude_epi_endo$inf_endo<-FALSE
  merged_filter_exclude_epi_endo$inf_endo[which(apply(merged_filter_exclude_epi_endo@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')], 1, function(x) {sum(x<thresh)>0}))]<-T
  merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_endo)[which(merged_filter_exclude_epi_endo$inf_endo)])]<-'inf_endo'
  res_endo_list[[as.character(thresh)]]<-DE_wilcox(as.matrix(merged_filter_exclude_epi_endo@assays$RNA@data[endo_marker,]),merged_filter_exclude_epi_endo$inf_endo,TRUE,FALSE)
}
res_endo_list2<-lapply(endo_marker, function(gene){
  as.data.frame(do.call(rbind,lapply(res_endo_list, function(x) {x[gene,]})))->df
  rownames(df)->df$thresh
  df
})
names(res_endo_list2)<-endo_marker

lapply(res_endo_list2, function(df) {
  maxY<-max(abs(df$log2FC))
  as.numeric(df$thresh)->df$thresh
  ggplot(df, aes(x = thresh/4, y = log2FC, color = -log10(pvalue))) +
    geom_point(size = 1.5) +
    
    # draw the x-axis at y=0
    geom_hline(yintercept = 0, color='black',linetype='dashed', linewidth = 0.6) +
    scale_color_viridis_c() +
    # symmetric y-axis
    scale_y_continuous(limits = c(-maxY, maxY)) +
    geom_point(data = df['132', ], aes(x = thresh/4, y = log2FC),
               color = "orange", size = 1.4) +
    theme_classic()->p
  ggsave(paste0('~/test_diff_thresh/endo/p_',unique(df$gene),'.png'),p,width = 5,height = 4,dpi = 300)
})




