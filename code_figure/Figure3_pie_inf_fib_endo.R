load("/ix1/wchen/liutianhao/work_rdata/work_10_05.RData")
#dissect the cell type composition in the peri-fibroblasts regions
test_dist<-function(cell_type){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  lapply(dist_matrix_all, function(x) {
    as.data.frame(apply(x[,colnames(x)%in%names],1,min))
  })->all
  do.call(rbind,all)->all
  return(all[match(colnames(merged_filter_exclude_epi),gsub('.*[.]','',rownames(all))),1])
}

merged_filter_exclude_epi2 <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
#merged_filter_exclude_epi<-merged_filter_exclude_epi[,grep('1',merged_filter_exclude_epi$orig.ident)]
thresh_prop<-data.frame()
marker_inf_fib_list<-list()
marker_inf_fib<-FindMarkers(merged_filter_exclude_epi,group.by = 'cell_type_l1_all',ident.1 = 'inf_fib',ident.2 = 'Fibroblasts')
features <- c("LYZ","C3","CD74","B2M","TMSB4X","COL1A1","COL1A2","COL3A1","CXCL14","SPARC")
for (thresh in 80:240) {
  thresh<-132
  merged_filter_exclude_epi2->merged_filter_exclude_epi
  merged_filter_exclude_epi_fib<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts')]
  merged_filter_exclude_epi_fib$inf_fib<-FALSE
  merged_filter_exclude_epi_fib$inf_fib[which(apply(merged_filter_exclude_epi_fib@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')], 1, function(x) {sum(x<thresh)>0}))]<-T
  merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_fib)[which(merged_filter_exclude_epi_fib$inf_fib)])]<-'inf_fib'
  test_dist('inf_fib')->merged_filter_exclude_epi$min_dis_to_inf_fib
  test_dist('Fibroblasts')->merged_filter_exclude_epi$min_dis_to_nor_fib
  merged_filter_exclude_epi[,which(!merged_filter_exclude_epi$cell_type_l1_all%in%c('Fibroblasts','inf_fib'))]->exclude_target
  exclude_target$cell_type_l1_all[which(exclude_target$min_dis_to_inf_fib<thresh)]->inf_composite
  exclude_target$cell_type_l1_all[which(exclude_target$min_dis_to_nor_fib<thresh)]->nor_composite
  
  nor_composite<-nor_composite[which(nor_composite!='Neutrophil')]
  inf_composite<-factor(inf_composite,levels = c('Macrophage','T_cell','B_Plasma','Neutrophil','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
  nor_composite<-factor(nor_composite,levels = c('Macrophage','T_cell','B_Plasma','Neutrophil','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
  as.data.frame(prop.table(table(inf_composite)))->df
  #sum(df$Freq[which(df$inf_composite%in%c('T_cell','B_Plasma','Macrophage','Neutrophil'))])->thresh_prop[thresh,1]
  marker_inf_fib<-FindMarkers(merged_filter_exclude_epi,group.by = 'cell_type_l1_all',ident.1 = 'inf_fib',ident.2 = 'Fibroblasts',features = feature,logfc.threshold = 0,min.pct = 0)
  marker_inf_fib<-marker_inf_fib[feature,]
  marker_inf_fib->marker_inf_fib_list[[as.character(thresh)]]
}
saveRDS(marker_inf_fib_list,'/ix1/wchen/liutianhao/marker_inf_fib_list.rds')

write.csv(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all=='inf_fib'|(merged_filter_exclude_epi$cell_type_l1_all%in%c('Neutrophil','T_cell','B_Plasma','Macrophage')&merged_filter_exclude_epi$min_dis_to_inf_fib<132)),c('x','y','cell_type_l1_all')],'~/immune_fib/result_all.csv')

colors <- c(
  "Fibroblasts"="#00BBFF",  # deep teal
  "Endothelium"="#D95F02",  # orange
  "Glandular_epithelium"="#7570B3",  # indigo
  "Keratinocyte"="#8DD3C7",  # cyan (replaces gray)
  "Neuron"="#E7298A",  # magenta
  "Smooth_muscle_cell"="#A6761D",  # brown
  "Melanocyte"="#B2DF8A",  # green
  "B_Plasma"="#A6CEE3",  # sky blue
  "Macrophage"="#E6AB02",  # gold
  "T_cell"="#66A61E",
  "Neutrophil" = "#CAB2D6"
)



ggplot(as.data.frame(inf_composite), aes(x = "", fill = inf_composite)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()

ggplot(as.data.frame(nor_composite), aes(x = "", fill = nor_composite)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()


#find differential genes
marker_inf_fib<-FindMarkers(merged_filter_exclude_epi,group.by = 'cell_type_l1_all',ident.1 = 'inf_fib',ident.2 = 'Fibroblasts')
marker_inf_fib<-marker_inf_fib[which(marker_inf_fib$p_val_adj<0.01),]
marker_inf_fib<-marker_inf_fib[order(0-abs(marker_inf_fib$avg_log2FC)),]
merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all%in%c('Fibroblasts','inf_fib'))]->object
dot_plot(object,c(rownames(marker_inf_fib)[which(marker_inf_fib$avg_log2FC>0)][1:5],rownames(marker_inf_fib)[which(marker_inf_fib$avg_log2FC<0)][c(1:3,5,7)]),'cell_type_l1_all')+coord_flip()+theme(axis.text= element_text(size = 15))->p
ggsave('p.png',p,width = 7.5,height = 3.5,dpi = 300)


#dissect the cell type composition in the peri-vascular regions
merged_filter_exclude_epi2->merged_filter_exclude_epi
thresh<-132
merged_filter_exclude_epi_endo<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')]
merged_filter_exclude_epi_endo$inf_endo<-FALSE
merged_filter_exclude_epi_endo$inf_endo[which(apply(merged_filter_exclude_epi_endo@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')], 1, function(x) {sum(x<thresh)>0}))]<-T
merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_endo)[which(merged_filter_exclude_epi_endo$inf_endo)])]<-'inf_endo'
test_dist('inf_endo')->merged_filter_exclude_epi$min_dis_to_inf_endo
test_dist('Endothelium')->merged_filter_exclude_epi$min_dis_to_nor_endo
merged_filter_exclude_epi[,which(!merged_filter_exclude_epi$cell_type_l1_all%in%c('Endothelium','inf_endo'))]->exclude_target
exclude_target$cell_type_l1_all[which(exclude_target$min_dis_to_inf_endo<thresh)]->inf_composite
exclude_target$cell_type_l1_all[which(exclude_target$min_dis_to_nor_endo<thresh)]->nor_composite

write.csv(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all=='inf_endo'|(merged_filter_exclude_epi$cell_type_l1_all%in%c('Neutrophil','T_cell','B_Plasma','Macrophage')&merged_filter_exclude_epi$min_dis_to_inf_endo<132)),c('x','y','cell_type_l1_all')],'~/immune_endo/result_all.csv')

nor_composite<-nor_composite[which(nor_composite!='Neutrophil')]
inf_composite<-factor(inf_composite,levels = c('Macrophage','T_cell','B_Plasma','Neutrophil','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
nor_composite<-factor(nor_composite,levels = c('Macrophage','T_cell','B_Plasma','Neutrophil','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))

colors <- c(
  "Fibroblasts"="#00BBFF",  # deep teal
  "Endothelium"="#D95F02",  # orange
  "Glandular_epithelium"="#7570B3",  # indigo
  "Keratinocyte"="#8DD3C7",  # cyan (replaces gray)
  "Neuron"="#E7298A",  # magenta
  "Smooth_muscle_cell"="#A6761D",  # brown
  "Melanocyte"="#B2DF8A",  # green
  "B_Plasma"="#A6CEE3",  # sky blue
  "Macrophage"="#E6AB02",  # gold
  "T_cell"="#66A61E",
  "Neutrophil" = "#CAB2D6"
)

ggplot(as.data.frame(inf_composite), aes(x = "", fill = inf_composite)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()

ggplot(as.data.frame(nor_composite), aes(x = "", fill = nor_composite)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()

#find differential genes
marker_inf_endo<-FindMarkers(merged_filter_exclude_epi,group.by = 'cell_type_l1_all',ident.1 = 'inf_endo',ident.2 = 'Endothelium',recorrect_umi=F)
marker_inf_endo<-marker_inf_endo[which(marker_inf_endo$p_val<0.05),]
marker_inf_endo<-marker_inf_endo[order(0-abs(marker_inf_endo$avg_log2FC)),]
marker_inf_endo<-marker_inf_endo[order(marker_inf_endo$p_val),]
merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all%in%c('Endothelium','inf_endo'))]->object
dot_plot(object,rownames(marker_inf_endo)[which(marker_inf_endo$avg_log2FC<0)],'cell_type_l1_all')
dot_plot(object,rownames(marker_inf_endo)[which(marker_inf_endo$avg_log2FC>0)],'cell_type_l1_all')
dot_plot(object,c(rownames(marker_inf_endo)[which(marker_inf_endo$avg_log2FC>0)][1:5],rownames(marker_inf_endo)[which(marker_inf_endo$avg_log2FC<0)][c(1:5)]),'cell_type_l1_all')+coord_flip()+theme(axis.text= element_text(size = 15))->p

merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
#plot out the inflammed endothelium and the nearby immune cells
as.data.frame(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all%in%c('B_Plasma','inf_endo','Macrophage','Neutrophil','T_cell')),c('x','y','cell_type_l1_all')])->out
colnames(out)[3]<-'cell_type'
write.csv(out,'~/immune_endo/result_all.csv')
lapply(unique(merged_filter_exclude_epi$orig.ident), function(x) {
  write.csv(out[grep(x,rownames(out)),],paste0('~/immune_endo/result_',x,'.csv') )
})



#Fibrotic endothelium (Only these two cell types)
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
as.data.frame(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi$cell_type_l1_all%in%c('Fibroblasts','Endothelium')),c('x','y','cell_type_l1_all')])->out
colnames(out)[3]<-'cell_type'
write.csv(out,'~/fib_endo/result_all.csv')
lapply(unique(merged_filter_exclude_epi$orig.ident), function(x) {
  write.csv(out[grep(x,rownames(out)),],paste0('~/fib_endo/result_',x,'.csv') )
})

#correlation between marker gene expression and distance to immune cells
merged_filter_exclude_epi_endo_fib <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_endo_fib.rds")
merged_filter_exclude_epi_endo_fib<-merged_filter_exclude_epi_endo_fib[,which(merged_filter_exclude_epi_endo_fib$orig.ident!='G_2')]

DE_sc_fib_inter_LS <- readRDS("~/DE_sc_fib_inter_LS.rds")
DE_sc_fib_inter_SSC <- readRDS("~/DE_sc_fib_inter_SSC.rds")
DE_sc_endo_inter_LS <- readRDS("~/DE_sc_endo_inter_LS.rds")
DE_sc_endo_inter_SSC <- readRDS("~/DE_sc_endo_inter_SSC.rds")

endo<-merged_filter_exclude_epi_endo_fib[,which(merged_filter_exclude_epi_endo_fib$cell_type_l1_all%in%c('Endothelium','inf_endo'))]
AddModuleScore(endo[,which(endo$condition=='LS')],features = list(rownames(DE_sc_endo_inter_LS)),name = 'LS_endo_marker')->endo_LS
cor.test(endo_LS$min_dis_to_Macrophage,endo_LS$LS_endo_marker1)->test_LS_macro
cor.test(endo_LS$min_dis_to_T,endo_LS$LS_endo_marker1)->test_LS_T
AddModuleScore(endo[,which(endo$condition=='SSC')],features = list(rownames(DE_sc_endo_inter_SSC)),name = 'SSC_endo_marker')->endo_SSC
cor.test(endo_SSC$min_dis_to_Macrophage,endo_SSC$SSC_endo_marker1)->test_SSC_macro
cor.test(endo_SSC$min_dis_to_T,endo_SSC$SSC_endo_marker1)->test_SSC_T

matrix(c(test_LS_macro$estimate,test_SSC_macro$estimate,test_LS_T$estimate,test_SSC_T$estimate),nrow = 2,ncol = 2)->mtx
pheatmap::pheatmap(abs(mtx),cluster_rows = F,cluster_cols = F,border_color = NA)->p
ggsave('p.png',p,width = 4,height = 3,dpi = 300)

library(ggpubr)
library(ggplot2)
ggplot(endo_LS@meta.data[which(endo_LS$min_dis_to_Macrophage<2000),], aes(x =min_dis_to_Macrophage/4, y = LS_endo_marker1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.y = max(endo_LS$LS_endo_marker1[which(endo_LS$min_dis_to_Macrophage<2000)]*1.1)) +
  theme_classic(base_size = 14) +
  geom_vline(xintercept = 33,color='cyan',linetype = "dashed")+
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p
ggsave('p.png',p,width = 5,height = 5,dpi = 300)

ggplot(endo_SSC@meta.data[which(endo_SSC$min_dis_to_Macrophage<2000),], aes(x =min_dis_to_Macrophage/4, y = SSC_endo_marker1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.y = max(endo_SSC$SSC_endo_marker1[which(endo_SSC$min_dis_to_Macrophage<2000)]*1.1)) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")

fibro<-merged_filter_exclude_epi_endo_fib[,which(merged_filter_exclude_epi_endo_fib$cell_type_l1_all%in%c('Fibroblasts','inf_fib'))]
AddModuleScore(fibro[,which(fibro$condition=='LS')],features = list(rownames(DE_sc_fib_inter_LS)),name = 'LS_fib_marker')->fibro_LS
cor.test(fibro_LS$min_dis_to_Macrophage,fibro_LS$LS_fib_marker1)->test_LS_macro
cor.test(fibro_LS$min_dis_to_T,fibro_LS$LS_fib_marker1)->test_LS_T
AddModuleScore(fibro[,which(fibro$condition=='SSC')],features = list(rownames(DE_sc_fib_inter_SSC)),name = 'SSC_fib_marker')->fibro_SSC
cor.test(fibro_SSC$min_dis_to_Macrophage,fibro_SSC$SSC_fib_marker1)->test_SSC_macro
cor.test(fibro_SSC$min_dis_to_T,fibro_SSC$SSC_fib_marker1)->test_SSC_T

matrix(c(test_LS_macro$estimate,test_SSC_macro$estimate,test_LS_T$estimate,test_SSC_T$estimate),nrow = 2,ncol = 2)->mtx
pheatmap::pheatmap(abs(mtx),cluster_rows = F,cluster_cols = F,border_color = NA)->p
ggsave('p.png',p,width = 4,height = 3,dpi = 300)

ggplot(fibro_LS@meta.data[which(fibro_LS$min_dis_to_Macrophage<2000),], aes(x =min_dis_to_Macrophage/4, y = LS_fib_marker1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.y =  max(fibro_LS$LS_fib_marker1[which(fibro_LS$min_dis_to_Macrophage<2000)]*1.1)) +
  theme_classic(base_size = 14) +
  geom_vline(xintercept = 33,color='cyan',linetype = "dashed")+
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")->p
ggsave('p.png',p,width = 5,height = 5,dpi = 300)

ggplot(fibro_SSC@meta.data[which(fibro_SSC$min_dis_to_Macrophage<2000),], aes(x =min_dis_to_Macrophage/4, y = SSC_fib_marker1)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  stat_cor(method = "pearson", label.y = max(fibro_SSC$SSC_fib_marker1[which(fibro_SSC$min_dis_to_Macrophage<2000)]*1.1)) +
  theme_classic(base_size = 14) +
  labs(x = "X variable", y = "Y variable", title = "Negative correlation between X and Y")



#dissect the cell type composition for each fibroblast or endothelial cells
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident!='G_1')]
test_dist<-function(cell_type){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  as.data.frame(merged_filter_exclude_epi@meta.data[!colnames(merged_filter_exclude_epi)%in%names,'cell_type_l1_all'])->meta
  rownames(meta)<-colnames(merged_filter_exclude_epi)[!colnames(merged_filter_exclude_epi)%in%names]
  colnames(meta)<-'cell_type'
  all_labels<-unique(merged_filter_exclude_epi$cell_type_l1_all)
  lapply(dist_matrix_all, function(x) {
    apply(x[!rownames(x)%in%names,colnames(x)%in%names],2,function(y) {
      meta[rownames(x)[which(y<thresh)],'cell_type']->labels
      table(factor(labels, levels = all_labels))->counts
      as.numeric(counts) / length(labels)
    })
  })->all
  
  do.call(cbind,all)->all
  as.data.frame(t(all))->all
  colnames(all)<-all_labels
  return(all)
}
thresh<-220
test_dist('Fibroblasts')->neighbor_fibroblast
test_dist('Endothelium')->neighbor_endo
neighbor_fibroblast[is.na(neighbor_fibroblast)] <- 0
neighbor_endo[is.na(neighbor_endo)] <- 0
neighbor_fibroblast_filter<-neighbor_fibroblast[which(rowSums(neighbor_fibroblast)>0),]
neighbor_endo_filter<-neighbor_endo[which(rowSums(neighbor_endo)>0),]


dist_cols <- dist(neighbor_fibroblast_filter, method = "euclidean")
hc <- hclust(dist_cols, method = "ward.D2")
plot(hc, labels = FALSE)
clusters <- cutree(hc, k = 8)
neighbor_fibroblast_filter$cluster<-clusters
apply(neighbor_fibroblast_filter[,-dim(neighbor_fibroblast_filter)[2]], 1, function(x) {
  colnames(neighbor_fibroblast_filter)[which.max(x)]
})->neighbor_fibroblast_filter$max
table(neighbor_fibroblast_filter$max,neighbor_fibroblast_filter$cluster)

neighbor_fibroblast_filter$immune_tag<-apply(neighbor_fibroblast_filter[,c('B_Plasma','T_cell','Macrophage','Neutrophil')], 1,function(x) {sum(x>0)>0})
table(neighbor_fibroblast_filter$immune_tag,neighbor_fibroblast_filter$cluster)
rownames(neighbor_fibroblast_filter)[which(neighbor_fibroblast_filter$cluster%in%c(2,3,8))]->inf_fib_index


dist_cols <- dist(neighbor_endo_filter, method = "euclidean")
hc <- hclust(dist_cols, method = "ward.D2")
plot(hc, labels = FALSE)
clusters <- cutree(hc, k = 8)
neighbor_endo_filter$cluster<-clusters
apply(neighbor_endo_filter[,-dim(neighbor_endo_filter)[2]], 1, function(x) {
  colnames(neighbor_endo_filter)[which.max(x)]
})->neighbor_endo_filter$max
table(neighbor_endo_filter$max,neighbor_endo_filter$cluster)

neighbor_endo_filter$immune_tag<-apply(neighbor_endo_filter[,c('B_Plasma','T_cell','Macrophage','Neutrophil')], 1,function(x) {sum(x>0)>0})
table(neighbor_endo_filter$immune_tag,neighbor_endo_filter$cluster)
rownames(neighbor_endo_filter)[which(neighbor_endo_filter$cluster%in%c(1,2,4,7))]->inf_endo_index


