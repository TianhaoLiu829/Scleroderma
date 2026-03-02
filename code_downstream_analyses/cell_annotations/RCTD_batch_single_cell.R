#RCTD deconvolution
list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/')->sample
reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Smooth_muscle_cell", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
#merged <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/merged.rds")
merged <- readRDS('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_unfilter.rds')
merged$condition[which(merged$condition=='HC')]<-'LS'
RCTD_DECON<-function(object,target,ref,thresh){
  object[,which(object$orig.ident==target)]->object
  PercentageFeatureSet(object,pattern = 'MT-')->object$mt_percentage
  object<-object[,which(object$nCount_RNA>50&object$mt_percentage<20)]
  ref<-ref[,which(ref$condition==unique(object$condition))]
  cluster <- as.factor(ref$cell_type_l1)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  refe <- Reference(ref@assays[["RNA"]]@counts, cluster, nUMI)
  counts <-object@assays$RNA@counts
  rownames(counts)<-rownames(object)
  colnames(counts)<-colnames(object)
  
  coords <- read.csv(paste0("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/",target,"/","center_",target,".csv"))
  coords<-coords[,c("x","y")]
  rownames(coords)<-as.character(0:(dim(coords)[1]-1))
  rownames(coords)<-paste0(target,'_',rownames(coords))
  coords<-coords[colnames(object),]
  query <- SpatialRNA(coords, counts, colSums(counts))
  RCTD <- create.RCTD(query, refe, max_cores = 30,UMI_min = 50,UMI_min_sigma = 200,CELL_MIN_INSTANCE=5)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  object <- AddMetaData(object, metadata = RCTD@results$results_df)
  object<-object[,rownames(RCTD@results[["weights"]])]
  object[["conv"]] <- CreateAssayObject(counts =Transpose(RCTD@results[["weights"]]))
  coords<-coords[colnames(object),]
  cbind(colnames(object),coords[,c('x','y')],t(as.matrix(object[['conv']]$counts)),object$first_type)->result
  colnames(result)[c(1,dim(result)[2])]<-c('index','cell_type')
  write.csv(result,paste0('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/all_reference/doublet/unfilter_thresh50/result_',target,'.csv'),row.names = FALSE)
  return(object)
}

.libPaths(c("/software/rhel9/manual/install/r/4.5.0/lib64/R/library","/ihome/wchen/tianhao/R/x86_64-pc-linux-gnu-library/4.5"))
library(Seurat)
library(SeuratData)
library(SeuratDisk)
.libPaths(c("/ihome/wchen/tianhao/R/x86_64-pc-linux-gnu-library/4.5","/software/rhel9/manual/install/r/4.5.0/lib64/R/library"))
library(spacexr)

#lapply(sample[which(nchar(sample)==3)][c(9,1,2,3,4,5,6,7,8,10)], function(x) {RCTD_DECON(merged,x,reference,5)})->result
lapply(unique(merged$orig.ident), function(x) {RCTD_DECON(merged,x,reference,5)})->result
#names(result)<-sample[1:10]
names(result)<-unique(merged$orig.ident)
saveRDS(result,'/ix1/wchen/liutianhao/result_l1_2.rds')

wqw
#calculate proportion
list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/')->sample
result<-list()
for(i in sample[which(nchar(sample)==3)][-c(2,14)]){
  read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/all_reference/doublet/unfilter_thresh50/result_',i,'.csv'))->result[[i]]
} 
do.call(rbind,result)->result
merged_filter$cell_type_l1_all<-result$cell_type[match(colnames(merged_filter),result$index)]
#find marker genes
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter.rds")
merged_filter_base<-merged_filter[,-grep('2',merged_filter$orig.ident)]
merged_filter_base@meta.data<-as.data.frame(cbind(merged_filter_base@meta.data,result[match(colnames(merged_filter_base),result$index),-c(1,2,3,14)]))
merged_filter_base<-sc_pro(merged_filter_base,FALSE)
#marker and heatmap
merged_filter_base_exclude_epi<-merged_filter_base[,which(!merged_filter_base$epidermis)]
sc_pro(merged_filter_base_exclude_epi,FALSE)->merged_filter_base_exclude_epi
FindAllMarkers(merged_filter_base_exclude_epi,only.pos = T)->marker_base_exclude_epi
library(dplyr)
marker_base_exclude_epi %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(merged_filter_base_exclude_epi,group.by = 'cell_type_l1_all',features = top10$gene)
merged_filter_base_exclude_epi@assays$RNA@var.features<-union(merged_filter_base_exclude_epi@assays$RNA@var.features,top10$gene)
merged_filter_base_exclude_epi<-ScaleData(merged_filter_base_exclude_epi)
merged_filter_base_exclude_epi$cell_type_l1_all<-factor(merged_filter_base_exclude_epi$cell_type_l1_all,levels = levels(top10$cluster))
DoHeatmap(merged_filter_base_exclude_epi,group.by = 'cell_type_l1_all',features = top10$gene)->p
ggsave('p2.png',p,width = 14,height = 9,dpi = 300)

#find marker
SSC_marker_list<-list()
for (celltype in unique(merged_filter_base$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi[,which(merged_filter_base$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'wilcox',latent.vars = colnames(merged_filter_base@meta.data)[8:17])->SSC_marker_list[[celltype]]
}
#write down the marker genes
library(openxlsx)
wb <- createWorkbook()
for (name in names(SSC_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = SSC_marker_list[[name]])
}
saveWorkbook(wb, file = "multiple_sheets.xlsx", overwrite = TRUE)
#markers for cell types
Idents(merged_filter_base)<-merged_filter_base$cell_type_l1_all
FindAllMarkers(merged_filter_base,only.pos = T)->marker_cell_type
df_list <- split(marker_cell_type, marker_cell_type$cluster)
wb <- createWorkbook()
for (label in names(df_list)) {
  addWorksheet(wb, sheetName = label)
  writeData(wb, sheet = label, x = df_list[[label]],rowNames = T)
}
saveWorkbook(wb, "split_by_cell_type.xlsx", overwrite = TRUE)

#do the excluding epidermis analysis
merged_filter_base<-merged_filter[,-grep('2',merged_filter$orig.ident)]
merged_filter_base_exclude_epi<-merged_filter_base[,which(!merged_filter_base$epidermis)]
merged_filter_base_exclude_epi@meta.data<-as.data.frame(cbind(merged_filter_base_exclude_epi@meta.data,result[match(colnames(merged_filter_base_exclude_epi),result$index),-c(1,2,3,14)]))
merged_filter_base_exclude_epi@meta.data[which(colnames(merged_filter_base_exclude_epi)%in%colnames(macro)[which(macro$seurat_clusters==3)]),'cell_type_l1_all']<-'Neutrophil'
SSC_marker_list<-list()
for (celltype in unique(merged_filter_base_exclude_epi$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi[,which(merged_filter_base_exclude_epi$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'wilcox',latent.vars = colnames(merged_filter_base_exclude_epi@meta.data)[10:19])->SSC_marker_list[[celltype]]
}
#write down the marker genes
library(openxlsx)
wb <- createWorkbook()
for (name in names(SSC_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = SSC_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/multiple_sheets_exclude_epi2.xlsx", overwrite = TRUE)

#markers for cell types
Idents(merged_filter_base)<-merged_filter_base$cell_type_l1_all
FindAllMarkers(merged_filter_base,only.pos = T)->marker_cell_type
df_list <- split(marker_cell_type, marker_cell_type$cluster)
wb <- createWorkbook()
for (label in names(df_list)) {
  addWorksheet(wb, sheetName = label)
  writeData(wb, sheet = label, x = df_list[[label]])
}
saveWorkbook(wb, "split_by_cell_type_exclude_epi.xlsx", overwrite = TRUE)

#markers for pre post treatment
merged_filter_treat<-merged_filter[,which(merged_filter$orig.ident%in%c('E_1','E_2','K_1','K_2','H_1','H_2')) ]
sc_pro(merged_filter_treat,FALSE)->merged_filter_treat
merged_filter_treat_exclude_epi<-merged_filter_treat[,which(!merged_filter_treat$epidermis)]
merged_filter_treat_exclude_epi@meta.data<-as.data.frame(cbind(merged_filter_treat_exclude_epi@meta.data,result[match(colnames(merged_filter_treat_exclude_epi),result$index),-c(1,2,3,14)]))
sc_pro(merged_filter_treat_exclude_epi,FALSE)->merged_filter_treat_exclude_epi

merged_filter_treat_exclude_epi@meta.data[which(colnames(merged_filter_treat_exclude_epi)%in%colnames(macro)[which(macro$seurat_clusters==3)]),'cell_type_l1_all']<-'Neutrophil'

merged_filter_treat_exclude_epi$pre_post<-'pre'
merged_filter_treat_exclude_epi$pre_post[grep('2',merged_filter_treat_exclude_epi$orig.ident)]<-'post'
merged_filter_treat$pre_post<-'pre'
merged_filter_treat$pre_post[grep('2',merged_filter_treat$orig.ident)]<-'post'

treat_marker_list<-list()
for (celltype in unique(merged_filter_treat_exclude_epi$cell_type_l1_all)) {
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('E_1','E_2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_E
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('H_1','H_2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_H
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('K_1','K_2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_K
  inte<-intersect(intersect(rownames(treat_marker_K),rownames(treat_marker_H)),rownames(treat_marker_E))
  cbind(treat_marker_E[inte,c(1,2,5)],treat_marker_H[inte,c(1,2,5)],treat_marker_K[inte,c(1,2,5)])->treat_marker_list[[celltype]]
}
library(openxlsx)
wb <- createWorkbook()
for (name in names(treat_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = treat_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/treat_exclude_epi.xlsx", overwrite = TRUE)




#compare the cell type proportion between different samples (treatment and LS/SSc)

library(dplyr)
df_agg <- merged_filter_base@meta.data %>%
  group_by(orig.ident, condition, cell_type_l1_all) %>%
  summarise(count = n()) %>%
  ungroup()
df_agg <- df_agg %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count))
library(ggplot2)
df_agg$orig.ident<-factor(df_agg$orig.ident,levels = c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))
ggplot(df_agg, aes(x = orig.ident, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ cell_type_l1_all, scales = "free_y") +  # One plot for each cell type
  labs(x = "Sample", y = "Proportion", fill = "Disease Type") +
  theme_minimal()->p
ggsave('p.png',p,width = 12,height = 10,dpi = 300)

df_agg <- merged_filter_treat@meta.data %>%
  group_by(orig.ident, pre_post, cell_type_l1_all) %>%
  summarise(count = n()) %>%
  ungroup()
df_agg <- df_agg %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count))
library(ggplot2)
df_agg$orig.ident<-factor(df_agg$orig.ident,levels = c('E_1','E_2','H_1','H_2','K_1','K_2'))
df_agg$pre_post<-factor(df_agg$pre_post,levels = c('pre','post'))
ggplot(df_agg, aes(x = orig.ident, y = proportion, fill = pre_post)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ cell_type_l1_all, scales = "free_y") +  # One plot for each cell type
  labs(x = "Sample", y = "Proportion", fill = "Disease Type") +
  theme_minimal()->p
ggsave('p.png',p,width = 12,height = 10,dpi = 300)

# compare the cell type proportion with flex samples
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter.rds")
readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K/combine_all_base.rds")->combine_all_base
combine_flex_only <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_flex_only.rds")
library(dplyr)
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
combine_flex_only@meta.data$cell_type_l1<-array[combine_flex_only$cell_type_high]
combine_all_base@meta.data$cell_type_l1<-array[combine_all_base$cell_type_high]

#exclude krt
combine_flex_only<-combine_flex_only[,which(combine_flex_only$cell_type_l1!='Keratinocyte')]
combine_flex_only@meta.data %>% count(cell_type_l1,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_flex
combine_all_base<-combine_all_base[,which(combine_all_base$cell_type_l1!='Keratinocyte')]
combine_all_base@meta.data %>% count(cell_type_l1,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_regular
merged_filter@meta.data[which(gsub('_','',merged_filter$orig.ident)%in%unique(combine_flex_only$orig.ident)),]->target
target[which(target$cell_type_l1_all!='Keratinocyte'),] %>% count(cell_type_l1_all,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_ST
combine_flex<-cbind(prop_ST,prop_flex[match(paste0(prop_ST$cell_type_l1_all,gsub('_','',prop_ST$orig.ident) ) ,paste0(prop_flex$cell_type_l1,prop_flex$orig.ident)),])
merged_filter@meta.data[which(gsub('_','',merged_filter$orig.ident)%in%unique(combine_all_base$orig.ident)),]->target
target[which(target$cell_type_l1_all!='Keratinocyte'),] %>% count(cell_type_l1_all,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_ST
combine_regular<-cbind(prop_ST,prop_regular[match(paste0(prop_ST$cell_type_l1_all,gsub('_','',prop_ST$orig.ident) ) ,paste0(prop_regular$cell_type_l1,prop_regular$orig.ident)),])

for (i in c('E_1','F_1','H_1','K_1','J_1')) {
  p<- ggplot(
    combine_regular[which(combine_regular$orig.ident...2 == i & combine_regular$cell_type_l1_all != 'Keratinocyte'),
                    c(1,4,8)],
    aes(x = proportion...4, y = proportion...8,
        color = cell_type_l1_all, label = cell_type_l1_all)
  ) +
    geom_point(size = 3) +                                # scatter points
    geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
    theme_minimal() +
    theme(legend.position = "right") +                    # keep legend
    labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
    scale_color_brewer(palette = "Set1")                  # more contrastive colors
  ggsave(paste0('p_',i,'_regular.png'),p,width=8,height=6,dpi=300)
}

for (i in c('I_1','H_2','G_2','M_1')) {
  p<- ggplot(
    combine_flex[which(combine_flex$orig.ident...2 == i & combine_flex$cell_type_l1_all != 'Keratinocyte'),
                 c(1,4,8)],
    aes(x = proportion...4, y = proportion...8,
        color = cell_type_l1_all, label = cell_type_l1_all)
  ) +
    geom_point(size = 3) +                                # scatter points
    geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
    theme_minimal() +
    theme(legend.position = "right") +                    # keep legend
    labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
    scale_color_brewer(palette = "Set1")                  # more contrastive colors
  ggsave(paste0('p_',i,'_flex.png'),p,width=8,height=6,dpi=300)
}
combine_regular<-combine_regular[which(!is.na(combine_regular$proportion...8)),]
cor(combine_regular$proportion...4,combine_regular$proportion...8)
combine_flex<-combine_flex[which(!is.na(combine_flex$proportion...8)),]
cor(combine_flex$proportion...4,combine_flex$proportion...8)

#include krt
combine_flex_only@meta.data %>% count(cell_type_l1,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_flex
combine_all_base@meta.data %>% count(cell_type_l1,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_regular
merged_filter@meta.data[which(gsub('_','',merged_filter$orig.ident)%in%unique(combine_flex_only$orig.ident)),]->target
target%>% count(cell_type_l1_all,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_ST
combine_flex<-cbind(prop_ST,prop_flex[match(paste0(prop_ST$cell_type_l1_all,gsub('_','',prop_ST$orig.ident) ) ,paste0(prop_flex$cell_type_l1,prop_flex$orig.ident)),])
merged_filter@meta.data[which(gsub('_','',merged_filter$orig.ident)%in%unique(combine_all_base$orig.ident)),]->target
target%>% count(cell_type_l1_all,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_ST
combine_regular<-cbind(prop_ST,prop_regular[match(paste0(prop_ST$cell_type_l1_all,gsub('_','',prop_ST$orig.ident) ) ,paste0(prop_regular$cell_type_l1,prop_regular$orig.ident)),])

for (i in c('E_1','F_1','H_1','K_1','J_1')) {
  p<- ggplot(
    combine_regular[which(combine_regular$orig.ident...2 == i & combine_regular$cell_type_l1_all != 'Keratinocyte'),
                    c(1,4,8)],
    aes(x = proportion...4, y = proportion...8,
        color = cell_type_l1_all, label = cell_type_l1_all)
  ) +
    geom_point(size = 3) +                                # scatter points
    geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
    theme_minimal() +
    theme(legend.position = "right") +                    # keep legend
    labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
    scale_color_brewer(palette = "Set1")                  # more contrastive colors
  ggsave(paste0('p_',i,'_regular_all.png'),p,width=8,height=6,dpi=300)
}

for (i in c('I_1','H_2','G_2','M_1')) {
  p<- ggplot(
    combine_flex[which(combine_flex$orig.ident...2 == i & combine_flex$cell_type_l1_all != 'Keratinocyte'),
                 c(1,4,8)],
    aes(x = proportion...4, y = proportion...8,
        color = cell_type_l1_all, label = cell_type_l1_all)
  ) +
    geom_point(size = 3) +                                # scatter points
    geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
    theme_minimal() +
    theme(legend.position = "right") +                    # keep legend
    labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
    scale_color_brewer(palette = "Set1")                  # more contrastive colors
  ggsave(paste0('p_',i,'_flex_all.png'),p,width=8,height=6,dpi=300)
}
combine_regular2<-combine_regular[which(!is.na(combine_regular$proportion...8)&combine_regular$cell_type_l1_all != 'Keratinocyte'),]
combine_regular2<-combine_regular[which(!is.na(combine_regular$proportion...8)),]
cor(combine_regular2$proportion...4,combine_regular2$proportion...8)
combine_flex2<-combine_flex[which(!is.na(combine_flex$proportion...8)&combine_flex$cell_type_l1_all != 'Keratinocyte'),]
combine_flex2<-combine_flex[which(!is.na(combine_flex$proportion...8)),]
cor(combine_flex2$proportion...4,combine_flex2$proportion...8)

#permutation all
readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds")->combine_all
combine_all<-combine_all[,which(grepl('1',combine_all$orig.ident)|combine_all$orig.ident%in%c('H2','G2'))]
combine_all<-combine_all[,which(combine_all$orig.ident!='G1')]
library(dplyr)
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
combine_all@meta.data$cell_type_l1<-array[combine_all$cell_type_high]

combine_all@meta.data %>% count(cell_type_l1,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_all
merged_filter@meta.data->target
target%>% count(cell_type_l1_all,orig.ident) %>% group_by(orig.ident) %>% mutate(proportion = n / sum(n))->prop_ST
prop_ST<-prop_ST[which(!is.na(prop_ST$cell_type_l1_all)),]
gsub('_','',prop_ST$orig.ident)->prop_ST$orig.ident
intersect(unique(prop_all$orig.ident),unique(prop_ST$orig.ident))->inter
cor_mtx<-matrix(ncol=length(inter),nrow=length(inter))
rownames(cor_mtx)<-inter
colnames(cor_mtx)<-inter
for (i in inter){
  for (j in inter){
    prop_ST[which(prop_ST$orig.ident==i),]->ST
    prop_all[which(prop_all$orig.ident==j),]->regular
    cbind(regular,ST[match(regular$cell_type_l1,ST$cell_type_l1_all),])->combine
    combine<-combine[which(combine$cell_type_l1_all!='Keratinocyte'),]
    cor(combine[,4],combine[,8],method = 'spearman')->cor_mtx[i,j]
  }
}

inter2<-unique(prop_ST$cell_type_l1_all)
cor_mtx<-matrix(ncol=length(inter2),nrow=length(inter2))
rownames(cor_mtx)<-inter2
colnames(cor_mtx)<-inter2
for (i in inter2){
  for (j in inter2){
    prop_ST[which(prop_ST$cell_type_l1_all==i&prop_ST$orig.ident%in%inter),]->ST
    prop_all[which(prop_all$cell_type_l1==j&prop_all$orig.ident%in%inter),]->regular
    cbind(regular,ST[match(regular$orig.ident,ST$orig.ident),])->combine
    cor(combine[,4],combine[,8])->cor_mtx[i,j]
  }
}

pheatmap::pheatmap(cor_mtx,cluster_cols=FALSE,cluster_rows=FALSE)->p
ggsave('p_cell.png',p,width=9,height=9)

i<-'M1'
j<-'M1'
prop_ST[which(prop_ST$orig.ident==i),]->ST
prop_all[which(prop_all$orig.ident==j),]->regular
cbind(regular,ST[match(regular$cell_type_l1,ST$cell_type_l1_all),])->combine
combine<-combine[which(!combine$cell_type_l1_all%in%c('Keratinocyte','Glandular_epithelium','Fibroblasts')),]
#combine<-combine[which(!combine$cell_type_l1_all%in%c('Keratinocyte')),]
cor(combine[,4],combine[,8],method = 'spearman')
combine[,c(4,8)]<-apply(combine[,c(4,8)],2,function(x) log(x/(1-x)) )
lims <- range(c(combine[,4], combine[,8]))
ggplot(
  combine[,c(1,4,8)],
  aes(x = proportion...4, y = proportion...8,color = cell_type_l1)
) +
  scale_color_manual(values = my_colors) +
  geom_text_repel(aes(label = cell_type_l1), 
                  show.legend = FALSE,   # avoid duplicating "a" in legend
                  max.overlaps = Inf,size=4.5) + 
  geom_point(size = 3) +                                # scatter points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
  theme_classic() +
  xlim(lims)+
  ylim(lims)+
  theme(legend.position = "right",axis.text.x=element_text(size=13),axis.text.y=element_text(size=13)) +                    # keep legend
  labs(x = "proportion in Flex", y = "proportion in ST", color = "Cell type") +
  scale_color_brewer(palette = "Set2")  ->p                # more contrastive colors
my_colors <- c(
  "#1f77b4",  # blue
  "#ff7f0e",  # orange
  "#2ca02c",  # green
  "#d62728",  # red
  "#9467bd",  # purple
  "#8c564b",  # brown
  "#e377c2",  # pink
  "#7f7f7f",  # gray
  "#bcbd22",  # olive
  "#17becf"   # cyan
)


combine_all_m1<-combine_all[,which(combine_all$orig.ident=='M1')]
merged_filter_m1<-merged_filter[,which(merged_filter$orig.ident=='M_1')]
Idents(combine_all_m1)<-combine_all_m1$cell_type_l1
FindMarkers(combine_all_m1,only.pos = T,ident.1 = 'Lymphocyte')->marker_sc_t
FindMarkers(combine_all_m1,only.pos = T,ident.1 = 'Fibroblasts')->marker_sc_fib
FindMarkers(combine_all_m1,only.pos = T,ident.1 = 'Myeloid')->marker_sc_mye
Idents(merged_filter_m1)<-merged_filter_m1$cell_type_l1_all
FindMarkers(merged_filter_m1,only.pos = T,ident.1 = 'Lymphocyte')->marker_st_t
FindMarkers(merged_filter_m1,only.pos = T,ident.1 = 'Fibroblasts')->marker_st_fib
FindMarkers(merged_filter_m1,only.pos = T,ident.1 = 'Myeloid')->marker_st_mye
inter_marker_t<-intersect(rownames(marker_sc_t)[which(marker_sc_t$p_val_adj<0.05)],rownames(marker_st_t)[which(marker_st_t$p_val_adj<0.05)])
inter_marker_fib<-intersect(rownames(marker_sc_fib)[which(marker_sc_fib$p_val_adj<0.05)],rownames(marker_st_fib)[which(marker_st_fib$p_val_adj<0.05)])
inter_marker_mye<-intersect(rownames(marker_sc_mye)[which(marker_sc_mye$p_val_adj<0.05)],rownames(marker_st_mye)[which(marker_st_mye$p_val_adj<0.05)])

a<-data.frame(ST=rowMeans(merged_filter_m1@assays$RNA@counts[inter_marker_t,]),SC=rowMeans(combine_all_m1@assays$RNA@counts[inter_marker_t,]))
a_t<-data.frame(ST=rowMeans(merged_filter_m1@assays$RNA@counts[inter_marker_t,which(merged_filter_m1$cell_type_l1_all=='Lymphocyte')]),SC=rowMeans(combine_all_m1@assays$RNA@counts[inter_marker_t,which(combine_all_m1$cell_type_l1=='Lymphocyte')]))
a_fib<-data.frame(ST=rowMeans(merged_filter_m1@assays$RNA@counts[inter_marker_fib,which(merged_filter_m1$cell_type_l1_all=='Fibroblasts')]),SC=rowMeans(combine_all_m1@assays$RNA@counts[inter_marker_fib,which(combine_all_m1$cell_type_l1=='Fibroblasts')]))
a_mye<-data.frame(ST=rowMeans(merged_filter_m1@assays$RNA@counts[inter_marker_mye,which(merged_filter_m1$cell_type_l1_all=='Myeloid')]),SC=rowMeans(combine_all_m1@assays$RNA@counts[inter_marker_mye,which(combine_all_m1$cell_type_l1=='Myeloid')]))

ggplot(a_fib,aes(x = log(SC), y = log(ST))) +
  geom_point(size = 1.5) +                                # scatter points
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, color = "blue",size = 0.5)+
  theme(legend.position = "right",axis.text.x=element_text(size=13),axis.text.y=element_text(size=13)) +                    # keep legend
  labs(x = "proportion in Flex", y = "proportion in ST", color = "Cell type") +
  scale_color_brewer(palette = "Set2")  ->p                # more contrastive colors
ggsave('p_fib.png',p,width = 6,height = 5,dpi = 300)

library(ggrepel)
p_E_1<- ggplot(
  combine[which(combine$orig.ident...2 == 'E_1' & combine$cell_type_l1_all != 'Keratinocyte'),
          c(1,4,8)],
  aes(x = proportion...4, y = proportion...8,
      color = cell_type_l1_all, label = cell_type_l1_all)
) +
  geom_point(size = 3) +                                # scatter points
  geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
  theme_minimal() +
  theme(legend.position = "right") +                    # keep legend
  labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
  scale_color_brewer(palette = "Set1")                  # more contrastive colors
ggsave('p_E_1.png',p_E_1,width=8,height=6,dpi=300)

p_H_2 <- ggplot(
  combine[which(combine$orig.ident...2 == 'H_2' & combine$cell_type_l1_all != 'Keratinocyte'),
          c(1,4,8)],
  aes(x = proportion...4, y = proportion...8,
      color = cell_type_l1_all, label = cell_type_l1_all)
) +
  geom_point(size = 3) +                                # scatter points
  geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
  theme_minimal() +
  theme(legend.position = "right") +                    # keep legend
  labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
  scale_color_brewer(palette = "Set1")                  # more contrastive colors

p_I_1 <- ggplot(
  combine[which(combine$orig.ident...2 == 'I_1' & combine$cell_type_l1_all != 'Keratinocyte'),
          c(1,4,8)],
  aes(x = proportion...4, y = proportion...8,
      color = cell_type_l1_all, label = cell_type_l1_all)
) +
  geom_point(size = 3) +                                # scatter points
  geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
  theme_minimal() +
  theme(legend.position = "right") +                    # keep legend
  labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
  scale_color_brewer(palette = "Set1")                  # more contrastive colors

p_M_1 <- ggplot(
  combine[which(combine$orig.ident...2 == 'M_1' & combine$cell_type_l1_all != 'Keratinocyte'),
          c(1,4,8)],
  aes(x = proportion...4, y = proportion...8,
      color = cell_type_l1_all, label = cell_type_l1_all)
) +
  geom_point(size = 3) +                                # scatter points
  geom_text(vjust = -0.5, hjust = 0.5) +                # text labels
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # y = x line
  theme_minimal() +
  theme(legend.position = "right") +                    # keep legend
  labs(x = "proportion in ST", y = "proportion in flex", color = "Cell type") +
  scale_color_brewer(palette = "Set1")                  # more contrastive colors


lapply(result, function(x) {prop.table(table(x$cell_type))})->prop_result
prop_result<-lapply(prop_result, function(x) as.data.frame(x))
# Add an identifier to each
for (i in seq_along(prop_result)) {
  prop_result[[i]]$sample <- names(result)[i]
}

combined <- do.call(rbind, prop_result)
colnames(combined) <- c("cell_type", "proportion", "sample")

#get the cell count from Visium data and E_1+F_1
get_prop<-function(object,sample){
  sweep(object@assays$conv@counts, 2, colSums(object@assays$conv@counts), FUN = "/")->mat
  mat<-as.matrix(mat)
  apply(mat,1,function(x) sum(x*object$nCount_RNA)/sum(object$nCount_RNA))->prop_result
  return(data.frame(cell_type=names(prop_result),proportion=as.numeric(prop_result),sample=sample))
}
combined <-combined [which(!combined$sample%in%c('E_1','F_1')),]

A_2_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/A_2_RCTD.rds")
A_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/A_1_RCTD.rds")
H_2_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/H_2_RCTD.rds")
H_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/H_1_RCTD.rds")
E_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/E_1_RCTD.rds")
E_2_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/SSC/E_2_RCTD.rds")
F_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/F_1_RCTD.rds")
G_2_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/G_2_RCTD.rds")
J_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/J_1_RCTD.rds")
G_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/G_1_RCTD.rds")
I_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/I_1_RCTD.rds")
B_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/B_1_RCTD.rds")
C_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/C_1_RCTD.rds")
D_1_RCTD <- readRDS("/ix1/wchen/liutianhao/result/skin_ST/HD/high_thresh/final_deconvolution/high/D_1_RCTD.rds")

for (i in c('A_1','B_1','C_1','D_1','E_1','F_1')) {
  get_prop(get(paste0(i,'_RCTD')),i)->p
  p$cell_type<-gsub('-','_',p$cell_type)
  rbind(combined,p)->combined
}

for (i in c('A_1','B_1','C_1','D_1')) {
  get_prop(get(paste0(i,'_RCTD')),i)->p
  p$cell_type<-gsub('-','_',p$cell_type)
  rbind(combined,p)->combined
}

library(ggplot2)
combined_LS_SSC<-combined[which(combined$sample%in%c('G_2','A_1','E_1','H_1','B_1','C_1','D_1','F_1','G_1','I_1','J_1')),]
combined_treat<-combined[which(combined$sample%in%c('A_1','A_2','E_1','E_2','H_1','H_2','K_1','K_2')),]
combined_treat$sample<-factor(combined_treat$sample,levels = c('A_1','A_2','E_1','E_2','H_1','H_2'))


as.numeric(table(E_1_RCTD$first_type)/sum(table(E_1_RCTD$first_type)))->E_1
combined$proportion[which(combined$sample=='E_1')]<-E_1
as.numeric(table(F_1_RCTD$first_type)/sum(table(F_1_RCTD$first_type)))->F_1
combined$proportion[which(combined$sample=='F_1')]<-F_1



#plot
combined_LS_SSC$condition<-'LS'
combined_LS_SSC$condition[which(combined_LS_SSC$sample%in%c('A_1','E_1','H_1'))]<-'SSc'
combined_LS_SSC$condition[which(combined_LS_SSC$sample%in%c('G_2'))]<-'HC'
combined_LS_SSC$sample[which(combined_LS_SSC$sample=='G_2')]<-'G_healthy'
combined_LS_SSC$sample[which(combined_LS_SSC$condition!='HC')]<-gsub('_.*','',combined_LS_SSC$sample[which(combined_LS_SSC$condition!='HC')])
combined_LS_SSC$sample<-factor(combined_LS_SSC$sample,levels = c('G_healthy','A','E','H','B','C','D','F','G','I','J'))
ggplot(combined_LS_SSC, aes(x = sample, y = proportion,fill=condition)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cell_type, scales = c("free")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # show and rotate x labels
    axis.ticks.x = element_line()                     # optional, improves layout
  )+
  ylab("Proportion") +
  xlab("Sample") +
  ggtitle("Proportion of Each Label Across Samples")->p
ggsave('p.png',p,width = 12,height = 9,dpi = 300)

library(ggplot2)
combined_treat$patient<-gsub('_.*','',combined_treat$sample)
combined_treat$condition<-gsub('.*_','',combined_treat$sample)

ggplot(combined_treat[which(combined_treat$cell_type=='T_cell'),], aes(x = condition, y = proportion)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_point(aes(color = patient), size = 3, position = position_dodge(width = 0.3)) +
  geom_line(aes(group = patient), color = "black", linetype = "dashed", position = position_dodge(width = 0.3)) +
  theme_minimal() +
  labs(x = "Patient", y = "Value", title = "Pre vs Post Treatment Paired Boxplot") +
  scale_color_manual(values = c("Pre" = "#1f77b4", "Post" = "#ff7f0e"))

# Plot
ggplot(combined_treat, aes(x = condition, y = proportion, group = patient)) +
  geom_point(size = 3, aes(color = patient)) +
  geom_line(color = "gray40") +
  theme_minimal() +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(x = "Treatment", y = "Value", title = "Paired Dot Plot (Pre vs Post Treatment)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave('p_treatment.png',p,width = 12,height = 10,dpi = 300)

pre<-function(object) reshape(object, idvar = "patient",timevar = "condition", direction = "wide")
for (i in unique(combined_treat$cell_type)) {
  combined_treat[which(combined_treat$cell_type==i),]->target
  pre(target)->target
  wilcox.test(target$proportion.2,target$proportion.1,paired = F)->a
  print(a$p.value)
}

#single cell
merged <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/merged.rds")
do.call(rbind,result)->result_df
merged$cell_type<-result_df$cell_type[match(colnames(merged),result_df$index)]
DimPlot(merged,group.by = 'cell_type')
DimPlot(merged,label = T)


