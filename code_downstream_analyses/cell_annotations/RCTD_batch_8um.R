#RCTD deconvolution
.libPaths(c("/ihome/wchen/tianhao/R/x86_64-pc-linux-gnu-library/4.5","/software/rhel9/manual/install/r/4.5.0/lib64/R/library"))
library(spacexr)
.libPaths(c("/software/rhel9/manual/install/r/4.5.0/lib64/R/library","/ihome/wchen/tianhao/R/x86_64-pc-linux-gnu-library/4.5"))
list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/')->sample
sample<-sample[which(sample!='crc')]
library(Seurat)

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
#merged <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/merged.rds")
#read the 8um spots with UMI>10
merge <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spot_level/all_processed_8um_spot_thresh_50.rds")
substr(colnames(merge),1,2)->merge$orig.ident
merge$condition[which(merge$condition=='HC')]<-'LS'

RCTD_DECON<-function(object,target,ref){
  object[,which(object$orig.ident==target)]->object
  ref<-ref[,which(ref$condition==unique(object$condition))]
  cluster <- as.factor(ref$cell_type_l1)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  refe <- Reference(ref@assays[["RNA"]]@counts, cluster, nUMI)
  counts <-object@assays$RNA@counts
  rownames(counts)<-rownames(object)
  colnames(counts)<-colnames(object)
  
  coords <- read.csv(paste0("/ix1/wchen/liutianhao/result/pathology_ST/position_8um/",target,"_positions.csv"))
  coords<-coords[which(coords$in_tissue==1),]
  rownames(coords)<-paste0(target,'_',coords$barcode)
  coords<-coords[,c('pxl_row_in_fullres','pxl_col_in_fullres')]
  coords<-coords[colnames(counts),]
  query <- SpatialRNA(coords, counts, colSums(counts))
  RCTD <- create.RCTD(query, refe, max_cores = 40,UMI_min = 50,UMI_min_sigma = 100,CELL_MIN_INSTANCE=5)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  object <- AddMetaData(object, metadata = RCTD@results$results_df)
  object<-object[,rownames(RCTD@results[["weights"]])]
  object[["conv"]] <- CreateAssayObject(counts =Transpose(RCTD@results[["weights"]]))
  coords<-coords[colnames(object),]
  cbind(colnames(object),coords[,c('pxl_row_in_fullres','pxl_col_in_fullres')],t(as.matrix(object[['conv']]$counts)),object$first_type)->result
  colnames(result)[c(1,dim(result)[2])]<-c('index','cell_type')
  write.csv(result,paste0('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/spot_level/l1/result_',target,'.csv'),row.names = FALSE)
  return(object)
}

library(SeuratData)
library(SeuratDisk)


#lapply(sample[which(nchar(sample)==3)][c(9,1,2,3,4,5,6,7,8,10)], function(x) {RCTD_DECON(merged,x,reference,5)})->result
lapply(gsub('_','',sample[which(nchar(sample)==3)]) , function(x) {RCTD_DECON(merge,x,reference)})->result
#names(result)<-sample[1:10]
names(result)<-c(gsub('_','',sample[which(nchar(sample)==3)]))
saveRDS(result,'/ix1/wchen/liutianhao/result_l1_8um.rds')
-asd-asd-asd
spot_merged <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spot_level/all_processed_8um_spot_thresh_50.rds")
result<-list()
for(i in gsub('_','',sample[which(nchar(sample)==3)])[c(1,3:13,15)]){
  read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/spot_level/l1/result_',i,'.csv'))->result[[i]]
} 
do.call(rbind,result)->result
spot_merged$cell_type_l1_all<-result$cell_type[match(colnames(spot_merged),result$index)]
DimPlot(spot_merged[,which(spot_merged$orig.ident!='L1')],group.by = 'cell_type_l1_all')
sc_pro(spot_merged[,which(spot_merged$orig.ident!='L1')],FALSE)->spot_spot_merged_filtereded
DimPlot(spot_spot_merged_filtereded,group.by = 'cell_type_l1_all')
saveRDS(spot_spot_merged_filtereded,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spot_level/all_processed_exclude_L1_thresh_50.rds')

#find marker genes (LS, SSC)
spot_merged_filtered <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spot_level/all_processed_exclude_L1_thresh_50.rds")
spot_merged_filtered_base<-spot_merged_filtered[,-grep('2',spot_merged_filtered$orig.ident)]
spot_merged_filtered_base@meta.data<-as.data.frame(cbind(spot_merged_filtered_base@meta.data,result[match(colnames(spot_merged_filtered_base),result$index),-c(1,2,3,14)]))
SSC_marker_list<-list()
for (celltype in unique(spot_merged_filtered_base$cell_type_l1_all)) {
  FindMarkers(spot_merged_filtered_base[,which(spot_merged_filtered_base$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'LR',latent.vars = colnames(spot_merged_filtered_base@meta.data)[9:18])->SSC_marker_list[[celltype]]
}
#write down the marker genes
library(openxlsx)
wb <- createWorkbook()
for (name in names(SSC_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = SSC_marker_list[[name]])
}
saveWorkbook(wb, file = "multiple_sheets_spot.xlsx", overwrite = TRUE)
#find cell type markers
Idents(spot_merged_filtered_base)<-spot_merged_filtered_base$cell_type_l1_all
FindAllMarkers(spot_merged_filtered_base,only.pos = T)->marker_cell_type
df_list <- split(marker_cell_type, marker_cell_type$cluster)
wb <- createWorkbook()
for (label in names(df_list)) {
  addWorksheet(wb, sheetName = label)
  writeData(wb, sheet = label, x = df_list[[label]])
}
saveWorkbook(wb, "split_by_cell_type_spot.xlsx", overwrite = TRUE)


#markers for treatment
#markers for pre post treatment
merged_filter<-readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spot_level/all_processed_exclude_L1_thresh_50.rds")
merged_filter_treat<-merged_filter[,which(merged_filter$orig.ident%in%c('E1','E2','K1','K2','H1','H2')) ]
sc_pro(merged_filter_treat,FALSE)->merged_filter_treat
merged_filter_treat_exclude_epi<-merged_filter_treat[,which(!merged_filter_treat$epidermis)]
merged_filter_treat_exclude_epi@meta.data<-as.data.frame(cbind(merged_filter_treat_exclude_epi@meta.data,result[match(colnames(merged_filter_treat_exclude_epi),result$index),-c(1,2,3,14)]))
sc_pro(merged_filter_treat_exclude_epi,FALSE)->merged_filter_treat_exclude_epi

merged_filter_treat_exclude_epi$pre_post<-'pre'
merged_filter_treat_exclude_epi$pre_post[grep('2',merged_filter_treat_exclude_epi$orig.ident)]<-'post'
merged_filter_treat$pre_post<-'pre'
merged_filter_treat$pre_post[grep('2',merged_filter_treat$orig.ident)]<-'post'

treat_marker_list<-list()
for (celltype in unique(merged_filter_treat_exclude_epi$cell_type_l1_all)[-2]) {
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('E1','E2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_E
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('H1','H2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_H
  FindMarkers(merged_filter_treat_exclude_epi[,which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('K1','K2'))],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'LR',latent.vars = colnames(merged_filter_treat_exclude_epi@meta.data)[10:19])->treat_marker_K
  inte<-intersect(intersect(rownames(treat_marker_K),rownames(treat_marker_H)),rownames(treat_marker_E))
  cbind(treat_marker_E[inte,c(1,2,5)],treat_marker_H[inte,c(1,2,5)],treat_marker_K[inte,c(1,2,5)])->treat_marker_list[[celltype]]
}
library(openxlsx)
wb <- createWorkbook()
for (name in names(treat_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = treat_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "treat_exclude_epi_adj.xlsx", overwrite = TRUE)

#HEATMAP
library(Seurat)
DoHeatmap(all_processed_exclude_L1_thresh_50_filter_epi,group.by = 'cell_type_l1_all',features = top10$gene)
all_processed_exclude_L1_thresh_50_filter_epi@assays$RNA@var.features<-union(all_processed_exclude_L1_thresh_50_filter_epi@assays$RNA@var.features,top10$gene)
all_processed_exclude_L1_thresh_50_filter_epi<-ScaleData(all_processed_exclude_L1_thresh_50_filter_epi)
all_processed_exclude_L1_thresh_50_filter_epi$cell_type_l1_all<-factor(all_processed_exclude_L1_thresh_50_filter_epi$cell_type_l1_all,levels = levels(top10$cluster))
DoHeatmap(all_processed_exclude_L1_thresh_50_filter_epi[,which(!is.na(all_processed_exclude_L1_thresh_50_filter_epi$cell_type_l1_all))],group.by = 'cell_type_l1_all',features = top10$gene)->p
ggsave('p.png',p,width = 14,height = 9,dpi = 300)

