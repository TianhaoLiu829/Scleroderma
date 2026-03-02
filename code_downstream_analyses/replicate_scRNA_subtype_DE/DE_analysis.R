library(Seurat, lib.loc = "/ihome/wchen/tianhao/R/x86_64-pc-linux-gnu-library/4.5")
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter_exclude_epi.rds")
#marker and heatmap
merged_filter_base_exclude_epi<-merged_filter_exclude_epi[,grep('1',merged_filter_exclude_epi$orig.ident)]
#find marker
SSC_marker_list_mast<-list()
for (celltype in c('Fibroblasts','Endothelium')) {
  FindMarkers(merged_filter_base_exclude_epi[,which(merged_filter_base_exclude_epi$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',latent.vars = 'orig.ident',test.use = 'MAST',logfc.threshold=0)->SSC_marker_list_mast[[celltype]]
}
saveRDS(SSC_marker_list_mast,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_mast_final.rds')
SSC_marker_list_mast_final<-readRDS('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_mast_final.rds')
#write down the marker genes
library(openxlsx)
wb <- createWorkbook()
for (name in names(SSC_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = SSC_marker_list[[name]])
}
saveWorkbook(wb, file = "multiple_sheets.xlsx", overwrite = TRUE)

#plot the DE dot plot
cell_type<-'Endothelium'
object<-merged_filter_base_exclude_epi[,which(merged_filter_base_exclude_epi$cell_type_l1_all==cell_type)]
features<-rownames(SSC_marker_list_mast_final[[cell_type]])[which(SSC_marker_list_mast_final[[cell_type]]$avg_log2FC>0.5)][1:10]
DotPlot(object,features = features,group.by = 'condition',scale = F)+coord_flip()
ggsave('p_LS_SSC_endo.png',p,width = 6,height = 8,dpi = 300)

VlnPlot(object,features = features,split.by = 'condition')

#get intersect with scRNA data
reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Smooth_muscle_cell", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]

FindMarkers(reference[,which(reference$cell_type_high=='Vascular_endothelium'&grepl('1',reference$orig.ident))],group.by = 'condition',ident.1 = 'SSC',ident.2 = 'LS',latent.vars = 'orig.ident',test.use = 'MAST')->DE_sc_endo
FindMarkers(reference[,which(reference$cell_type_high=='Fibroblasts'&grepl('1',reference$orig.ident))],group.by = 'condition',ident.1 = 'SSC',ident.2 = 'LS',latent.vars = 'orig.ident',test.use = 'MAST')->DE_sc_fib

SSC_marker_list_mast <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_mast_final.rds")
SSC_endo_marker<-rownames(SSC_marker_list_mast$Endothelium)[which(SSC_marker_list_mast$Endothelium$p_val_adj<0.05&SSC_marker_list_mast$Endothelium$avg_log2FC>0)]
DE_sc_endo[SSC_endo_marker,]->DE_sc_endo_inter
DE_sc_endo_inter_SSC<-DE_sc_endo_inter[which(DE_sc_endo_inter$p_val<0.05&DE_sc_endo_inter$avg_log2FC>0),]

LS_endo_marker<-rownames(SSC_marker_list_mast$Endothelium)[which(SSC_marker_list_mast$Endothelium$p_val_adj<0.05&SSC_marker_list_mast$Endothelium$avg_log2FC<0)]
DE_sc_endo[LS_endo_marker,]->DE_sc_endo_inter
DE_sc_endo_inter_LS<-DE_sc_endo_inter[which(DE_sc_endo_inter$p_val<0.05&DE_sc_endo_inter$avg_log2FC<0),]

SSC_fib_marker<-rownames(SSC_marker_list_mast$fibthelium)[which(SSC_marker_list_mast$fibthelium$p_val_adj<0.05&SSC_marker_list_mast$Fibroblasts$avg_log2FC>0)]
DE_sc_fib[SSC_fib_marker,]->DE_sc_fib_inter
DE_sc_fib_inter_SSC<-DE_sc_fib_inter[which(DE_sc_fib_inter$p_val<0.05&DE_sc_fib_inter$avg_log2FC>0),]

LS_fib_marker<-rownames(SSC_marker_list_mast$fibthelium)[which(SSC_marker_list_mast$fibthelium$p_val_adj<0.05&SSC_marker_list_mast$Fibroblasts$avg_log2FC<0)]
DE_sc_fib[LS_fib_marker,]->DE_sc_fib_inter
DE_sc_fib_inter_LS<-DE_sc_fib_inter[which(DE_sc_fib_inter$p_val<0.05&DE_sc_fib_inter$avg_log2FC<0),]



DE_sc_fib_inter_LS <- readRDS("~/DE_sc_fib_inter_LS.rds")
DE_sc_fib_inter_SSC <- readRDS("~/DE_sc_fib_inter_SSC.rds")
DE_sc_endo_inter_LS <- readRDS("~/DE_sc_endo_inter_LS.rds")
DE_sc_endo_inter_SSC <- readRDS("~/DE_sc_endo_inter_SSC.rds")

SSC_fib_marker<-rownames(DE_sc_fib_inter_SSC)
SSC_endo_marker<-rownames(DE_sc_endo_inter_SSC)

LS_fib_marker<-rownames(DE_sc_fib_inter_LS)
LS_endo_marker<-rownames(DE_sc_endo_inter_LS)

merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter_exclude_epi.rds")
object_endo<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium'&grepl('1',merged_filter_exclude_epi$orig.ident))]
dot_plot(object_endo,rev(rownames(DE_sc_endo_inter_SSC)[c(1:15,16,17,21,23,27,38)]),'condition',11)->p
ggsave('p.png',p,width = 3.5,height = 5.5,dpi = 300)
dot_plot(object_endo,rev(rownames(res_endo_inter))[-c(2,6,7,11,16,17)],'condition',11)->p
DotPlot(object_endo,features = rownames(DE_sc_endo_inter_LS),group.by = 'condition')
aggregate(object_endo@assays$RNA@counts[genes,],by=list(object_endo$condition),mean)
ggsave('p.png',p,width = 3.5,height = 4.5,dpi = 300)

object_fib<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts'&grepl('1',merged_filter_exclude_epi$orig.ident))]
dot_plot(object_fib,rev(rownames(DE_sc_fib_inter_LS))[-17],'condition',11)->p
ggsave('p.png',p,width = 3.5,height = 5.5,dpi = 300)
dot_plot(object_fib,rownames(DE_sc_fib_inter_LS)[which(DE_sc_fib_inter_LS$p_val_sc)],'condition',11)->p
ggsave('p.png',p,width = 3.5,height = 5.5,dpi = 300)


####
LS_fib_marker[c(1:8,11,14)]->a
merged_filter_base_exclude_epi_LS_SSC$condition<-factor(merged_filter_base_exclude_epi_LS_SSC$condition,levels = c('LS','SSC'))
DotPlot(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all=='Fibroblasts')],features = a[c(10,9,8,7,6,5,4,3,2,1)],group.by = 'condition',scale = T)+coord_flip()->p

LS_endo_marker[c(1:10)]->a
merged_filter_base_exclude_epi_LS_SSC$condition<-factor(merged_filter_base_exclude_epi_LS_SSC$condition,levels = c('LS','SSC'))
DotPlot(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all=='Endothelium')],features = a[c(10,9,8,7,6,5,4,3,2,1)],group.by = 'condition',scale = T)+coord_flip()->p

ggsave('p.png',p,width = 5,height = 5,dpi = 300)





