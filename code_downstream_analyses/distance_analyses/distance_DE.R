

merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi_update <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi_update.rds")
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds")

library(Seurat)
write_xlsx<-function(file_in,out_dir){
  library(openxlsx)
  wb <- createWorkbook()
  for (name in names(file_in)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = file_in[[name]],rowNames = T)
  }
  saveWorkbook(wb, file = out_dir, overwrite = TRUE)
}
#do the excluding epidermis analysis for all cells
merged_filter_base_exclude_epi_LS_SSC<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident%in%c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))]
merged_filter_base_exclude_epi_LS_SSC<-merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$orig.ident!='G_1')]

SSC_marker_list_wilcox<-list()
for (celltype in unique(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'wilcox')->SSC_marker_list_wilcox[[celltype]]
}

SSC_marker_list_mast<-list()
for (celltype in unique(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',latent.vars = 'orig.ident',test.use = 'MAST')->SSC_marker_list_mast[[celltype]]
}
saveRDS(SSC_marker_list_mast,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_mast.rds')



FeaturePlot(endo_harmony,'LS_endo1')
VlnPlot(endo_harmony,group.by = 'seurat_clusters',features ='LS_endo1')
VlnPlot(endo_harmony,group.by = 'seurat_clusters',features ='nCount_RNA')
VlnPlot(endo_harmony[,which(endo_harmony$seurat_clusters==1)],group.by = 'orig.ident',features ='nCount_RNA')
t.test(endo_harmony$SSC_endo1[which(endo_harmony$seurat_clusters==1)],endo_harmony$SSC_endo1[which(endo_harmony$seurat_clusters!=6)])



cell_type<-'B_Plasma'
object<-merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==cell_type)]
features<-c(rownames(SSC_marker_list_mast[[cell_type]])[which(SSC_marker_list_mast[[cell_type]]$avg_log2FC>0)][1:10],rownames(SSC_marker_list_mast[[cell_type]])[which(SSC_marker_list_mast[[cell_type]]$avg_log2FC<0)][1:10])
DotPlot(object,features = features,group.by = 'condition')+coord_flip()->p
ggsave('p_LS_SSC_endo.png',p,width = 6,height = 8,dpi = 300)

cell_type<-'Fibroblasts'
merged_filter_min_Macrophage_80 <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_min_Macrophage_80.rds")
object<-merged_filter_min_Macrophage_80[,which(merged_filter_min_Macrophage_80$cell_type_l1_all==cell_type)]
SSC_marker_list_distance_mast$min_Macrophage_80[[cell_type]]->mtx
mtx<-mtx[which(mtx$p_val_adj<0.05),]
features<-c(rownames(mtx)[which(mtx$avg_log2FC>0)][1:10],rownames(mtx)[which(mtx$avg_log2FC<0)][1:10])
DotPlot(object,features = features,group.by = 'condition')+coord_flip()->p
ggsave('p_LS_SSC_endo.png',p,width = 6,height = 8,dpi = 300)

SSC_marker_list_negbinom<-list()
for (celltype in unique(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',latent.vars = 'orig.ident',test.use = 'negbinom')->SSC_marker_list_negbinom[[celltype]]
}

SSC_marker_list_poisson<-list()
for (celltype in unique(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',latent.vars = 'orig.ident',test.use = 'poisson')->SSC_marker_list_poisson[[celltype]]
}


fit <- zlm(~ condition + patient, sca)  # all fixed effects
lrt <- lrTest(fit, "conditiondisease")  # test the disease vs control coefficient
tab <- MAST::summary(lrt, doLRT = "conditiondisease")

SSC_marker_list_LR<-list()
for (celltype in unique(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all)) {
  FindMarkers(merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',latent.vars = 'orig.ident',test.use = 'LR')->SSC_marker_list_LR[[celltype]]
}

write_xlsx(SSC_marker_list_mast,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/SSC_LS_exclude_epi.xlsx')

write_xlsx(markers_ccc_all_recei,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/markers_ccc_all_recei.xlsx')
write_xlsx(markers_ccc_all_send,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/markers_ccc_all_send.xlsx')

# see the spatial distribution of those markers
cor_mtx_list<-list()
p_value_list<-list()
for (i in 1:length(SSC_marker_list)) {
  SSC_marker_list[[i]]->mtx
  merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all==names(SSC_marker_list)[i]&merged_filter_exclude_epi$orig.ident!='G_1')]->target
  rownames(mtx)[which(mtx$p_val_adj<0.05&mtx$avg_log2FC>0)]->marker_pos
  rownames(mtx)[which(mtx$p_val_adj<0.05&mtx$avg_log2FC<0)]->marker_neg
  if (length(marker_pos)>1&length(marker_neg)>1) {
    AddModuleScore(target,features = list(marker_pos),name = paste0(names(SSC_marker_list)[i],'_pos') )->target
    AddModuleScore(target,features = list(marker_neg),name = paste0(names(SSC_marker_list)[i],'_neg') )->target
    dis_mtx<-target@meta.data[,grep('min_dis_to',colnames(target@meta.data))]
    gep_mtx<-target@meta.data[,c(paste0(names(SSC_marker_list)[i],'_pos1'),paste0(names(SSC_marker_list)[i],'_neg1'))]
    cor_mtx<-cor(gep_mtx,dis_mtx,method = 'pearson')
    cor_mtx<-as.data.frame(cor_mtx)
    cor_mtx->cor_mtx_list[[i]]
    p_value<-data.frame(matrix(nrow = 2,ncol = 6))
    for (k in 1:6) {
      as.data.frame(cbind(gep_mtx[,c(1,2)],dis_mtx[,k]))->data
      colnames(data) <- c('gep_pos','gep_neg','dis')
      if(sum(data$dis)!=0) {
        summary(lm(gep_pos~dis,data))->model_pos
        summary(lm(gep_neg~dis,data))->model_neg
        model_pos$coefficients[2,4]->p_value[1,k]
        model_neg$coefficients[2,4]->p_value[2,k]
        
      } else {
        NA->p_value[1,k]
        NA->p_value[2,k]
      }
      
    }
    colnames(p_value)<-colnames(dis_mtx)
    rownames(p_value)<-c(paste0(names(SSC_marker_list)[i],'_pos1'),paste0(names(SSC_marker_list)[i],'_neg1'))
    p_value->p_value_list[[i]]
    
  } else if (length(marker_pos)>1&length(marker_neg)<2) {
    AddModuleScore(target,features = list(marker_pos),name = paste0(names(SSC_marker_list)[i],'_pos') )->target
    dis_mtx<-target@meta.data[,grep('min_dis_to',colnames(target@meta.data))]
    gep_mtx<-target@meta.data[,paste0(names(SSC_marker_list)[i],'_pos1')]
    cor_mtx<-cor(gep_mtx,dis_mtx,method = 'pearson')
    cor_mtx<-as.data.frame(cor_mtx)
    rownames(cor_mtx)<-c(paste0(names(SSC_marker_list)[i],'_pos1'))
    cor_mtx->cor_mtx_list[[i]]
    p_value<-data.frame(matrix(nrow = 1,ncol = 6))
    for (k in 1:6) {
      as.data.frame(cbind(gep_mtx,dis_mtx[,k]))->data
      colnames(data) <- c('gep_pos','dis')
      if(sum(data$dis)!=0) {
        summary(lm(gep_pos~dis,data))->model_pos
        model_pos$coefficients[2,4]->p_value[1,k]
      } else {
        NA->p_value[1,k]
      }
      
    }
    colnames(p_value)<-colnames(target@meta.data)[grep('min_dis_to',colnames(target@meta.data))]
    rownames(p_value)<-c(paste0(names(SSC_marker_list)[i],'_pos1'))
    p_value->p_value_list[[i]]
    
  } else if (length(marker_neg)<=1) {
    NA->p_value_list[[i]]
    NA->cor_mtx_list[[i]]
  }
}
names(cor_mtx_list)<-names(SSC_marker_list)
names(p_value_list)<-names(SSC_marker_list)

lapply(1:10,function(x) {
  which(as.matrix(p_value_list[[x]])<0.05)
})

pheatmap::pheatmap(cor_mtx_list$Endothelium[,-2])
pheatmap::pheatmap(cor_mtx_list$Fibroblasts[,-1])
pheatmap::pheatmap(cor_mtx_list$Myeloid[,-4],cluster_rows = FALSE)
pheatmap::pheatmap(cor_mtx_list$T_cell[,-5],cluster_rows = FALSE)
p_value_list$Myeloid

#distance wise or region specific analysis
target<-c('min_T_80','min_Macrophage_80','min_kera_80','min_endo_80','min_fib_80')
SSC_marker_list_distance_wilcox<-list()
for (label in 1:5) {
  readRDS(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_',target[label],'.rds'))->a
  a_base<-a[,which(a$orig.ident%in%c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))]
  SSC_marker_list<-list()
  as.data.frame(table(a_base$cell_type_l1_all,a_base$condition))->tab
  tab<-tab[which(tab$Freq>5),]
  for (celltype in tab$Var1[which(table(tab$Var1)==2)]) {
    FindMarkers(a_base[,which(a_base$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'wilcox')->SSC_marker_list[[celltype]]
  }
  SSC_marker_list->SSC_marker_list_distance_wilcox[[target[label]]]
  write_xlsx(SSC_marker_list,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/wilcox/SSC_LS_exclude_epi',target[label],'.xlsx'))
}
saveRDS(SSC_marker_list_distance_wilcox,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_distance_wilcox.rds')

SSC_marker_list_distance_mast<-list()
for (label in 1:5) {
  readRDS(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_',target[label],'.rds'))->a
  a_base<-a[,which(a$orig.ident%in%c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))]
  SSC_marker_list<-list()
  as.data.frame(table(a_base$cell_type_l1_all,a_base$condition))->tab
  tab<-tab[which(tab$Freq>5),]
  for (celltype in tab$Var1[which(table(tab$Var1)==2)]) {
    FindMarkers(a_base[,which(a_base$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'MAST',latent.vars = 'orig.ident')->SSC_marker_list[[celltype]]
  }
  SSC_marker_list->SSC_marker_list_distance_mast[[target[label]]]
  write_xlsx(SSC_marker_list,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/mast/SSC_LS_exclude_epi',target[label],'.xlsx'))
}
saveRDS(SSC_marker_list_distance_mast,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/SSC_marker_list_distance_mast.rds')



#markers for pre post treatment
merged_filter_treat_exclude_epi<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident%in%c('E_1','H_1','K_1','E_2','H_2','K_2'))]
merged_filter_treat_exclude_epi$pre_post<-'pre'
merged_filter_treat_exclude_epi$pre_post[grep('2',merged_filter_treat_exclude_epi$orig.ident)]<-'post'
merged_filter_treat_exclude_epi$patient<-substr(merged_filter_treat_exclude_epi$orig.ident,1,1)
treat_marker_list<-list()
treat_marker_list_filter<-list()
for (celltype in unique(merged_filter_treat_exclude_epi$cell_type_l1_all)) {
  index_E<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('E_1','E_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_E],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_E
  index_H<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('H_1','H_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_H],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_H
  index_K<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('K_1','K_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_K],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_K
  inte<-intersect(intersect(rownames(treat_marker_K),rownames(treat_marker_H)),rownames(treat_marker_E))
  cbind(treat_marker_E[inte,c(1,2,5)],treat_marker_H[inte,c(1,2,5)],treat_marker_K[inte,c(1,2,5)])->a
  paste0(colnames(a),c(rep('_E_1',3),rep('_H_1',3),rep('_K_1',3)))->colnames(a)
  FindMarkers(merged_filter_treat_exclude_epi,ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'MAST',logfc.threshold=0,min.pct=0.05,latent.vars="patient")->all_mast
  inter<-intersect(rownames(a),rownames(all_mast))
  cbind(a[inter,],all_mast[inter,])->a
  a<-a[order(a$p_val),]
  a->treat_marker_list[[celltype]]
  a<-a[which(a$p_val_adj<0.05&apply(a[,c(2,5,8)],1,function(x){sum(sign(x))%in%c(3,0-3)})),]
  a->treat_marker_list_filter[[celltype]]
  
}

write_xlsx(treat_marker_list_filter,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/differential/gene/post_pre_treat/post_pre_treat_exclude_epi_filter.xlsx'))
write_xlsx(treat_marker_list,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/post_pre_treat/post_pre_treat_exclude_epi.xlsx'))


for (label in 1:5) {
  readRDS(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter_',target[label],'.rds'))->s
  merged_filter_treat_exclude_epi<-s[,which(s$orig.ident%in%c('E_1','H_1','K_1','E_2','H_2','K_2'))]
  merged_filter_treat_exclude_epi$pre_post<-'pre'
  merged_filter_treat_exclude_epi$pre_post[grep('2',merged_filter_treat_exclude_epi$orig.ident)]<-'post'
  merged_filter_treat_exclude_epi$patient<-substr(merged_filter_treat_exclude_epi$orig.ident,1,1)
  treat_marker_list<-list()
  treat_marker_list_filter<-list()
  as.data.frame(table(merged_filter_treat_exclude_epi$orig.ident,merged_filter_treat_exclude_epi$cell_type_l1_all))->tab
  tab<-tab[which(tab$Freq>5),]
  
  for (celltype in names(table(tab$Var2))[which(table(tab$Var2)==6)]) {
    index_E<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('E_1','E_2'))
    FindMarkers(merged_filter_treat_exclude_epi[,index_E],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_E
    index_H<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('H_1','H_2'))
    FindMarkers(merged_filter_treat_exclude_epi[,index_H],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_H
    index_K<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('K_1','K_2'))
    FindMarkers(merged_filter_treat_exclude_epi[,index_K],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0.05)->treat_marker_K
    inte<-intersect(intersect(rownames(treat_marker_K),rownames(treat_marker_H)),rownames(treat_marker_E))
    cbind(treat_marker_E[inte,c(1,2,5)],treat_marker_H[inte,c(1,2,5)],treat_marker_K[inte,c(1,2,5)])->a
    paste0(colnames(a),c(rep('_E_1',3),rep('_H_1',3),rep('_K_1',3)))->colnames(a)
    FindMarkers(merged_filter_treat_exclude_epi,ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'MAST',logfc.threshold=0,min.pct=0.05,latent.vars="patient")->all_mast
    inter<-intersect(rownames(a),rownames(all_mast))
    cbind(a[inter,],all_mast[inter,])->a
    a<-a[order(a$p_val),]
    a->treat_marker_list[[celltype]]
    a<-a[which(a$p_val_adj<0.05&apply(a[,c(2,5,8)],1,function(x){sum(sign(x))%in%c(3,0-3)})),]
    a->treat_marker_list_filter[[celltype]]
    
  }
  write_xlsx(treat_marker_list_filter,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/post_pre_treat/post_pre_treat_exclude_epi_filter_',target[label],'.xlsx'))
  write_xlsx(treat_marker_list,paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/post_pre_treat/post_pre_treat_exclude_epi_',target[label],'.xlsx'))
  
}

#compare the proportion in the neaby cells
labels<-merged_filter_min_endo_80@meta.data[which(merged_filter_min_endo_80$cell_type_l1_all!='Endothelium'),c('cell_type_l1_all','orig.ident','condition')]
library(dplyr)
df_agg <- labels %>%
  group_by(orig.ident, condition, cell_type_l1_all) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n))
library(ggplot2)
df_agg<-df_agg[which(df_agg$orig.ident%in%c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1')),]
df_agg$orig.ident<-factor(df_agg$orig.ident,levels = c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))
ggplot(df_agg, aes(x = orig.ident, y = prop, fill = condition)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ cell_type_l1_all, scales = "free_y") +  # One plot for each cell type
  labs(x = "Sample", y = "Proportion", fill = "Disease Type") +
  theme_minimal()->p




if (nrow(a)>=2) {
  a_meta<-meta(a[,c(2,5,8)],a[,c(1,4,7)],c(length(index_E),length(index_H),length(index_K)))
  a<-cbind(a,a_meta)
} 
a->treat_marker_list[[celltype]]
a<-a[which(apply(a[,c(2,5,8)],1,function(x) {abs(sum(sign(x)))==3} )),]

library(openxlsx)
wb <- createWorkbook()
for (name in names(treat_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = treat_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/treat_exclude_epi_near_fib.xlsx", overwrite = TRUE)




#differential cell signals (SSC vs LS)
SSC_marker_list<-list()
DefaultAssay(merged_filter_min_fib_80)<-'send'
for (celltype in unique(merged_filter_min_fib_80$cell_type_l1_all)) {
  FindMarkers(merged_filter_min_fib_80[,which(merged_filter_min_fib_80$cell_type_l1_all==celltype)],ident.1 = 'SSC',ident.2 = 'LS',group.by = 'condition',test.use = 'wilcox',logfc.threshold=0,min.pct=0)->SSC_marker_list[[celltype]]
}
SSC_marker_list<-lapply(SSC_marker_list, function(x) {x[order(x$p_val),]})
SSC_marker_list<-lapply(SSC_marker_list, function(x) {x[which(x$p_val_adj<0.05),]})
#write down the marker genes
library(openxlsx)
wb <- createWorkbook()
for (name in names(SSC_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = SSC_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/SSC_LS_exclude_epi_near_fib.xlsx", overwrite = TRUE)

#differential cell signals (post vs pre)
merged_filter_treat_exclude_epi<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident%in%c('E_1','H_1','K_1','E_2','H_2','K_2'))]
merged_filter_treat_exclude_epi$pre_post<-'pre'
merged_filter_treat_exclude_epi$pre_post[grep('2',merged_filter_treat_exclude_epi$orig.ident)]<-'post'
DefaultAssay(merged_filter_treat_exclude_epi)<-'send'
treat_marker_list<-list()
treat_marker_list_filter<-list()
for (celltype in unique(merged_filter_treat_exclude_epi$cell_type_l1_all)) {
  index_E<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('E_1','E_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_E],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0)->treat_marker_E
  index_H<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('H_1','H_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_H],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0)->treat_marker_H
  index_K<-which(merged_filter_treat_exclude_epi$cell_type_l1_all==celltype&merged_filter_treat_exclude_epi$orig.ident%in%c('K_1','K_2'))
  FindMarkers(merged_filter_treat_exclude_epi[,index_K],ident.1 = 'post',ident.2 = 'pre',group.by = 'pre_post',test.use = 'wilcox',logfc.threshold=0,min.pct=0)->treat_marker_K
  inte<-intersect(intersect(rownames(treat_marker_K),rownames(treat_marker_H)),rownames(treat_marker_E))
  cbind(treat_marker_E[inte,c(1,2,5)],treat_marker_H[inte,c(1,2,5)],treat_marker_K[inte,c(1,2,5)])->a
  paste0(colnames(a),c(rep('_E_1',3),rep('_H_1',3),rep('_K_1',3)))->colnames(a)
  if (nrow(a)>=2) {
    a_meta<-meta(a[,c(2,5,8)],a[,c(1,4,7)],c(length(index_E),length(index_H),length(index_K)))
    a<-cbind(a,a_meta)
  } 
  a->treat_marker_list[[celltype]]
  a<-a[which(apply(a[,c(2,5,8)],1,function(x) {abs(sum(sign(x)))==3} )),]
  a<-a[which(a$FDR_meta<0.05),]
  dim(a)
  a->treat_marker_list_filter[[celltype]]
}

library(openxlsx)
wb <- createWorkbook()
for (name in names(treat_marker_list)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = treat_marker_list[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/specify_recei_send_cell/treat_exclude_epi_all_send.xlsx", overwrite = TRUE)




#distance DE analysis
merged_filter_exclude_epi$near_fib<-'far'
merged_filter_exclude_epi$near_fib[which(merged_filter_exclude_epi$min_dis_to_fib<80)]<-'near'

#find markers for near and away
marker_near<-list()
for (celltype in unique(merged_filter_exclude_epi$cell_type_l1_all)[-5]) {
  FindMarkers(merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all==celltype)],ident.1='near',ident.2='far',group.by='near_fib')->marker_near[[celltype]]
}
library(openxlsx)
wb <- createWorkbook()
for (name in names(marker_near)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = marker_near[[name]],rowNames = T)
}
saveWorkbook(wb, file = "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/cell_chat/differential/gene/near_vs_far_fib_exclude_epi.xlsx", overwrite = TRUE)


#plot the proportion of <80um
target<-'min_dis_to_Lymphocyte'
exclude<-'T_cell'

target<-'min_dis_to_Myeloid'
exclude<-'Myeloid'

target<-'min_dis_to_fib'
exclude<-'Fibroblasts'

target<-'min_dis_to_endo'
exclude<-'Endothelium'

target<-'min_dis_to_kera'
exclude<-'Keratinocyte'

mtx<-as.data.frame(merged_filter_exclude_epi@meta.data[,c('orig.ident','cell_type_l1_all',target,'condition')])
mtx$near<-'far'
mtx$near[which(mtx[,target]<80)]<-'near'
mtx<-mtx[which(mtx$orig.ident%in%c('E_1','H_1','K_1','F_1','I_1','J_1','M_1')),]
library(dplyr)
df_prop <- mtx %>%
  group_by(cell_type_l1_all, condition, near) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type_l1_all, condition) %>%
  mutate(prop = n / sum(n))
#df_prop_filter<-df_prop[which(df_prop$orig.ident%in%c('E_1','H_1','K_1','F_1','I_1','J_1','M_1')),]
#df_prop_filter$orig.ident<-factor(df_prop_filter$orig.ident,levels = c('E_1','H_1','K_1','F_1','I_1','J_1','M_1'))
df_prop<-df_prop[which(df_prop$cell_type_l1_all!=exclude),]
aggregate(df_prop$prop[which(df_prop$near=='near')],by=list(df_prop$cell_type_l1_all[which(df_prop$near=='near')]),sum)->order_index
df_prop$cell_type_l1_all<-factor(df_prop$cell_type_l1_all,levels = order_index$Group.1[order(0-order_index$x)])

ggplot(df_prop, aes(x = cell_type_l1_all, y = prop, fill = near)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~condition) +   # optional: separate panels for subgroup
  labs(x = "Main Label", y = "Proportion", fill = "Binary Label") +
  theme_minimal()



df_prop <- mtx %>%
  group_by(cell_type_l1_all, orig.ident, near,condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type_l1_all, orig.ident) %>%
  mutate(prop = n / sum(n))
df_prop$orig.ident<-factor(df_prop$orig.ident,levels = c('E_1','H_1','K_1','F_1','I_1','J_1','M_1'))
df_prop<-df_prop[which(df_prop$cell_type_l1_all!=exclude),]
aggregate(df_prop$prop[which(df_prop$near=='near')],by=list(df_prop$cell_type_l1_all[which(df_prop$near=='near')]),sum)->order_index
df_prop$cell_type_l1_all<-factor(df_prop$cell_type_l1_all,levels = order_index$Group.1[order(0-order_index$x)])

df_prop<-df_prop[which(df_prop$near=='near'),]
ggplot(df_prop, aes(x = cell_type_l1_all, y = prop, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.7) +
  theme_minimal()







