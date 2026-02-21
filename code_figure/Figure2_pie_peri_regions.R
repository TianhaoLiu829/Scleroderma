merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
test_dist<-function(cell_type){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  lapply(dist_matrix_all, function(x) {
    as.data.frame(apply(x[,colnames(x)%in%names],1,min))
  })->all
  do.call(rbind,all)->all
  return(all[match(colnames(merged_filter_exclude_epi),gsub('.*[.]','',rownames(all))),1])
}
test_dist('Glandular_epithelium')->merged_filter_exclude_epi$min_dis_to_gland


merged_filter_exclude_epi<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident!='G_1')]
thresh<-100
margin_list<-list()
for (celltype in c('Macrophage','T_cell','B_Plasma')) {
  merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all==celltype)]->macro
  as.data.frame(prop.table(table(macro$orig.ident,macro$min_dis_to_fib<thresh),margin = 1))->margin_fib
  as.data.frame(prop.table(table(macro$orig.ident,macro$min_dis_to_endo<thresh),margin = 1))->margin_endo
  margin_fib<-margin_fib[which(margin_fib$Var2==TRUE),]
  margin_endo<-margin_endo[which(margin_endo$Var2==TRUE),]
  margin_fib$peri<-'fib'
  margin_endo$peri<-'endo'
  as.data.frame(rbind(margin_endo,margin_fib))->margin
  margin$celltype<-celltype
  margin->margin_list[[celltype]]
}
do.call(rbind,margin_list)->margin_all
t.test(margin_all$Freq[which(margin_all$celltype=='T_cell'),3],margin_all[which(margin_all$celltype=='T_cell'),6])

library(tidyverse)
df_summary <- margin_all %>%
  group_by(peri, celltype) %>%
  summarise(
    mean_prop = mean(Freq),
    sd_prop = sd(Freq),
    se_prop = sd(Freq) / sqrt(n()),
    .groups = "drop"
  )
ggplot(df_summary, aes(x = peri, y = mean_prop, fill = celltype)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Proportion (%)", fill = "Group") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))



margin_list<-list()
for (dis in c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')) {
  merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts')]->fib
  merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')]->endo
  as.data.frame(prop.table(table(fib$orig.ident,fib@meta.data[,dis]<thresh),margin = 1))->margin_fib
  as.data.frame(prop.table(table(endo$orig.ident,endo@meta.data[,dis]<thresh),margin = 1))->margin_endo
  margin_fib<-margin_fib[which(margin_fib$Var2==TRUE),]
  margin_endo<-margin_endo[which(margin_endo$Var2==TRUE),]
  margin_fib$peri<-'fib'
  margin_endo$peri<-'endo'
  as.data.frame(rbind(margin_endo,margin_fib))->margin
  margin$dis<-dis
  margin->margin_list[[dis]]
}
do.call(rbind,margin_list)->margin_all
t.test(margin_all$Freq[which(margin_all$peri=='endo'&margin_all$dis=='min_dis_to_Macrophage')],margin_all$Freq[which(margin_all$peri=='fib'&margin_all$dis=='min_dis_to_Macrophage')],paired = T)
t.test(margin_all$Freq[which(margin_all$peri=='endo'&margin_all$dis=='min_dis_to_T')],margin_all$Freq[which(margin_all$peri=='fib'&margin_all$dis=='min_dis_to_T')],paired = T)
t.test(margin_all$Freq[which(margin_all$peri=='endo'&margin_all$dis=='min_dis_to_B_Plasma')],margin_all$Freq[which(margin_all$peri=='fib'&margin_all$dis=='min_dis_to_B_Plasma')],paired = T)


