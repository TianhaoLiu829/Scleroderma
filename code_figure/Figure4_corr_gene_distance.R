endo_marker_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/endo_marker_LS.rds")
DE_sc_endo_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_endo_inter_SSC.rds")
DE_sc_fib_inter_LS <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_LS.rds")
DE_sc_fib_inter_SSC <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/DE_sc_fib_inter_SSC.rds")


merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")

library(Seurat)
AddModuleScore(merged_filter_exclude_epi,features = list(rownames(endo_marker_LS)),name = 'endo_LS_marker')->merged_filter_exclude_epi
AddModuleScore(merged_filter_exclude_epi,features = list(rownames(DE_sc_fib_inter_LS)),name = 'fib_LS_marker')->merged_filter_exclude_epi
AddModuleScore(merged_filter_exclude_epi,features = list(rownames(DE_sc_fib_inter_SSC)),name = 'fib_SSC_marker')->merged_filter_exclude_epi
AddModuleScore(merged_filter_exclude_epi,features = list(rownames(DE_sc_endo_inter_SSC)),name = 'endo_SSC_marker')->merged_filter_exclude_epi



endo<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium'&grepl('1',merged_filter_exclude_epi$orig.ident))]
endo_LS<-endo[,which(endo$condition=='LS')]
endo_SSC<-endo[,which(endo$condition=='SSC')]
cor.test(endo_LS$endo_LS_marker1,endo_LS$min_dis_to_Macrophage,method = 'spearman')->a1
cor.test(endo_LS$endo_LS_marker1,endo_LS$min_dis_to_T,method = 'spearman')->a2
cor.test(endo_SSC$endo_SSC_marker1,endo_SSC$min_dis_to_Macrophage,method = 'spearman')->a3
cor.test(endo_SSC$endo_SSC_marker1,endo_SSC$min_dis_to_T,method = 'spearman')->a4
matrix(c(a1$estimate,a3$estimate,a2$estimate,a4$estimate),nrow = 2,ncol = 2)->mtx_endo
breaks <- c(
  seq(-max(abs(mtx_endo))-0.01,0, length.out = 51),
  seq(0,max(abs(mtx_endo)+0.01), length.out = 51)[-1]
)

pheatmap::pheatmap(0-mtx_endo,cluster_rows = F,cluster_cols = F,breaks =breaks,color =  colorRampPalette(c("#B40426","#FFFFFF","#3B4CC0"))(101))->p
pheatmap::pheatmap(0-mtx_endo,cluster_rows = F,cluster_cols = F,breaks =breaks)->p
matrix(c(a1$p.value,a3$p.value,a2$p.value,a4$p.value),nrow = 2,ncol = 2)->mtx_endo_p
ggsave('p_endo.png',p,width = 4,height = 3,dpi = 300)
data.frame(GEP=endo_LS$endo_LS_marker1,dis=endo_LS$min_dis_to_Macrophage/4)->df_endo
for (thresh in c(100,200,300,400,500,600,700,800,10^10)) {
  cor.test(df_endo$GEP[which(df_endo$dis<thresh)],df_endo$dis[which(df_endo$dis<thresh)])->a
  print(a$estimate)
  print(a$p.value)
}
ggplot(df_endo[which(df_endo$dis<400),], aes(x =dis, y = GEP)) +
  geom_point(color = "black",alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm",se = T,color = "red") +
  #geom_vline(xintercept = 33,linetype='dashed',linewidth = 1,color='cyan') +
  theme_classic()+
  theme(axis.text = element_text(size = 16))->p
ggsave('p_endo_cor.png',p,width = 5,height = 5)



fib<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts'&grepl('1',merged_filter_exclude_epi$orig.ident))]
fib_LS<-fib[,which(fib$condition=='LS')]
fib_SSC<-fib[,which(fib$condition=='SSC')]
cor.test(fib_LS$fib_LS_marker1,fib_LS$min_dis_to_Macrophage,method = 'spearman')->a1
cor.test(fib_LS$fib_LS_marker1,fib_LS$min_dis_to_T,method = 'spearman')->a2
cor.test(fib_SSC$fib_SSC_marker1,fib_SSC$min_dis_to_Macrophage,method = 'spearman')->a3
cor.test(fib_SSC$fib_SSC_marker1,fib_SSC$min_dis_to_T,method = 'spearman')->a4
matrix(c(a1$estimate,a3$estimate,a2$estimate,a4$estimate),nrow = 2,ncol = 2)->mtx_fib
breaks <- c(
  seq(-max(abs(mtx_fib))-0.01,0, length.out = 51),
  seq(0,max(abs(mtx_fib)+0.01), length.out = 51)[-1]
)
pheatmap::pheatmap(0-mtx_fib,cluster_rows = F,cluster_cols = F,breaks = breaks)->p
matrix(c(a1$p.value,a3$p.value,a2$p.value,a4$p.value),nrow = 2,ncol = 2)->mtx_fib_p
ggsave('p_fib.png',p,width = 4,height = 3,dpi = 300)
data.frame(GEP=fib_LS$fib_LS_marker1,dis=fib_LS$min_dis_to_Macrophage/4)->df_fib
for (thresh in c(100,200,300,400,500,600,700,800)) {
  cor.test(df_fib$GEP[which(df_fib$dis<thresh)],df_fib$dis[which(df_fib$dis<thresh)])->a
  print(a$estimate)
  print(a$p.value)
}
ggplot(df_fib[which(df_fib$dis<400),], aes(x =dis, y = GEP)) +
  geom_point(color = "black",alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm",se = T,color = "red") +
  #geom_vline(xintercept = 33,linetype='dashed',linewidth = 1,color='cyan') +
  theme_classic()+
  theme(axis.text = element_text(size = 16))->p
ggsave('p_fib_cor.png',p,width = 5,height = 5)
