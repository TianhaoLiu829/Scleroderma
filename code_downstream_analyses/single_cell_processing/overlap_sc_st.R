library(Seurat)
reference <-readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds")
reference[,which(reference$percent.mt<20&reference$nCount_RNA<20000&reference$nCount_RNA>200)]->reference

Idents(reference)<-reference$cell_type_l1
FindAllMarkers(reference,only.pos = T)->reference_marker
saveRDS(reference_marker,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/single_cell_DE/reference_marker.rds')

Idents(merged_filter)<-merged_filter$cell_type_l1_all
FindAllMarkers(merged_filter,only.pos = T)->merged_filter_marker
saveRDS(merged_filter_marker,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/single_cell_DE/merged_filter_marker.rds')

-asd-asd

merged_filter_marker <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/single_cell_DE/merged_filter_marker.rds")
reference_marker <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/single_cell_DE/reference_marker.rds")
reference_marker<-reference_marker[which(reference_marker$gene%in%rownames(merged_filter)),]
reference_marker$p_val_adj<-p.adjust(reference_marker$p_val,method = 'bonferroni')
reference_marker<-reference_marker[which(reference_marker$p_val_adj<0.01),]
reference_marker<-reference_marker[which(reference_marker$gene!='CCL19'),]
merged_filter_marker<-merged_filter_marker[which(merged_filter_marker$p_val_adj<0.01),]
merged_filter_marker<-merged_filter_marker[which(merged_filter_marker$gene!='CCL19'),]

p_pvalue_list<-list()
p_enrich_list<-list()
for (thresh in seq(5,300,5)) {
  thresh<-20
  p_mtx<-data.frame(row.names = unique(reference_marker$cluster))
  p_mtx->p_mtx_fisher
  p_mtx->p_mtx_ratio
  for (singlecell in unique(reference_marker$cluster)) {
    for (spatial in unique(merged_filter_marker$cluster)) {
      reference_marker$gene[which(reference_marker$cluster==singlecell)]->sc_marker
      sc_marker<-sc_marker[which(!sc_marker%in%reference_marker$gene[which(reference_marker$cluster!=singlecell)])][1:thresh]
      merged_filter_marker$gene[which(merged_filter_marker$cluster==spatial)]->st_marker
      st_marker<-st_marker[which(!st_marker%in%merged_filter_marker$gene[which(merged_filter_marker$cluster!=spatial)])][1:thresh]
      length(intersect(sc_marker,st_marker))->k
      length(sc_marker)->n1
      18085->N
      length(st_marker)->n2
      (k/n2)/(n1/N)->p_mtx_ratio[singlecell,spatial]
      phyper(k-1,n1,N-n1,n2,lower.tail = FALSE)->p_mtx[singlecell,spatial]
      mat <- matrix(c(k, n1 - k, n2 - k, N - n1 - n2 + k), nrow = 2)
      p_mtx_fisher[singlecell,spatial] <- fisher.test(mat, alternative = "greater")$p.value
    }
  }
  #mtx and mtx_fisher are equivalent
  p_mtx <- matrix(
    p.adjust(as.vector(as.matrix(p_mtx)), method = "bonferroni"),
    nrow = nrow(p_mtx),
    ncol = ncol(p_mtx),
    dimnames = dimnames(p_mtx)
  )
  p_mtx<-0-log10(p_mtx)
  pheatmap::pheatmap(p_mtx)->p_pvalue_list[[as.character(thresh)]]
  as.matrix(p_mtx_ratio)->p_mtx_ratio
  p_mtx_ratio[which(p_mtx_ratio>200)]<-200
  pheatmap::pheatmap(p_mtx_ratio)->p
  ggsave('p.png',p,width=7,height = 7,dpi = 300)
  pheatmap::pheatmap(log2(p_mtx_ratio+1))->p_enrich_list[[as.character(thresh)]]
  pheatmap::pheatmap(log2(p_mtx_ratio+1))->p
  log2(p_mtx_ratio+1)[p$tree_row$order[c(2:7,1,8:10)],p$tree_col$order]->mtx
  png("/ihome/wchen/tianhao/heatmap.png", width = 2800, height = 2400, res = 300)
  pheatmap::pheatmap(mtx,cluster_rows = F,cluster_cols = F,fontsize = 17)->p
  dev.off()
  
  ggsave('p.png',p,width=7,height=6,dpi=300)
}

cor_mtx<-data.frame(row.names = unique(reference_marker$cluster))
for (singlecell in unique(reference_marker$cluster)) {
  for (spatial in unique(merged_filter_marker$cluster)) {
    reference_marker$gene[which(reference_marker$cluster==singlecell)][1:10]->sc_marker
    intersect(sc_marker,rownames(merged_filter))->share
    mean_marker<-data.frame(st_marker=rowMeans(as.matrix(merged_filter@assays$RNA@data[share,which(merged_filter$cell_type_l1_all==spatial)])),sc_marker=rowMeans(as.matrix(reference@assays$RNA@data[share,which(reference$cell_type_l1==singlecell)])))
    cor(mean_marker)[2,1]->cor_mtx[singlecell,spatial]
  }
}

pheatmap::pheatmap(cor_mtx)->p
ggsave('p.png',p,width=8,height=8,dpi=300)

cor_mtx<-data.frame(row.names = unique(reference_marker$cluster))
for (singlecell in unique(reference_marker$cluster)) {
  for (spatial in unique(merged_filter_marker$cluster)) {
    reference_marker$gene[which(reference_marker$cluster==singlecell)]->sc_marker
    sc_marker<-sc_marker[which(!sc_marker%in%reference_marker$gene[which(reference_marker$cluster!=singlecell)])][1:30]
    merged_filter_marker$gene[which(merged_filter_marker$cluster==spatial)]->st_marker
    st_marker<-st_marker[which(!st_marker%in%merged_filter_marker$gene[which(merged_filter_marker$cluster!=spatial)])][1:30]
    intersect(sc_marker,st_marker)->share
    length(share)->k
    if (k>=6) {
      mean_marker<-data.frame(st_marker=rowMeans(as.matrix(merged_filter@assays$RNA@data[share,which(merged_filter$cell_type_l1_all==spatial)])),sc_marker=rowMeans(as.matrix(reference@assays$RNA@data[share,which(reference$cell_type_l1==singlecell)])))
      cor(mean_marker)[2,1]->cor_mtx[singlecell,spatial]
    } else {
      cor_mtx[singlecell,spatial]<-0
    }
  }
}
pheatmap(cor_mtx)

#correlate cell proportions in ST and SC
library(dplyr)
library(ggplot2)
library(scales)
reference <-readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds")
merged_filter$cell_type_l1_all[which(merged_filter$cell_type_l1_all%in%c('Macrophage','Neutrophil'))]<-'Myeloid'

prop_ssc_base_sc<-as.data.frame(prop.table(table(reference$cell_type_l1[which(grepl('1',reference$orig.ident)&reference$condition=='SSC')])))
prop_ls_base_sc<-as.data.frame(prop.table(table(reference$cell_type_l1[which(grepl('1',reference$orig.ident)&reference$condition=='LS')])))
prop_ssc_post_sc<-as.data.frame(prop.table(table(reference$cell_type_l1[which(grepl('2',reference$orig.ident)&reference$condition=='SSC')]))) 

prop_ssc_base_st<-as.data.frame(prop.table(table(merged_filter$cell_type_l1_all[which(grepl('1',merged_filter$orig.ident)&merged_filter$condition=='SSC')])))
prop_ls_base_st<-as.data.frame(prop.table(table(merged_filter$cell_type_l1_all[which(grepl('1',merged_filter$orig.ident)&merged_filter$condition=='LS')])))
prop_ssc_post_st<-as.data.frame(prop.table(table(merged_filter$cell_type_l1_all[which(grepl('2',merged_filter$orig.ident)&merged_filter$condition=='SSC')]))) 

colors <- c(
  "Fibroblasts"="#00BBFF",  # deep teal
  "Endothelium"="#D95F02",  # orange
  "Glandular_epithelium"="#7570B3",  # indigo
  "Keratinocyte"="#8DD3C7",  # cyan (replaces gray)
  "Neuron"="#E7298A",  # magenta
  "Smooth_muscle_cell"="#A6761D",  # brown
  "Melanocyte"="#B2DF8A",  # green
  "B_Plasma"="#A6CEE3",  # sky blue
  "Myeloid"="#E6AB02",  # gold
  "T_cell"="#66A61E"
)

as.data.frame(rbind(prop_ssc_base_sc,prop_ssc_base_st))->plot_df
as.data.frame(rbind(prop_ls_base_sc,prop_ls_base_st))->plot_df

as.data.frame(rbind(prop_ssc_base_st,prop_ssc_post_st))->plot_df
as.data.frame(rbind(prop_ssc_base_sc,prop_ssc_post_sc))->plot_df

as.data.frame(rbind(prop_ssc_base_st,prop_ls_base_st))->plot_df
as.data.frame(rbind(prop_ssc_base_sc,prop_ls_base_sc))->plot_df

level<-c("Myeloid","T_cell","B_Plasma","Fibroblasts","Endothelium","Keratinocyte","Melanocyte","Glandular_epithelium","Smooth_muscle_cell" ,"Neuron"  )
plot_df$Var1<-factor(plot_df$Var1,levels = level[seq(10,1,-1)])
plot_df$group<-c(rep('single_cell',10),rep('spatial',10))
ggplot(plot_df, aes(x = group, y = Freq, fill = Var1)) +
  geom_col(width = 0.8, color = NA) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  scale_fill_manual(
    values = colors,
    drop = FALSE        # keep same legend/order even if missing
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    legend.position = "right"
  )->p
ggsave('p.png',p,width = 10,height = 5,dpi = 300)


prop_all_sc<-as.data.frame(prop.table(table(reference$cell_type_l1)))
prop_all_st<-as.data.frame(prop.table(table(merged_filter$cell_type_l1_all)))
#prop_all_sc<-as.data.frame(prop.table(table(reference$cell_type_l1[which(reference$orig.ident%in%c('G2','H2','M1','I1'))])))
#prop_all_st<-as.data.frame(prop.table(table(merged_filter$cell_type_l1_all[which(merged_filter$orig.ident%in%c('G_2','H_2','M_1','I_1'))])))
as.data.frame(cbind(prop_all_sc,prop_all_st[-1]))->df
colnames(df)<-c('group','sc','st')
cor.test(df[,2],df[,3],method = 'spearman')
library(ggrepel)
df[5,c('sc','st')]-0.35->df[5,c('sc','st')]
ggplot(df, aes(x =sc, y = st, color = group)) +
  geom_point(size = 1.5) +                              # scatter dots
  geom_abline(slope = 1, intercept = 0,               # y = x line
              linetype = "dashed", size = 1) + 
  scale_color_manual(values = colors) +  
  theme_classic() +
  scale_x_continuous(limits = c(0,0.25)) +
  scale_y_continuous(limits = c(0,0.25)) +
  theme(axis.text = element_text(size = 15))+
  geom_text_repel(aes(label = group), size = 4)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  labs(x = "SC", y = "ST", color = "Group")->p
ggsave('p1.png',p,width=9,height = 7)

#bin distance
library(dplyr)
merged_filter_exclude_epi<- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident!='G_1')]

df_prop <- as.data.frame(merged_filter_exclude_epi@meta.data) %>%
  mutate(bin = cut(min_dis_to_fib, breaks = seq(0, max(min_dis_to_fib), by = 10))) %>%
  group_by(cell_type_l1_all, bin) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(total = sum(count)) %>%          # total per cell type
  mutate(proportion = count / total) %>%  # fraction within bin
  ungroup()

library(ggplot2)

ggplot(as.data.frame(merged_filter_exclude_epi@meta.data[which(merged_filter_exclude_epi@meta.data$min_dis_to_fib<500&merged_filter_exclude_epi@meta.data$cell_type_l1_all=='Keratinocyte'),]), aes(x = min_dis_to_fib, color = cell_type_l1_all, fill = cell_type_l1_all)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 80, color = "red")+
  labs(x = "Distance", y = "Proportion (density)") +
  theme_minimal()->p

merged_filter_exclude_epi$fib<-FALSE
merged_filter_exclude_epi$fib[which(merged_filter_exclude_epi$min_dis_to_fib<80)]<-TRUE

