
reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
readarray<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Fibroblasts", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]

#reference<-reference[,which(reference$cell_type_l1!='Keratinocyte')]
#merged_filter<-merged_filter[,which(merged_filter$cell_type_l1_all!='Keratinocyte')]

sc_SSC<-as.data.frame(table(reference$cell_type_l1[which(reference$orig.ident%in%c('A1','E1','H1','K1'))])/sum(table(reference$cell_type_l1[which(reference$orig.ident%in%c('A1','E1','H1','K1'))])))
sc_LS<-as.data.frame(table(reference$cell_type_l1[which(reference$orig.ident%in%c('B1','C1','D1','F1','G1','I1','J1','M1'))])/sum(table(reference$cell_type_l1[which(reference$orig.ident%in%c('B1','C1','D1','F1','G1','I1','J1','M1'))])))
sc_all<-as.data.frame(table(reference$cell_type_l1)/sum(table(reference$cell_type_l1)))
st_SSC<-as.data.frame(table(merged_filter$cell_type_l1_all[which(merged_filter$orig.ident%in%c('E_1','H_1','K_1'))])/sum(table(merged_filter$cell_type_l1_all[which(merged_filter$orig.ident%in%c('E_1','H_1','K_1'))])))
st_LS<-as.data.frame(table(merged_filter$cell_type_l1_all[which(merged_filter$orig.ident%in%c('F_1','G_1','I_1','J_1','M_1'))])/sum(table(merged_filter$cell_type_l1_all[which(merged_filter$orig.ident%in%c('F_1','G_1','I_1','J_1','M_1'))])))
st_all<-as.data.frame(table(merged_filter$cell_type_l1_all)/sum(table(merged_filter$cell_type_l1_all)))


sum(st_SSC[which(st_SSC$Var1 %in% c('Macrophage','Neutrophil')),2])->st_SSC[which(st_SSC$Var1==c('Macrophage')),2]
sum(st_LS[which(st_LS$Var1 %in% c('Macrophage','Neutrophil')),2])->st_LS[which(st_LS$Var1==c('Macrophage')),2]
sum(st_all[which(st_all$Var1 %in% c('Macrophage','Neutrophil')),2])->st_all[which(st_all$Var1==c('Macrophage')),2]
st_SSC<-st_SSC[-which(st_SSC$Var1==c('Neutrophil')),]
st_LS<-st_LS[-which(st_LS$Var1==c('Neutrophil')),]
st_all<-st_all[-which(st_all$Var1==c('Neutrophil')),]
st_SSC[,1]<-as.character(st_SSC[,1])
st_LS[,1]<-as.character(st_LS[,1])
st_all[,1]<-as.character(st_all[,1])
st_SSC[which(st_SSC$Var1==c('Macrophage')),1]<-'Myeloid'
st_LS[which(st_LS$Var1==c('Macrophage')),1]<-"Myeloid"
st_all[which(st_all$Var1==c('Macrophage')),1]<-"Myeloid"
cbind(sc_SSC[,2],sc_LS[,2],st_SSC[,2],st_LS[,2])->all_sep
cbind(sc_all[,2],st_all[,2])->all_all


colnames(all_sep)<-c('sc_SSC','sc_LS','st_SSC','st_LS')
colnames(all_all)<-c('sc','st')
rownames(all_sep)<-sc_SSC[,1]
rownames(all_all)<-sc_all[,1]
all_sep<-as.data.frame(all_sep)
all_all<-as.data.frame(all_all)

#correlation plot
rownames(all_all)[6]<-'T_cell'
prop_table$group<-rownames(prop_table)
library(ggbreak)
prop_table['Keratinocyte',c(1,2)]<-prop_table['Keratinocyte',c(1,2)]-0.35
all_all$group<-rownames(all_all)
ggplot(all_all, aes(x =sc, y = st, color = group)) +
  geom_point(size = 1.5) +                              # scatter dots
  geom_abline(slope = 1, intercept = 0,               # y = x line
              linetype = "dashed", size = 1) + 
  scale_color_manual(values = colors) +  
  theme_classic() +
  theme(axis.text = element_text(size = 15))+
  geom_text_repel(aes(label = group), size = 4)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  labs(x = "X value", y = "Y value", color = "Group")->p
ggsave('p.png',p,width = 7.5,height = 6,dpi = 300)

inter<-intersect(unique(reference$orig.ident),gsub('_','',unique(merged_filter_exclude_epi$orig.ident)))
reference$cell_type_l1[which(reference$cell_type_l1=='Myeloid')]<-'Macrophage'
reference$cell_type_l1[which(reference$cell_type_l1=='Lymphocyte')]<-'T_cell'

prop_table_list<-list()
for (sample in inter) {
  data.frame(matrix(nrow = 10,ncol = 2))->prop_table
  rownames(prop_table)<-unique(reference@meta.data$cell_type_l1)
  colnames(prop_table)<-c('ST','SC')
  paste0(substr(sample,1,1),'_',substr(sample,2,2))->sample1
  sc_sample<-reference@meta.data$cell_type_l1[which(reference$orig.ident==sample)]
  st_sample<-merged_filter@meta.data$cell_type_l1_all[which(merged_filter$orig.ident==sample1)]
  table(sc_sample)/sum(table(sc_sample))->sc_prop
  table(st_sample)/sum(table(st_sample))->st_prop
  prop_table$ST<-st_prop[rownames(prop_table)]
  prop_table$SC<-sc_prop[rownames(prop_table)]
  prop_table->prop_table_list[[sample]]
}
lapply(prop_table_list, function(x) {
  print(cor(x)[,2])
})

library(tidyr)
library(dplyr)
all_all->all
all$row<-rownames(all)
all <- pivot_longer(
  all,
  cols = -row,                 # pivot all original columns
  names_to = "Column",             # these are your original column names
  values_to = "Proportion"
)
all$row[which(all$row=='Lymphocyte')]<-'T_cell'
all$row<-factor(all$row,levels = c('Myeloid','T_cell','B_Plasma','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
ggplot(all, aes(y = Column, x = Proportion, fill = row)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = colors) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))->p
ggsave('p.png',p,width = 8,height = 3)

# stacked bars (proportions per column)
all_sep->all
all$row<-rownames(all)
all <- pivot_longer(
  all,
  cols = -row,                 # pivot all original columns
  names_to = "Column",             # these are your original column names
  values_to = "Proportion"
)

all$row[which(all$row=='Lymphocyte')]<-'T_cell'
all$row<-factor(all$row,levels = c('Myeloid','T_cell','B_Plasma','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
colors <- c(
  "Fibroblasts"="#7D0CEA",  # deep teal
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
ggplot(all[grep('LS',all$Column),], aes(y = Column, x = Proportion, fill = row)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = colors) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))->p
ggsave('p.png',p,width = 8,height = 3)


ggplot(all[grep('SSC',all$Column),], aes(y = Column, x = Proportion, fill = row)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = colors) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))->p
ggsave('p1.png',p,width = 8,height = 3)

#pieplot for ST excluding epidermis
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi@meta.data->meta
meta$cell_type_l1_all[which(meta$cell_type_l1_all%in%c('Macrophage','Neutrophil'))]<-'Myeloid'
meta$cell_type_l1_all<-factor(meta$cell_type_l1_all,levels = c('Myeloid','T_cell','B_Plasma','Fibroblasts','Keratinocyte','Glandular_epithelium','Endothelium','Smooth_muscle_cell','Neuron','Melanocyte'))
ggplot(as.data.frame(meta), aes(x = "", fill = cell_type_l1_all)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()->p
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")


#dimplot
PercentageFeatureSet(merged_filter,pattern = 'MT-')->merged_filter$percent.mt
DimPlot(merged_filter[,which(merged_filter$percent.mt<20)],group.by = 'cell_type_l1_all',cols = colors)
DimPlot(merged_filter[,which(merged_filter$percent.mt<20)],group.by = 'cell_type_l1_all')->p
ggsave('p.png',p,width=8,height = 6.5,dpi = 300)

reference$cell_type_l1[which(reference$cell_type_l1=='Lymphocyte')]<-'T_cell'
reference$cell_type_l1[which(reference$seurat_clusters==26)]<-'Neuron'
reference$cell_type_l1[which(reference$cell_type_high=='Myofibroblast')]<-'Smooth_muscle_cell'
DimPlot(reference,group.by = 'cell_type_l1')->p
ggsave('p_sc.png',p,width=8,height = 6.5,dpi = 300)

#plot heatmap
library(Seurat)
Idents(merged_filter)<-merged_filter$cell_type_l1_all
FindAllMarkers(merged_filter,only.pos = T)->marker_all




