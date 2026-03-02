merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds")

#DE gene list (Panel G)
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

merged_filter_exclude_epi$cell_type_l1_all2<-as.character(merged_filter_exclude_epi$cell_type_l1_all)
merged_filter_exclude_epi$cell_type_l1_all2[which(merged_filter_exclude_epi$cell_type_l1_all2%in%c('Neutrophil','Macrophage'))]<-'Myeloid'
Idents(merged_filter_exclude_epi)<-merged_filter_exclude_epi$cell_type_l1_all2
FindAllMarkers(merged_filter_exclude_epi,only.pos = T)->marker_all
marker_all_sig<-marker_all[which(marker_all$p_val_adj<0.05&marker_all$avg_log2FC>1),]
marker_all_sig %>% group_by(cluster) %>%  slice_head(n = 4) %>% ungroup() -> top5
FindMarkers(merged_filter_exclude_epi,ident.1 = 'Endothelium',group.by = 'cell_type_l1_all',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker1
FindMarkers(merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident=='F_1')],ident.1 = 'Endothelium',group.by = 'cell_type_l1_all',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker2
FindMarkers(merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident=='M_1')],ident.1 = 'Endothelium',group.by = 'cell_type_l1_all',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker3
FindMarkers(merged_filter_exclude_epi[,which(merged_filter_exclude_epi$condition=='LS')],ident.1 = 'Endothelium',group.by = 'cell_type_l1_all',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker3
FindMarkers(merged_filter_exclude_epi[,which(grepl('1',merged_filter_exclude_epi$orig.ident))],ident.1 = 'LS',ident.2='SSC',group.by = 'condition',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker4
FindMarkers(merged_filter_exclude_epi[,which(grepl('1',merged_filter_exclude_epi$orig.ident)&merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')],ident.1 = 'LS',ident.2='SSC',group.by = 'condition',logfc.threshold = 0,min.pct = 0,features = c('PECAM1','VCAM1','ICAM1','ICAM2'))->marker5

merged_filter_exclude_epi$endo<-F
merged_filter_exclude_epi$endo[which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')]<-T
dot_plot(merged_filter_exclude_epi,features = c('PECAM1','ICAM1','VCAM1','ICAM2'),'endo',10)->p1
dot_plot(merged_filter_exclude_epi[,which(grepl('1',merged_filter_exclude_epi$orig.ident))],features = c('PECAM1','ICAM1','VCAM1','ICAM2'),'condition',10)->p2
dot_plot(merged_filter_exclude_epi[,which(grepl('1',merged_filter_exclude_epi$orig.ident)&merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')],features = c('PECAM1','ICAM1','VCAM1','ICAM2'),'condition',10)->p3
dot_plot(merged_filter_exclude_epi[,which(merged_filter_exclude_epi$condition=='LS')],features = c('PECAM1','ICAM1','VCAM1','ICAM2'),'endo',10)->p4


DoHeatmap(merged_filter_exclude_epi,features = top5$gene,group.by = 'cell_type_l1_all2')
merged_filter_exclude_epi$cell_type_l1_all2<-factor(merged_filter_exclude_epi$cell_type_l1_all2,levels = names(sort(table(merged_filter_exclude_epi$cell_type_l1_all2), decreasing = TRUE)))
top5<-as.data.frame(do.call(rbind,lapply(levels(merged_filter_exclude_epi$cell_type_l1_all2), function(x) {top5[which(top5$cluster==x),]})))
DotPlot(merged_filter_exclude_epi,features = top5$gene,group.by = 'cell_type_l1_all2')+theme(axis.text.x=element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.text.y = element_text(size=20))->p
ggsave('p1.png',p,width = 15,height = 5,dpi = 300)
p_bar <- ggplot(top5,
                aes(x = factor(gene,levels =top5$gene ),
                    y = 1, fill = cluster)) +
  geom_tile() +
  scale_fill_manual(values = colors)+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(drop = FALSE) +
  # choose your own colors here if you like
  # scale_fill_manual(values = c("Kera"="red","Fibro"="blue","Endo"="green")) +
  theme_void() +
  theme(
    legend.position = "bottom",          # or "none" if you don’t want a legend
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  ) +
  guides(fill = guide_legend(title = "Gene group"))

p / p_bar + plot_layout(heights = c(10, 1))


ggplot(top5, aes(x = factor(gene,levels = top5$gene), y = 0)) +
  geom_text(
    aes(label = gene, color = cluster),
    angle = 60,
    hjust = 1,
    vjust = 1,
    size  = 5,
    fontface = "bold"
  ) +
  scale_color_manual(values = colors) +  # or omit this to use ggplot defaults
  coord_cartesian(clip = "off") +
  theme_void() +  # remove axes, grid, etc.
  theme(
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10),
    legend.position = "right"   # or "none" if you don’t want a legend
  )->p
ggsave('p.png',p,width = 15,height = 3.5,dpi = 300)



top5$color<-colors[top5$cluster]
  # custom colored text
  geom_text(
    data = top5,
    aes(
      x = gene,
      y = -0.4,           # position below zero (adjust as needed)
      label = gene,
      color = cluster
    ),
    angle = 60,
    hjust = 1,
    vjust = 1,
    size = 3
  ) +
  scale_color_manual(values = colors) +
  coord_cartesian(clip = "off") +        # allow text outside plot area
  theme(plot.margin = margin(b = 30))    # add room for labels

ggsave('p.png',p,width=15,heigh=6,dpi=300)
VariableFeatures(merged_filter_exclude_epi)<-union(VariableFeatures(merged_filter_exclude_epi),top5$gene)


#pieplot
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

merged_filter_exclude_epi$cell_type_l1_all<-factor(merged_filter_exclude_epi$cell_type_l1_all,levels = c('Macrophage','T_cell','B_Plasma','Neutrophil','Keratinocyte','Fibroblasts','Endothelium','Glandular_epithelium','Smooth_muscle_cell','Neuron','Melanocyte'))
ggplot(data.frame(cell_type=merged_filter_exclude_epi$cell_type_l1_all), aes(x = "", fill = cell_type)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Pie chart of labels") +
  theme_void()->p
ggsave('p.png',p,width = 5,height = 4,dpi = 300)
