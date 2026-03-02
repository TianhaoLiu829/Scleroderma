load("/ix1/wchen/liutianhao/work_10_05.RData")
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
skin_ST_combine_harmony <- readRDS("~/ST/ST_work/skin/skin_ST_combine_harmony.rds")
colors_dark11 <- c("blue", "#D95F02",  "#7570B3",  "#E7298A",  "#66A61E",  "#E6AB02","#A6761D",  "#666666","#1F78B4","#B15928",  "#6A3D9A" )

#peri_vascular
pvals<-data.frame()
for (thresh in seq(120,160,5)) {
  thresh<-132
  merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
  merged_filter_exclude_epi_endo<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Endothelium')]
  merged_filter_exclude_epi_endo$inf_endo<-FALSE
  #merged_filter_exclude_epi_endo$inf_endo[which(colnames(merged_filter_exclude_epi_endo)%in%inf_endo_index)]<-T
  merged_filter_exclude_epi_endo$inf_endo[which(apply(merged_filter_exclude_epi_endo@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T','min_dis_to_B_Plasma')], 1, function(x) {sum(x<thresh)>0}))]<-T
  as.data.frame(prop.table(table(merged_filter_exclude_epi_endo$inf_endo,merged_filter_exclude_epi_endo$orig.ident),margin = 2))->prop
  prop<-prop[which(prop$Var1==TRUE),]
  t.test(prop$Freq[which(prop$Var2%in%c('E_1','H_1','K_1'))],prop$Freq[which(prop$Var2%in%c('I_1','J_1','F_1','M_1'),)])->test
  test$p.value->pvals[as.character(thresh),'p']
  prop$condition<-merged_filter_exclude_epi$condition[match(prop$Var2,merged_filter_exclude_epi$orig.ident)]
  prop[grepl('1',prop$Var2)&prop$Var2!='G_1',]->prop
  ggplot(prop, aes(x = condition, y = Freq, fill = condition)) +
    geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
    geom_jitter(width = 0.1, size = 1.5, color = "black") + # raw points
    theme_classic() +
    labs(title = "Frequency of inflammed Endothelium")
  
}
skin_ST_combine_harmony@images$sliceB_1@coordinates->l
skin_ST_combine_harmony<-skin_ST_combine_harmony[,which(!colnames(skin_ST_combine_harmony)%in%rownames(l)[which(l$imagecol>7000)])]
apply(skin_ST_combine_harmony@assays$conv@counts,2,function(x) {rownames(skin_ST_combine_harmony@assays$conv@counts)[which.max(x)]})->k
#endo_visium<-skin_ST_combine_harmony[,which(skin_ST_combine_harmony$first_type=='Endothelial'|skin_ST_combine_harmony$second_type=='Endothelial')]
#endo_visium$annotation4<-'other'
#endo_visium$annotation4[which(colSums(as.matrix(endo_visium@assays$conv@counts[c('T cell','Macrophages','plasma cell','B cells'),]))>0.14)]<-'peri_vascular'
endo_visium<-skin_ST_combine_harmony[,which(skin_ST_combine_harmony$annotation4%in%c('Blood vessel','Immune_perivascular','Immune_perivascular/perineural')|skin_ST_combine_harmony$first_type=='Endothelial')]
endo_visium$annotation4[which(endo_visium$annotation4%in%c('Immune_perivascular','Immune_perivascular/perineural'))]<-'peri_vascular'
as.data.frame(table(endo_visium$annotation4,endo_visium$orig.ident))->a
table(endo_visium$orig.ident)[a$Var2]->a$total
a[grep('peri_vascular',a$Var1),]->a
a$Freq<-a$Freq/a$total
a<-a[,-4]
a->s

#as.data.frame(prop.table(table(endo_visium@meta.data$annotation4,endo_visium@meta.data$orig.ident),margin = 2))->s
#s[grep('peri_vascular',s$Var1),]->s
s$Var2<-c('A_1','B_1','C_1','D_1')
s$Var1<-T
s$condition<-c('SSC',rep('LS',3))
as.data.frame(rbind(prop,s))->prop2
t.test(prop2$Freq[which(prop2$Var2%in%c('A_1','E_1','H_1','K_1'))],prop2$Freq[which(prop2$Var2%in%c('B_1','C_1','D_1','I_1','J_1','F_1','M_1'),)])
prop2$Var2<-factor(prop2$Var2,levels = c('A_1','E_1','H_1','K_1','B_1','C_1','D_1','I_1','J_1','F_1','M_1'))

set.seed(50)
ggplot(prop2, aes(x = condition, y = Freq*100, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
  scale_fill_manual(values = c("LS" = "#E9967A", "SSC" = "#66C2A5")) +
  geom_jitter(width = 0.1, size = 2, aes(color=Var2)) + # raw points
  scale_color_manual(values = colors_dark11) +
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  labs(title = "Frequency of inflammed Endothelium")->p
set.seed(50)
ggsave('p.png',width = 5.5,height = 5,dpi = 300)


PGAD<-c(B_1=36,C_1=37,D_1=28,F_1=81,I_1=8,J_1=35,M_1=80)
PGAA<-c(B_1=33,C_1=43,D_1=82,F_1=78,I_1=25,J_1=76,M_1=0)
prop2$PGAD<-PGAD[as.character(prop2$Var2)]
prop2$PGAA<-PGAA[as.character(prop2$Var2)]
prop2_LS<-prop2[which(prop2$condition=='LS'),]

ggplot(prop2_LS, aes(x = Freq, y = PGAA)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave('p_pgaa.png',p,width=5,height = 5,dpi = 300)
ggplot(prop2_LS[which(prop2_LS$Var2!='M_1'),], aes(x = Freq, y = PGAD)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave('p_pgad.png',p,width=5,height = 5,dpi = 300)

cor.test(log2(prop2_LS$Freq+1),log2(prop2_LS$PGAA+1))
cor.test(prop2_LS$Freq,prop2_LS$PGAA)
cor.test(log2(prop2_LS$Freq+1),log2(prop2_LS$PGAA+1))
cor.test(prop2_LS$Freq,prop2_LS$PGAA)

prop2_LS<-prop2_LS[which(prop2_LS$Var2!='M_1'),]
cor.test(log2(prop2_LS$Freq+1),log2(prop2_LS$PGAD+1))
cor.test(prop2_LS$Freq,prop2_LS$PGAD)
cor.test(log2(prop2_LS$Freq+1),log2(prop2_LS$PGAD+1))
cor.test(prop2_LS$Freq,prop2_LS$PGAD)


lossi<-c(A_1=28,B_1=1,C_1=4,D_1=0,E_1=23,F_1=52,H_1=36,I_1=5,J_1=3,K_1=24,M_1=0)
lossd<-c(A_1=28,B_1=9,C_1=8,D_1=8,E_1=23,F_1=68,H_1=36,I_1=8,J_1=14,K_1=24,M_1=4)

prop2$lossi<-lossi[as.character(prop2$Var2)]
prop2$lossd<-lossd[as.character(prop2$Var2)]

cor.test(prop2$lossi,prop2$Freq,method = 'spearman')
cor.test(prop2$lossd,prop2$Freq,method = 'spearman')

cor.test(log2(prop2$lossi+1),prop2$Freq,method = 'pearson')
cor.test(log2(prop2$lossd+1),prop2$Freq,method = 'pearson')


#peri_fibroblast
pvals<-data.frame()
for (thresh in seq(20,400,20)) {
  thresh<-132
  merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
  merged_filter_exclude_epi_fib<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$cell_type_l1_all=='Fibroblasts')]
  merged_filter_exclude_epi_fib$inf_fib<-FALSE
  merged_filter_exclude_epi_fib$inf_fib[which(apply(merged_filter_exclude_epi_fib@meta.data[,c('min_dis_to_Macrophage','min_dis_to_T')], 1, function(x) {sum(x<thresh)>0}))]<-T
  as.data.frame(prop.table(table(merged_filter_exclude_epi_fib$inf_fib,merged_filter_exclude_epi_fib$orig.ident),margin = 2))->prop
  prop<-prop[which(prop$Var1==TRUE),]
  t.test(prop$Freq[which(prop$Var2%in%c('E_1','H_1','K_1'))],prop$Freq[which(prop$Var2%in%c('I_1','J_1','F_1','M_1'),)])->test
  test$p.value->pvals[as.character(thresh),'p']
  prop$condition<-merged_filter_exclude_epi$condition[match(prop$Var2,merged_filter_exclude_epi$orig.ident)]
  prop[grepl('1',prop$Var2)&prop$Var2!='G_1',]->prop
  ggplot(prop, aes(x = condition, y = Freq, fill = condition)) +
    geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
    geom_jitter(width = 0.1, size = 1.5, color = "black") + # raw points
    theme_classic() +
    scale_color_manual(values = colors_dark11) +
    labs(title = "Frequency of inflammed Endothelium")
  
}
pvals<-data.frame()
t<-1
for (thresh1 in seq(0,0.3,0.05)) {
  for (thresh2 in seq(0,0.3,0.05)) {
    #skin_sc <- readRDS("~/ST/ST_work/skin/result/skin_sc.rds")
    skin_ST_combine_harmony@images$sliceB_1@coordinates->l
    skin_ST_combine_harmony<-skin_ST_combine_harmony[,which(!colnames(skin_ST_combine_harmony)%in%rownames(l)[which(l$imagecol>7000)])]
    #as.data.frame(aggregate(skin_sc$nCount_RNA,by=list(skin_sc$cell_type),mean))->a
    #a$Group.1[12]<-'Myelin cells'
    #t(as.matrix(skin_ST_combine_harmony@assays$conv@counts))%*%a$x->mean
    #apply(skin_ST_combine_harmony@assays$conv@counts,2,function(x) {rownames(skin_ST_combine_harmony@assays$conv@counts)[which.max(x)]})->k
    #fib_visium<-skin_ST_combine_harmony[,which(skin_ST_combine_harmony@assays$conv@counts['Fibroblasts',]>thresh1)]
    fib_visium<-skin_ST_combine_harmony[,which(skin_ST_combine_harmony$first_type=='Fibroblasts'|skin_ST_combine_harmony$second_type=='Fibroblasts')]
    #fib_visium<-fib_visium[,which(fib_visium@assays$conv@counts['Fibroblasts',])]
    #fib_visium$annotation4[which(colSums(as.matrix(fib_visium@assays$conv@counts[c(1,9,10,14,15),]))>thresh2)]<-'peri_fib'
    fib_visium$annotation4[fib_visium$first_type%in%c('T cell','Macrophages','plasma cell','B cells')|fib_visium$second_type%in%c('T cell','Macrophages','plasma cell','B cells')]<-'peri_fib'
    
    as.data.frame(table(fib_visium$annotation4,fib_visium$orig.ident))->a
    table(fib_visium$orig.ident)[a$Var2]->a$total
    a[grep('peri_fib',a$Var1),]->a
    a$Freq<-a$Freq/a$total
    a<-a[,-4]
    a->s
    #as.data.frame(prop.table(table(endo_visium@meta.data$annotation4,endo_visium@meta.data$orig.ident),margin = 2))->s
    #s[grep('peri_vascular',s$Var1),]->s
    s$Var2<-c('A_1','B_1','C_1','D_1')
    s$Var1<-T
    s$condition<-c('SSC',rep('LS',3))
    as.data.frame(rbind(prop,s))->prop2
    t.test(prop2$Freq[which(prop2$Var2%in%c('A_1','E_1','H_1','K_1'))],prop2$Freq[which(prop2$Var2%in%c('B_1','C_1','D_1','I_1','J_1','F_1','M_1'),)])->test
    test$p.value->pvals[t,'p']
    thresh1->pvals[t,'thresh1']
    thresh2->pvals[t,'thresh2']
    prop2$Var2<-factor(prop2$Var2,levels = c('A_1','E_1','H_1','K_1','B_1','C_1','D_1','I_1','J_1','F_1','M_1'))
    set.seed(500)
    ggplot(prop2, aes(x = condition, y = Freq*100, fill = condition)) +
      geom_bar(stat = "summary", fun = mean, alpha = 0.6, width = 0.6) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # error bars
      geom_jitter(width = 0.1, size = 2, aes(color=Var2),alpha = 1.5) + # raw points
      scale_color_manual(values = colors_dark11) +
      theme_classic() +
      scale_fill_manual(values = c("LS" = "#E9967A", "SSC" = "#66C2A5")) +
      theme(axis.text = element_text(size = 15)) +
      labs(title = "Frequency of inflammed Endothelium")->p
    set.seed(500)
    ggsave("p1.png",p,width = 5.5,height = 5)
    t<-t+1
  }
}
pvals[which(pvals$thresh1==pvals$thresh2),]

PGAD<-c(B_1=36,C_1=37,D_1=28,F_1=81,I_1=8,J_1=35,M_1=80)
PGAA<-c(B_1=33,C_1=43,D_1=82,F_1=78,I_1=25,J_1=76,M_1=0)
prop2$PGAD<-PGAD[as.character(prop2$Var2)]
prop2$PGAA<-PGAA[as.character(prop2$Var2)]
prop2_LS<-prop2[which(prop2$condition=='LS'),]

ggplot(prop2_LS, aes(x = Freq, y = PGAA)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave('p_pgaa.png',p,width=5,height = 5,dpi = 300)
ggplot(prop2_LS[which(prop2_LS$Var2!='M_1'),], aes(x = Freq, y = PGAD)) +
  geom_point(size = 2) +                 # scatter points
  geom_smooth(method = "lm", se = F) +# fitted linear regression line
  theme_classic() +
  theme(axis.text = element_text(size = 15)) +
  labs(x = "X value", y = "Y value")->p
ggsave('p_pgad.png',p,width=5,height = 5,dpi = 300)

cor.test(prop2_LS$Freq,prop2_LS$PGAA,method = 'spearman')
prop2_LS<-prop2_LS[which(prop2_LS$Var2!='M_1'),]
cor.test(prop2_LS$Freq,prop2_LS$PGAD,method = 'spearman')



#test the proportion of immune cells in neighborhood of each individual fib/endo
test_dist<-function(cell_type,thresh){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  as.data.frame(merged_filter_exclude_epi@meta.data[!colnames(merged_filter_exclude_epi)%in%names,'cell_type_l1_all'])->meta
  rownames(meta)<-colnames(merged_filter_exclude_epi)[!colnames(merged_filter_exclude_epi)%in%names]
  colnames(meta)<-'cell_type'
  all_labels<-unique(merged_filter_exclude_epi$cell_type_l1_all)
  lapply(dist_matrix_all, function(x) {
    x[!rownames(x)%in%names,colnames(x)%in%names]->filtered_x
    apply(filtered_x,2,function(y) {
      meta[rownames(filtered_x)[which(y<thresh)],'cell_type']->labels
      table(factor(labels, levels = all_labels))->counts
      as.numeric(counts) / length(labels)
    })
  })->all
  
  do.call(cbind,all)->all
  all[which(is.na(all))]<-0
  as.data.frame(t(all))->all
  colnames(all)<-all_labels
  all<-all[rowSums(all)>0,]
  all$sample<-substr(rownames(all),1,3)
  all$immune<-rowSums(all[,c('Macrophage','T_cell','B_Plasma')])
  all$condition<-merged_filter_exclude_epi$condition[match(all$sample,merged_filter_exclude_epi$orig.ident)]
  all$treatment<-substr(all$sample,3,3)
  return(all)
}
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
test_dist('Fibroblasts',132)->inf_fib_prop
test_dist('Endothelium',132)->inf_endo_prop
test_dist('Endothelium',220)->inf_endo_prop3
merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_endo)[which(merged_filter_exclude_epi_endo$inf_endo)])]<-'inf_endo'
merged_filter_exclude_epi$cell_type_l1_all[which(colnames(merged_filter_exclude_epi)%in%colnames(merged_filter_exclude_epi_fib)[which(merged_filter_exclude_epi_fib$inf_fib)])]<-'inf_fib'
test_dist('inf_fib',132)->inf_fib_prop2
test_dist('inf_endo',132)->inf_endo_prop2


wilcox.test(immune~condition,data=inf_fib_prop[grepl('_1',inf_fib_prop$sample),])->test1
wilcox.test(immune~condition,data=inf_endo_prop[grepl('_1',inf_endo_prop$sample),])->test2
ggplot(inf_fib_prop[grepl('_1',inf_fib_prop$sample),], aes(x = condition, y =immune*100, fill =condition)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.05, size = 0.3, alpha = 0.55) +
  geom_violin(trim = T, width = 1.3, color = "black", alpha = 0.85) +
  scale_fill_manual(values = c("LS" = "#E9C2B3", "SSC" = "#9ED5C6")) +  # pick colors you like
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18)  # big LS / SSc labels
  )->p
ggsave('p.png',p,width = 5.3,height = 6,dpi = 300)

ggplot(inf_endo_prop[grepl('_1',inf_endo_prop$sample),], aes(x = condition, y =immune*100, fill =condition)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.05, size = 0.3, alpha = 0.55) +
  geom_violin(trim = T, width = 1.3, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("LS" = "#E9C2B3", "SSC" = "#9ED5C6")) +  # pick colors you like
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18)  # big LS / SSc labels
  )->p
ggsave('p1.png',p,width = 5.3,height = 6,dpi = 300)


wilcox.test(immune~sample,data=inf_endo_prop2[inf_endo_prop2$sample%in%c('E_1','E_2'),])
wilcox.test(immune~sample,data=inf_endo_prop2[inf_endo_prop2$sample%in%c('H_1','H_2'),])
wilcox.test(immune~sample,data=inf_endo_prop2[inf_endo_prop2$sample%in%c('K_1','K_2'),])

#wilcox.test(immune~sample,data=inf_fib_prop2[inf_fib_prop2$sample%in%c('E_1','E_2'),])
#wilcox.test(immune~sample,data=inf_fib_prop2[inf_fib_prop2$sample%in%c('H_1','H_2'),])
#wilcox.test(immune~sample,data=inf_fib_prop2[inf_fib_prop2$sample%in%c('K_1','K_2'),])


wilcox.test(immune~sample,data=inf_endo_prop[inf_endo_prop$sample%in%c('E_1','E_2'),])
wilcox.test(immune~sample,data=inf_endo_prop[inf_endo_prop$sample%in%c('H_1','H_2'),])
wilcox.test(immune~sample,data=inf_endo_prop[inf_endo_prop$sample%in%c('K_1','K_2'),])


ggplot(inf_endo_prop2[inf_endo_prop2$sample%in%c('E_1','E_2'),], aes(x = treatment, y =immune*100, fill =treatment)) +
  geom_violin(trim = T, width = 0.7, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.05, size = 0.4, alpha = 0.6) +
  scale_fill_manual(values = c("1" = "#E9C2B3", "2" = "#9ED5C6")) +  # pick colors you like
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 18)  # big LS / SSc labels
  )->p
ggsave('p_E.png',p,width = 4,height = 4,dpi = 300)

array(dim = length(sample))->area
names(area)<-sample
for (sample in names(area)) {
  read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/position_8um/',gsub('_','',sample),'_positions.csv') )->position
  sum(position$in_tissue)*64/10^6->area[sample]
}
lapply(names(area), function(x) {length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all=='inf_endo'))/area[x]})
lapply(names(area), function(x) {length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all=='inf_endo'))/length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all%in%c('inf_endo','Endothelium')))})

lapply(names(area), function(x) {length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all=='inf_fib'))/area[x]})
lapply(names(area), function(x) {length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all=='inf_fib'))/length(which(merged_filter_exclude_epi$orig.ident==x&merged_filter_exclude_epi$cell_type_l1_all%in%c('inf_fib','Endothelium')))})


endo_visium_A<-endo_visium[,which(endo_visium$orig.ident=='a_SSc')]
#endo_visium_A<-fib_visium[,which(fib_visium$orig.ident=='a_SSc')]
A_1<-colSums(as.matrix(endo_visium_A@assays$conv@counts[c('T cell','Macrophages','plasma cell','B cells'),]))
A_2<-inf_endo_prop3$immune[which(inf_endo_prop3$sample=='A_2')]
wilcox.test(A_1,A_2)
length(A_2)
ggplot(data.frame(immune=c(A_1,A_2),treatment=c(rep("1",length(A_1)),rep("2",length(A_2)))), aes(x = treatment, y =immune*100, fill =treatment)) +
  geom_violin(trim = T, width = 1, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.05, size = 0.4, alpha = 0.6) +
  scale_fill_manual(values = c("1" = "#E9C2B3", "2" = "#9ED5C6")) +  # pick colors you like
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 18)  # big LS / SSc labels
  )->p
ggsave('p_A.png',p,width = 4,height = 4,dpi = 300 )
