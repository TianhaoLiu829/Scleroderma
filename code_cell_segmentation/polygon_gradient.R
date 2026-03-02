library(Seurat)
list.files('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1')->sample
unique(substr(sample[grep('result',sample)],8,10))->sample
spot_to_polygon_list<-list()
for (i in sample) {
  spot_to_polygon<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spots_gradient/spots_gradient_',i,'_80.csv'))
  count_mtx<-Read10X_h5(paste0('/ix1/wchen/liutianhao/result/pathology_ST/count_mtx/2um/filtered_feature_bc_matrix_',gsub('_','',i),'.h5'))
  colSums(count_mtx)->sum
  intersect(names(sum),spot_to_polygon$barcode)->inter
  spot_to_polygon<-spot_to_polygon[match(inter,spot_to_polygon$barcode),]
  spot_to_polygon$UMI<-sum[inter]
  spot_to_polygon->spot_to_polygon_list[[i]]
}
spot_to_polygon_list_sum<-lapply(spot_to_polygon_list, function(x) aggregate(x$UMI,by=list(x$band),sum))
spot_to_polygon_list_sum->spot_to_polygon_list_sum2

for (i in sample) {
  area<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/script/','global_band_areas_um_',i,'_80.csv'))
  colnames(area)<-gsub('[.]','-',gsub('X','',colnames(area)))
  read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/position_8um/',gsub('_','',i),'_positions.csv') )->position
  area['>=80']<-sum(position$in_tissue)*64
  spot_to_polygon_list_sum[[i]]->x
  x<-x[match(names(area),x$Group.1),]
  x$area<-as.numeric(as.matrix(area))
  x$density<-x$x/x$area
  x->spot_to_polygon_list_sum[[i]]
}

do.call(rbind,lapply(spot_to_polygon_list_sum, function(x) x$x))->all_umi
do.call(rbind,lapply(spot_to_polygon_list_sum, function(x) x$area))->all_area
a<-colSums(all_umi[which(rownames(all_umi)!='G_1'),])/colSums(all_area[which(rownames(all_area)!='G_1'),])
plot(a)
df<-data.frame(inter=spot_to_polygon_list_sum$A_2$Group.1,density=a)
df$inter[grep('-',df$inter)]<-unlist(lapply(strsplit(df$inter[grep('-',df$inter)],'-'), function(x) paste(as.character(as.numeric(x)/4),collapse = '-')))
df$color<-'intermediate'
df$color[c(1,dim(df)[1])]<-c('inside','outside')
ggplot(df,aes(x=factor(inter,levels=df$inter),y=density,colour = color)) +  
  geom_point() +
  scale_color_manual(values = c('inside'='red','intermediate'='black','outside'='blue')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,size = 10,hjust = 1),axis.text.y=element_text(size=13))->p
ggsave('p.png',p,width = 7.5,height = 4,dpi = 300)
diff(a)

pheatmap::pheatmap(t(df$density),cluster_rows = F,cluster_cols = F,border_color = NA)->p
ggsave('p1.png',p,width = 6.5,height = 3,dpi = 300)


#check the quality
spot_to_polygon_list2<-lapply(names(spot_to_polygon_list), function(x) {
  spot_to_polygon_list[[x]]->mtx
  mtx$density<-spot_to_polygon_list_sum[[x]]$density[match(mtx$band,spot_to_polygon_list_sum[[x]]$Group.1)]
  return(mtx)
})
names(spot_to_polygon_list2)<-names(spot_to_polygon_list)

for (i in sample) {
  gsub('>=80','84',gsub('inside','0',gsub('.*-','',spot_to_polygon_list2[[i]]$band)))->spot_to_polygon_list2[[i]]$band
  write.csv(spot_to_polygon_list2[[i]],paste0(i,'.csv'))
}






