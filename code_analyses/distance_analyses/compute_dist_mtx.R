#distance between each 2 individual cells
merged_filter_exclude_epi<-readRDS('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds')
list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/')->sample
result<-list()
for(i in sample[which(nchar(sample)==3)][-c(2,14)]){
  read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/all_reference/doublet/unfilter_thresh50/result_',i,'.csv'))->result[[i]]
} 
do.call(rbind,result)->result
result_exclude_epi<-result[match(colnames(merged_filter_exclude_epi),result$index),]
result<-as.data.frame(merged_filter@meta.data[,c('x','y')])
dist_matrix_all<-list()
for (i in unique(substr(rownames(result_exclude_epi),1,3))) {
  result_exclude_epi[grep(i,rownames(result_exclude_epi)),]->a
  as.matrix(dist(a))->dist_matrix_all[[i]]
}
lapply(1:13, function(x) {sum(dist_matrix_all2[[x]]==dist_matrix_all[[x]])==dim(dist_matrix_all[[x]])[1]^2})
saveRDS(dist_matrix_all,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/dist_matrix_all.rds')

dist_matrix_all<-readRDS('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/dist_matrix_all.rds')
dist_matrix_all$G_1<-NULL
test_dist<-function(cell_type){
  names<-colnames(merged_filter_exclude_epi)[which(merged_filter_exclude_epi$cell_type_l1_all%in%cell_type)]
  lapply(dist_matrix_all, function(x) {
    which(colnames(x)%in%names)->idx
    if (length(idx)>1){
      as.data.frame(apply(x[,idx],1,min)) 
    }
  })->all
  do.call(rbind,all)->all
  return(all[match(colnames(merged_filter_exclude_epi),gsub('.*[.]','',rownames(all))),1])
}

test_dist('Fibroblasts')->merged_filter_exclude_epi$min_dis_to_fib
test_dist(c('Macrophage','Neutrophil','T_cell','B_Plasma'))->merged_filter_exclude_epi$min_dis_to_immune
test_dist('Fibroblasts')->merged_filter_exclude_epi$min_dis_to_fib
test_dist('Endothelium')->merged_filter_exclude_epi$min_dis_to_endo
test_dist('Keratinocyte')->merged_filter_exclude_epi$min_dis_to_kera
test_dist(c('Macrophage','Neutrophil'))->merged_filter_exclude_epi$min_dis_to_Myeloid
test_dist('Macrophage')->merged_filter_exclude_epi$min_dis_to_Macrophage
test_dist('Glandular_epithelium')->merged_filter_exclude_epi$min_dis_to_gland
test_dist('T_cell')->merged_filter_exclude_epi$min_dis_to_T
test_dist('endo_fib')->merged_filter_exclude_epi$min_dis_to_endo_fib
saveRDS(merged_filter_exclude_epi,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds')
write.csv(as.data.frame(merged_filter_exclude_epi@meta.data[,grep('min_dis_to',colnames(merged_filter_exclude_epi@meta.data))]),'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/distance_all.csv',row.names = T)
