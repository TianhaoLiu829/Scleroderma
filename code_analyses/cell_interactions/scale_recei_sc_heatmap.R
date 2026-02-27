#for the individual signal comparison (receiver_sc)
recei_list<-lapply(result, function(x){
  lapply(x, function(y){
    return(y[['recei_sc']])
  })
})

L_switched <- lapply(seq_along(recei_list[[1]]), function(i) {
  lapply(recei_list, function(x) x[[i]])
})
names(L_switched)<-names(recei_list$A_2)
recei_list<-lapply(L_switched,function(x) {
  return(x[which(names(x)!='G_1')])
})

Reduce(intersect,lapply(recei_list[["commot-cellchat-AGRN-DAG1"]],rownames))->common_cells
recei_list_filter<-lapply(recei_list, function(x) {
  lapply(x, function(y){y[common_cells,]})
})
recei_list_filter<-lapply(recei_list_filter, function(x) {
  as.matrix(do.call(cbind,x))->df
  df<-df[,colnames(merged_filter_base_exclude_epi_LS_SSC)]
  df
})

lapply(recei_list_filter,colnames)->name
all(sapply(name[-1], function(x) identical(x, name[[1]])))
sum(name$`commot-cellchat-AGRN-DAG1`==colnames(merged_filter_base_exclude_epi_LS_SSC))==dim(merged_filter_base_exclude_epi_LS_SSC)[2]

recei_cells<-rownames(recei_list_filter$`commot-cellchat-AGRN-DAG1`)
recei_list_filter<-lapply(recei_cells,function(x) {
  lapply(recei_list_filter, function(y) y[x,])->list
  as.data.frame(do.call(cbind,list))
})
names(recei_list_filter)<-recei_cells
all_list<-list()
for (recei in common_cells) {
  merged_filter_base_exclude_epi_LS_SSC@meta.data[which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all%in%recei),]->target
  which(target$condition=='SSC')->SSC
  which(target$condition=='LS')->LS
  send_list<-list()
  for (sender in common_cells) {
    recei_list_filter[[sender]]->recei_list_filter_endo_recei_LS_SSC
    
    pvals <- apply(recei_list_filter_endo_recei_LS_SSC[rownames(target),], 2, function(x) {
      wilcox.test(x[SSC], x[LS])->test
      mean(x[SSC])->ssc
      mean(x[LS])->ls
      length(which(x[LS]>0))/length(x[LS])->ls_pct
      length(which(x[SSC]>0))/length(x[SSC])->ssc_pct
      c(stat = test$statistic,p=test$p.value,ssc=ssc,ls=ls,ssc_pct=ssc_pct,ls_pct=ls_pct,diff=ssc-ls,log2fc=log2(ssc/ls) )
    })
    pvals<-as.data.frame(t(pvals))
    pvals<-pvals[which(abs(pvals$diff)>0),]
    pvals<-pvals[order(pvals$p),]
    pvals->send_list[[sender]]
  }
  all_list[[recei]]<-send_list
}
saveRDS(all_list,'all_list_inf_recei_sender.rds')
saveRDS(all_list,'all_list_recei_sender.rds')
all_list_inf_recei_sender<-readRDS('all_list_inf_recei_sender.rds')

all_list<-list()
for (recei in common_cells) {
  merged_filter_base_exclude_epi_LS_SSC@meta.data[which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all%in%recei),]->target
  which(target$condition=='SSC')->SSC
  which(target$condition=='LS')->LS
  send_list<-list()
  for (sender in common_cells) {
    recei_list_filter[[sender]]->recei_list_filter_endo_recei_LS_SSC
    recei_list_filter_endo_recei_LS_SSC[rownames(target),all_lr]->mtx
    mtx_scale<-scale(mtx)
    pvals <- lapply(1:ncol(mtx), function(x1) {
      x<-mtx[,x1]
      wilcox.test(x[SSC], x[LS])->test
      mean(x[SSC])->ssc
      mean(x[LS])->ls
      x_scale<-mtx_scale[,x1]
      mean(x_scale[SSC])->ssc_scale
      mean(x_scale[LS])->ls_scale
      length(which(x[LS]>0))/length(x[LS])->ls_pct
      length(which(x[SSC]>0))/length(x[SSC])->ssc_pct
      c(stat = test$statistic,p=test$p.value,ssc=ssc,ssc_scale=ssc_scale,ls=ls,ls_scale=ls_scale,ssc_pct=ssc_pct,ls_pct=ls_pct,diff=ssc-ls,log2fc=log2(ssc/ls) )
    })
    pvals<-as.data.frame(do.call(rbind,pvals))
    rownames(pvals)<-all_lr
    pvals<-pvals[which(abs(pvals$diff)>0),]
    pvals<-pvals[order(pvals$p),]
    pvals->send_list[[sender]]
  }
  all_list[[recei]]<-send_list
}
