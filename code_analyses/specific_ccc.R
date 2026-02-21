#read the merged data
#merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi <- readRDS("~/merged_filter_exclude_epi_endo_fib.rds")
merged_filter_base_exclude_epi_LS_SSC<-merged_filter_exclude_epi[,which(merged_filter_exclude_epi$orig.ident%in%c('E_1','H_1','K_1','F_1','G_1','I_1','J_1','M_1'))]
merged_filter_base_exclude_epi_LS_SSC<-merged_filter_base_exclude_epi_LS_SSC[,which(merged_filter_base_exclude_epi_LS_SSC$orig.ident!='G_1')]


#for the individual signal comparison (sender_sc)
recei_list<-lapply(result, function(x){
  lapply(x, function(y){
    return(t(y[['recei_sc']]))
  })
})

L_switched <- lapply(seq_along(recei_list[[1]]), function(i) {
  lapply(recei_list, function(x) x[[i]])
})
names(L_switched)<-names(recei_list$A_2)
recei_list<-lapply(L_switched,function(x) {
  return(x[which(names(x)!='G_1')])
})

Reduce(intersect,lapply(recei_list[["commot-cellchat-AGRN-DAG1"]],colnames))->common_cells
recei_list_filter<-lapply(recei_list, function(x) {
  lapply(x, function(y){y[,common_cells]})
})
recei_list_filter<-lapply(recei_list_filter, function(x) {
  as.matrix(do.call(rbind,x))->df
  df<-df[colnames(merged_filter_base_exclude_epi_LS_SSC),]
  df
})

lapply(recei_list_filter,rownames)->name
all(sapply(name[-1], function(x) identical(x, name[[1]])))
sum(name$`commot-cellchat-AGRN-DAG1`==colnames(merged_filter_base_exclude_epi_LS_SSC))==dim(merged_filter_base_exclude_epi_LS_SSC)[2]

recei_cells<-colnames(recei_list_filter$`commot-cellchat-AGRN-DAG1`)
recei_list_filter<-lapply(recei_cells,function(x) {
  lapply(recei_list_filter, function(y) y[,x])->list
  as.data.frame(do.call(cbind,list))
})
names(recei_list_filter)<-recei_cells
all_list<-list()
for (sender in common_cells) {
  merged_filter_base_exclude_epi_LS_SSC@meta.data[which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all%in%sender),]->target
  which(target$condition=='SSC')->SSC
  which(target$condition=='LS')->LS
  recei_list<-list()
  for (recei in common_cells) {
    recei_list_filter[[recei]]->recei_list_filter_endo_recei_LS_SSC
    pvals <- apply(recei_list_filter_endo_recei_LS_SSC[rownames(target),], 2, function(x) {
      wilcox.test(x[SSC], x[LS])->test
      ssc=mean(x[SSC])
      ls=mean(x[LS])
      length(which(x[LS]>0))/length(x[LS])->ls_pct
      length(which(x[SSC]>0))/length(x[SSC])->ssc_pct
      c(stat = test$statistic,p=test$p.value,ssc=ssc,ls=ls,ssc_pct=ssc_pct,ls_pct=ls_pct,diff=ssc-ls,log2fc=log2(ssc/ls) )
    })
    pvals<-as.data.frame(t(pvals))
    pvals<-pvals[which(abs(pvals$diff)>0),]
    pvals<-pvals[order(pvals$p),]
    pvals->recei_list[[recei]]
  }
  all_list[[sender]]<-recei_list
}
saveRDS(all_list,'all_list_inf_sender_recei.rds')
saveRDS(all_list,'all_list_sender_recei.rds')
all_list_inf_sender_recei<-readRDS('all_list_inf_sender_recei.rds')



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





aggregate_block_counts <- function(M, lab) {
  stopifnot(nrow(M) == ncol(M))
  
  # If lab isn't named, assume it's in the same order as colnames/rownames
  if (is.null(names(lab))) {
    if (!is.null(colnames(M))) names(lab) <- colnames(M)
    else if (!is.null(rownames(M))) names(lab) <- rownames(M)
    else names(lab) <- as.character(seq_len(ncol(M)))
  }
  
  # Long form without extra packages
  df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  # Map row/col to labels via names
  df$g_row <- lab[as.character(df$Var1)]
  df$g_col <- lab[as.character(df$Var2)]
  
  # Count non-NA entries per (g_row, g_col)
  counts <- xtabs(!is.na(df$Freq) ~ g_row + g_col, data = df)
  
  total_non_na <- sum(!is.na(M))
  props <- counts / total_non_na
  
  list(counts = counts, proportion = props, total_non_na = total_non_na)
}

#for the cell level comparison
result_cell_level<-lapply(result, function(x){
  lapply(x, function(y) {
    y[[2]]
  })
})
result_cell_level<-result_cell_level[-5]
result_cell_level<-lapply(result_cell_level, function(x) {
  x[[grep('total',names(x))]]
})


library(Matrix)
library(Matrix.utils)
library(dplyr)


#weight<-lapply(result, function(x) {
 # merged_filter_exclude_epi$cell_type_l1_all[match(colnames(x[[1]][[5]]),colnames(merged_filter_exclude_epi))]->lab
  #aggregate_block_counts(as.matrix(x[[1]][[5]]),lab)[['proportion']]
#})
#weight<-weight[-5]
#lapply(weight, function(x) {x[common_cells,common_cells]})->weight
common_cells<-rownames(result_cell_level[["I_1"]])
result_cell_level_filter<-lapply(result_cell_level, function(x) {x[common_cells,common_cells]})
#result_cell_level_filter<-lapply(1:length(result_cell_level_filter), function(x) {
  #as.matrix(result_cell_level_filter[[x]])%*%weight[[x]]
#})
#names(result_cell_level_filter)<-names(result)[-5]

#result_cell_level<-lapply(result_cell_level, function(x){Reduce("+", x)})

as.data.frame(Reduce("+", result_cell_level_filter[c(2,6,10)]))/3 ->result_cell_level_filter_SSC
as.data.frame(Reduce("+", result_cell_level_filter[c(4,8,9,12)]))/4->result_cell_level_filter_LS

c(3551,338,286,9573)/sum(c(3551,338,286,9573))->weights
as.data.frame(Reduce("+", Map(function(m, w) m * w, result_cell_level_filter[c(4,8,9,12)], weights)))->result_cell_level_filter_LS

c(2630,1341,1451)/sum(c(2630,1341,1451))->weights
as.data.frame(Reduce("+", Map(function(m, w) m * w, result_cell_level_filter[c(2,6,10)], weights)))->result_cell_level_filter_SSC

result_cell_level_filter_SSC[c('B_Plasma','Macrophage','T_cell','Endothelium','Fibroblasts'),c('B_Plasma','Macrophage','T_cell','Endothelium','Fibroblasts')]->plot_SSC
diag(plot_SSC)<-NA
pheatmap::pheatmap(log10(plot_SSC+1),cluster_rows = F,cluster_cols = F,scale = 'none')
pheatmap::pheatmap(log10(t(result_cell_level_filter_SSC[c('B_Plasma','Macrophage','T_cell'),c('Endothelium','Fibroblasts')]))+3,cluster_rows = F,breaks = breaks,cluster_cols = F,scale = 'none')

result_cell_level_filter_LS[c('B_Plasma','Macrophage','T_cell','Endothelium','Fibroblasts'),c('B_Plasma','Macrophage','T_cell','Endothelium','Fibroblasts')]->plot_LS
diag(plot_LS)<-NA
pheatmap::pheatmap(plot_LS,cluster_rows = F,cluster_cols = F,scale = 'none')
pheatmap::pheatmap(log10(result_cell_level_filter_LS[c('Endothelium','Fibroblasts'),c('B_Plasma','Macrophage','T_cell')])+3,cluster_rows = F,cluster_cols = F,breaks = breaks,scale = 'none')->p

breaks <- seq(0-0.7, 0.5, length.out = 100)


list_all_p<-list()
list_all_effect<-list()
unique(substr(list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/'),1,3) )->sample
list.files('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/',pattern = '.*p.csv')->all_file

for (id in sample) {
  list_sample_p<-list()
  list_sample_effect<-list()
  a<-all_file[grep(id,all_file)]
  for (LR in a) {
    read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/',gsub('_p','_effect',LR)))->effect
    rownames(effect)<-effect[,1]
    effect<-effect[,-1]
    effect->list_sample_effect[[gsub('_p.csv','',gsub('._[1-2]_','',LR))]]
    read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/',LR))->p
    rownames(p)<-p[,1]
    p<-p[,-1]
    p->list_sample_p[[gsub('_p.csv','',gsub('._[1-2]_','',LR))]]
  }
  list_sample_p->list_all_p[[id]]
  list_sample_effect->list_all_effect[[id]]
}

lapply(list_all_effect, function(x) {
  lapply(x, function(y){
    y['Endothelium','Myeloid']
  })
})->s
lapply(s, function(x) {
  unlist(x)
})->s1
as.data.frame(do.call(cbind,s1))->s1
out<-apply(s1, 1, function(x) {
  x[c('E_1','H_1','K_1')]->SSC
  x[c('F_1','I_1','J_1','M_1')]->LS
  wilcox.test(SSC,LS)->test
  return(c(p=test$p.value,diff=mean(SSC)-mean(LS)))
})
out<-as.data.frame(t(out))

out<-apply(s1, 1, function(x) {
  x[c('E_1','H_1','K_1')]->SSC
  x[c('E_2','H_2','K_2')]->LS
  wilcox.test(SSC,LS,paired = TRUE)->test
  return(c(p=test$p.value,diff=mean(SSC)-mean(LS),diff1=SSC[1]-LS[1],diff2=SSC[2]-LS[2],diff3=SSC[3]-LS[3]))
})
out<-as.data.frame(t(out))
out[which(sign(out$diff1.E_1)==sign(out$diff2.H_1)&sign(out$diff1.E_1)==sign(out$diff3.K_1)),]->out

target_cell<-c('Fibroblasts','Endothelium','Lymphocyte','Myeloid','B_Plasma')
list_all_effect$M_1[target_cell,target_cell]->plot
diag(plot)<-NA
pheatmap::pheatmap(plot, cluster_rows = F,cluster_cols = F)


a<-lapply(list_all_effect, function(x){
  as.data.frame(as.table(as.matrix(x)))
})
do.call(cbind,lapply(a,function(x) x$Freq) )->all
rownames(all)<-paste0(a$A_2$Var1,'_',a$A_2$Var2)
apply(all,1,function(x){
  x[c('E_1','H_1','K_1')]->SSC
  x[c('F_1','I_1','J_1','M_1')]->LS
  wilcox.test(SSC,LS)->test
  return(c(p=test$p.value,diff=mean(SSC)-mean(LS)))
})->pvalues
pvalues<-as.data.frame(t(pvalues))

all[which(rownames(all)%in%as.vector(outer(target_cell[1:2], target_cell[3:5], paste, sep = "_"))),]
pvalues[which(rownames(pvalues)%in%as.vector(outer(target_cell[3:5],target_cell[1:2], paste, sep = "_"))),]



