

#read data

#merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")
merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/immune_proximal/merged_filter_exclude_epi_endo_fib.rds")
library(Matrix)
library(parallel)

root <- "/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/"
all_file <- list.files(root)
all_file <- all_file[nchar(all_file) == 3]

# ---- Fast coercion helpers ----
to_dgC_fast <- function(M) {
  # Common cases first
  if (inherits(M, "dgCMatrix")) return(M)
  if (inherits(M, "dgTMatrix")) return(as(M, "dgCMatrix"))
  if (inherits(M, "dsCMatrix")) {
    # expand symmetric CSC minimally
    return(as(M + t(M) - Diagonal(x = diag(M)), "dgCMatrix"))
  }
  if (inherits(M, "dsTMatrix")) {
    Mt <- as(M, "dgTMatrix")
    Mt_full <- Mt + t(Mt) - Diagonal(x = diag(Mt))
    return(as(Mt_full, "dgCMatrix"))
  }
  as(M, "dgCMatrix")
}

read_adj_sparse <- function(path, nodes) {
  M <- readMM(path)
  if (length(M) == 0L || any(dim(M) == 0L)) {
    n <- length(nodes)
    M <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
  } else {
    M <- to_dgC_fast(M)
  }
  dimnames(M) <- list(nodes, nodes)
  M
}

# ---- Build G and counts once per folder ----
build_G_and_counts <- function(cell_type_vec) {
  # Ensure factor with stable levels (NAs -> "UNK" so rows aren’t dropped)
  f <- factor(replace(cell_type_vec, is.na(cell_type_vec), "UNK"))
  k <- nlevels(f)
  
  # Build G with i/j indices (faster than formula interface)
  n <- length(f)
  G <- sparseMatrix(
    i = seq_len(n),
    j = as.integer(f),
    x = 1.0,
    dims = c(n, k),
    dimnames = list(NULL, levels(f))
  )
  
  counts <- as.numeric(tabulate(as.integer(f), nbins = k))
  list(G = G, counts = counts, levels = levels(f))
}

# ---- Aggregate using prebuilt G ----
agg_sum_and_mean <- function(A, G, counts, group_names,cell_names) {
  # B = A %*% G  (n x k)
  B <- A %*% G
  dimnames(B)<-list(cell_names, group_names)
  # sum_mat = G' %*% B  (k x k) but via crossprod for speed
  sum_mat <- crossprod(G, B)
  dimnames(sum_mat) <- list(group_names, group_names)
  recei_sc <- crossprod(G, A)
  dimnames(recei_sc) <- list(group_names, cell_names)
  
  
  # mean_mat = D^{-1} * sum_mat * D^{-1}
  Dinv <- Diagonal(x = 1 / counts)
  mean_mat <- Dinv %*% sum_mat %*% Dinv
  dimnames(mean_mat) <- list(group_names, group_names)
  
  list(sum = sum_mat, mean = mean_mat,send_sc=B, recei_sc=recei_sc,original=A)
}

# ---- Main loop with per-folder precompute + per-file parallel ----
result <- vector("list", length(all_file))
names(result) <- all_file

# Precompute a name->celltype map once (used across folders)
ct_map <- merged_filter_exclude_epi$cell_type_l1_all
names(ct_map) <- colnames(merged_filter_exclude_epi)

# choose cores (adjust to your node)
n_cores <- max(1L, floor(detectCores() * 0.8))

for (i in seq_along(all_file)) {
  pa <- file.path(root, all_file[i])
  
  nodes <- scan(file.path(pa, "nodes.tsv"), what = character(), quiet = TRUE)
  mtx_files <- list.files(pa, pattern = "\\.mtx$", full.names = TRUE)
  
  # cell types aligned to nodes
  cell_type_vec <- unname(ct_map[nodes])
  
  # Build G, counts once
  GC <- build_G_and_counts(cell_type_vec)
  G <- GC$G
  counts <- GC$counts
  gnames <- GC$levels
  
  # parallel over files
  mats_i <- mclapply(
    mtx_files,
    mc.cores = n_cores,
    FUN = function(f) {
      A <- read_adj_sparse(f, nodes)
      ag <- agg_sum_and_mean(A, G, counts, gnames,rownames(A))
      # free A quickly by not returning it
      ag
    }
  )
  names(mats_i) <- sub("\\.mtx$", "", basename(mtx_files))
  result[[i]] <- mats_i
}

## speed-up version
recei_list<-lapply(result, function(x){
  lapply(x, function(y){
    return(y[['original']])
  })
})

#permute test
pair_all <- list()
all_list_recei_sender <- readRDS("~/all_list_recei_sender.rds")
LR_pairs<-Reduce(union,lapply(all_list_recei_sender, function(x) Reduce(union,lapply(x,rownames))))

samples<-unique(merged_filter_exclude_epi$orig.ident)
library(parallel)
res <- mclapply(LR_pairs, mc.cores=20,function(LR_pair) {
  pair1 <- vector("list", length(samples))
  names(pair1) <- samples
  for (sample in samples) {
    M <- recei_list[[sample]][[LR_pair]]
    group <- merged_filter_exclude_epi$cell_type_l1_all[merged_filter_exclude_epi$orig.ident == sample]
    
    # make sure lengths match
    if (length(group) != nrow(M)) {
      stop("Length of group does not match nrow(M) for sample ", sample)
    }
    
    ## unique cell types and precomputed indices
    celltype2 <- unique(group)
    f <- factor(group, levels = celltype2)
    idx_by_g <- split(seq_along(group), f)  # list of indices per cell type
    G <- length(celltype2)
    
    ## observed means matrix (G x G)
    obs <- matrix(NA_real_, nrow = G, ncol = G,
                  dimnames = list(celltype2, celltype2))
    for (i in seq_len(G)) {
      A <- idx_by_g[[i]]
      for (j in seq_len(G)) {
        B <- idx_by_g[[j]]
        obs[i, j] <- mean(M[A, B])
      }
    }
    
    ## permutation test: count how many permuted >= observed
    set.seed(100)
    nperm <- 10000L
    ge_count <- matrix(0L, nrow = G, ncol = G)
    
    n <- length(group)
    
    for (p in seq_len(nperm)) {
      # permute cell indices (equivalent to permuting labels under the null)
      perm <- sample.int(n)
      
      # permuted indices per group:
      # P_g[[k]] = permuted indices of cells in group k
      P_g <- lapply(idx_by_g, function(idx) perm[idx])
      
      # compute permuted means and update counts directly
      for (i in seq_len(G)) {
        A_perm <- P_g[[i]]
        for (j in seq_len(G)) {
          B_perm <- P_g[[j]]
          val <- mean(M[A_perm, B_perm])
          if (!is.na(val) && val >= obs[i, j]) {
            ge_count[i, j] <- ge_count[i, j] + 1L
          }
        }
      }
    }
    
    ## p-value matrix
    pval <- ge_count / nperm
    rownames(pval) <- celltype2
    colnames(pval) <- celltype2
    pair1[[sample]] <- pval
  }
  pair1
})
names(res) <- LR_pairs
saveRDS(res,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/intermediate_result/permutation_test/res_inf.rds')
-asd-asd-asd

#res <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/intermediate_result/permutation_test/res_G_1.rds")
res <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/intermediate_result/permutation_test/res_inf_G_1.rds")
combine_p_stouffer <- function(p, w = NULL) {
  p <- pmin(pmax(p, .Machine$double.xmin), 1 - 1e-16)  # avoid Inf
  z <- qnorm(1 - p)  # one-sided -> Z
  if (is.null(w)) w <- rep(1, length(p))
  Z <- sum(w * z) / sqrt(sum(w^2))
  p_combined <- 1 - pnorm(Z)  # one-sided combined p
  p_combined
}
fisher_combine_p <- function(p, na.rm = TRUE) {
  p <- as.numeric(p)
  if (na.rm) p <- p[!is.na(p)]
  if (length(p) == 0) stop("No p-values provided after removing NA.")
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("All p-values must be between 0 and 1.")
  
  # Avoid -Inf when p = 0
  p <- pmax(p, .Machine$double.xmin)
  
  stat <- -2 * sum(log(p))
  df <- 2 * length(p)
  p_comb <- pchisq(stat, df = df, lower.tail = FALSE)
  p_comb
}

lapply(res, function(x) {
  lapply(x, function(y){
    df_long <- as.data.frame(as.table(y))
    colnames(df_long) <- c("sender", "receiver", "p_value")
    return(df_long)
  })->list_signal
  Reduce(intersect,lapply(list_signal,function(y) paste0(y$sender,'-',y$receiver)))->common
  lapply(list_signal, function(y) {
    match(common,paste0(y$sender,'-',y$receiver))->order
    y[order,3]
  })->list_signal
  do.call(cbind,list_signal)->mtx_signal
  rownames(mtx_signal)<-common
  as.data.frame(mtx_signal)->mtx_signal
  apply(mtx_signal, 1, combine_p_stouffer)->mtx_signal$combine
  apply(mtx_signal[,c('E_1','H_1','K_1')], 1, combine_p_stouffer)->mtx_signal$combine_ssc
  apply(mtx_signal[,c('F_1','I_1','J_1','M_1')], 1, combine_p_stouffer)->mtx_signal$combine_ls
  apply(mtx_signal[,c('combine','combine_ssc','combine_ls')], 2, function(x) {p.adjust(x,method = 'fdr')})->mtx_signal[,c('combine_adj','combine_ssc_adj','combine_ls_adj')]
  return(mtx_signal)
})->res1_stouffer

lapply(res, function(x) {
  lapply(x, function(y){
    df_long <- as.data.frame(as.table(y))
    colnames(df_long) <- c("sender", "receiver", "p_value")
    return(df_long)
  })->list_signal
  Reduce(intersect,lapply(list_signal,function(y) paste0(y$sender,'-',y$receiver)))->common
  lapply(list_signal, function(y) {
    match(common,paste0(y$sender,'-',y$receiver))->order
    y[order,3]
  })->list_signal
  do.call(cbind,list_signal)->mtx_signal
  rownames(mtx_signal)<-common
  as.data.frame(mtx_signal)->mtx_signal
  apply(mtx_signal, 1, fisher_combine_p)->mtx_signal$combine
  apply(mtx_signal[,c('E_1','H_1','K_1')], 1, fisher_combine_p)->mtx_signal$combine_ssc
  apply(mtx_signal[,c('F_1','I_1','J_1','M_1')], 1, fisher_combine_p)->mtx_signal$combine_ls
  apply(mtx_signal[,c('combine','combine_ssc','combine_ls')], 2, function(x) {p.adjust(x,method = 'fdr')})->mtx_signal[,c('combine_adj','combine_ssc_adj','combine_ls_adj')]
  return(mtx_signal)
})->res1_fisher


res1_name<-lapply(res1_fisher, rownames)
all(sapply(res1_name[-1], function(x) identical(x, res1_name[[1]])))


#The first level of all_list is the first word
#the first word is the category which is on the cell level
all_list_sender_recei <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/intermediate_result/compare_LS_SSC/all_list_inf_sender_recei.rds")
all_list_recei_sender <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/intermediate_result/compare_LS_SSC/all_list_inf_recei_sender.rds")

#let's start with sender on single cell level
res1_stouffer->res1
lapply(names(all_list_sender_recei),function(x){
  x1<-all_list_sender_recei[[x]]
  lapply(names(x1), function(y){
    x1[[y]]->df
    df1<-do.call(rbind,lapply(res1[rownames(df)],function(z) z[paste(x,y,sep = '-'),c('combine','combine_ssc','combine_ls')]))
    df<-as.data.frame(cbind(df,df1))
    df$p_adj<-p.adjust(df$p,method = 'fdr')
    return(df)
  })->list_all
  names(list_all)<-names(x1)
  return(list_all)
})->all_list_sender_recei1
names(all_list_sender_recei1)<-names(all_list_sender_recei)

lapply(names(all_list_recei_sender),function(x){
  x1<-all_list_recei_sender[[x]]
  lapply(names(x1), function(y){
    x1[[y]]->df
    df1<-do.call(rbind,lapply(res1[rownames(df)],function(z) z[paste(y,x,sep = '-'),c('combine','combine_ssc','combine_ls')]))
    df<-as.data.frame(cbind(df,df1))
    df$p_adj<-p.adjust(df$p,method = 'fdr')
    return(df)
  })->list_all
  names(list_all)<-names(x1)
  return(list_all)
})->all_list_recei_sender1
names(all_list_recei_sender1)<-names(all_list_recei_sender)

#compare the numbers between LS and SSC
res1_stouffer->res1
lapply(colnames(res1[["commot-cellchat-CDH"]]), function(x){
  lapply(res1,function(y){
    y[,x]
  })->list
  do.call(rbind,list)->list
  colnames(list)<-rownames(res1[["commot-cellchat-CDH"]])
  return(list)
})->a
names(a)<-colnames(res1[["commot-cellchat-CDH"]])


lapply(a, function(x){
  apply(x,2, function(y) length(which(y<0.05)))
})->a_num
a_num_name<-lapply(a_num, names)
all(sapply(a_num_name[-1], function(x) identical(x, a_num_name[[1]])))
as.data.frame(do.call(rbind,a_num))->a_num
group <- sapply(strsplit(colnames(a_num), "-"), function(z) {
  paste(sort(z), collapse = "-")
})
s<-apply(a_num, 1, function(x) aggregate(x,by=list(group),sum))
cbind(s[["combine_ssc_adj"]],s[["combine_ls_adj"]][,-1])->ls_ssc
colnames(ls_ssc)<-c('cell_pair','ssc','ls')
ls_ssc$ls_ssc<-ls_ssc$ls-ls_ssc$ssc
df_num<-ls_ssc[which(ls_ssc$cell_pair%in%c('Fibroblasts-Macrophage','Fibroblasts-T_cell','Endothelium-Macrophage','Endothelium-T_cell','B_Plasma-Endothelium','B_Plasma-Fibroblasts')),]
df_num1<-ls_ssc[which(ls_ssc$cell_pair%in%c('inf_fib-Macrophage','inf_fib-T_cell','inf_endo-Macrophage','inf_endo-T_cell','B_Plasma-inf_endo','B_Plasma-inf_fib')),]
View(cbind(df_num,df_num1))
df_num_inf_ssc<-data.frame(inf_fib=c(30,4,3),inf_endo=c(25,5,2),row.names = c('Macrophage','T_cell','B_Plasma'))
df_num_inf_ls<-data.frame(inf_fib=c(40,23,14),inf_endo=c(16,7,3),row.names = c('Macrophage','T_cell','B_Plasma'))


result<-result[which(names(result)!='G_1')]
result_sum<-lapply(result,function(x) {Reduce(`+`,lapply(x,function(y) {as.matrix(y$sum)}))})
result_sum<-lapply(result_sum, function(x) {
  as.data.frame(as.table(x))->df
  rownames(df)<-paste(df$Var1,df$Var2,sep='-')
  return(df)
})
Reduce(intersect,lapply(result_sum,function(x) rownames(x)))->common
do.call(cbind,lapply(result_sum,function(x) x[common,3]))->result_sum1
rownames(result_sum1)<-common
key <- sapply(strsplit(rownames(result_sum1), "-"),
              function(z) paste(sort(z), collapse = "-"))
result_sum1_1 <- rowsum(result_sum1, group = key)
t(scale(t(result_sum1_1)))->result_sum2
df<-data.frame(SSC=rowMeans(result_sum2[,c('E_1','H_1','K_1')]),LS=rowMeans(result_sum2[,c('F_1','I_1','J_1','M_1')]))
#df<-data.frame(SSC=rowMeans(result_sum1[,c('E_1','H_1','K_1')]),LS=rowMeans(result_sum1[,c('F_1','I_1','J_1','M_1')]))
#df<-data.frame(
  #SSC=c((result_sum1[grep('inf_fib',rownames(result_sum1)),c('E_1','H_1','K_1')]%*%c(1/99,1/92,1/123))/3,(result_sum1[which(grepl('inf_endo',rownames(result_sum1))&!grepl('inf_fib',rownames(result_sum1))),c('E_1','H_1','K_1')]%*%c(1/64,1/63,1/77))/3),
  #LS=c((result_sum1[grep('inf_fib',rownames(result_sum1)),c('F_1','I_1','J_1','M_1')]%*%c(1/708,1/31,1/22,1/922))/3,(result_sum1[which(grepl('inf_endo',rownames(result_sum1))&!grepl('inf_fib',rownames(result_sum1))),c('F_1','I_1','J_1','M_1')]%*%c(1/206,1/11,1/15,1/293))/3),
  #row.names = c(rownames(result_sum1)[grep('inf_fib',rownames(result_sum1))],rownames(result_sum1)[which(grepl('inf_endo',rownames(result_sum1))&!grepl('inf_fib',rownames(result_sum1)))])
#)
#df$SSC<-df$SSC/3
#df$LS<-df$LS/4

df_num_inf_ssc_streng<-data.frame(inf_fib=df[index_inf,'SSC'][c(1,2,6)],inf_endo=df[index_inf,'SSC'][c(3,4,5)],row.names = c('Macrophage','T_cell','B_Plasma'))
df_num_inf_ls_streng<-data.frame(inf_fib=df[index_inf,'LS'][c(1,2,6)],inf_endo=df[index_inf,'LS'][c(3,4,5)],row.names = c('Macrophage','T_cell','B_Plasma'))


#later analyses (this is trash code)
df_para_immune1<-df[c('Fibroblasts-Macrophage','Fibroblasts-T_cell','Endothelium-Macrophage','Endothelium-T_cell','Endothelium-B_Plasma','Fibroblasts-B_Plasma'),]
df_immune_para1<-df[c('Macrophage-Fibroblasts','T_cell-Fibroblasts','Macrophage-Endothelium','T_cell-Endothelium','B_Plasma-Endothelium','B_Plasma-Fibroblasts'),]
df_para_immune2<-df[c('inf_fib-Macrophage','inf_fib-T_cell','inf_endo-Macrophage','inf_endo-T_cell','inf_endo-B_Plasma','inf_fib-B_Plasma'),]
df_immune_para2<-df[c('Macrophage-inf_fib','T_cell-inf_fib','Macrophage-inf_endo','T_cell-inf_endo','B_Plasma-inf_endo','B_Plasma-inf_fib'),]

(df_immune_para1+df_para_immune1)/2->a
(df_immune_para2+df_para_immune2)/2->b
as.data.frame(cbind(a,b))->mtx_num
colnames(mtx_num)<-c('SSc_nor','LS_nor','SSc_inf','LS_inf')
View(mtx_num)

pheatmap::pheatmap((df_immune_para+df_para_immune)/2)
pheatmap::pheatmap(df_immune_para)
pheatmap::pheatmap(df_para_immune)


#summarized comparison between inflammatory cells and other cells (this is trash code)
index_inf<-c('inf_fib-Macrophage','inf_fib-T_cell','inf_endo-Macrophage','inf_endo-T_cell','B_Plasma-inf_endo','B_Plasma-inf_fib')
index_non_inf<-c('Fibroblasts-Macrophage','Fibroblasts-T_cell','Endothelium-Macrophage','Endothelium-T_cell','B_Plasma-Endothelium','B_Plasma-Fibroblasts')
non_inf<-a[["combine_adj"]][,index_non_inf]
inf<-a[["combine_adj"]][,index_inf]
apply(non_inf, 2, function(x) length(which(x<0.05)))
df_non_inf<-data.frame(fib=c(1,1,0),endo=c(1,0,0),row.names = c('Macrophage','T_cell','B_Plasma'))
apply(inf, 2, function(x) length(which(x<0.05)))
df_inf<-data.frame(fib=c(31,20,14),endo=c(25,10,1),row.names = c('Macrophage','T_cell','B_Plasma'))
#pheatmap::pheatmap(df_inf,breaks = seq(0,31,length.out=100),cluster_rows = F,cluster_cols = F,border_color = NA)->p_inf
#pheatmap::pheatmap(df_non_inf,breaks = seq(0,31,length.out=100),cluster_rows = F,cluster_cols = F,border_color = NA)->p_non_inf
#ggsave('p_inf.png',p_inf,width = 3.5,height = 3.5,dpi = 300)
#ggsave('p_non_inf.png',p_non_inf,width = 3.5,height = 3.5,dpi = 300)

result_sum1->result_sum2
#df_inf_streng<-data.frame(fib=as.numeric(rowMeans(result_sum2[index_inf,])[c(1,2,6)]),endo=as.numeric(rowMeans(result_sum2[index_inf,])[c(3,4,5)]),row.names = c('Macrophage','T_cell','B_Plasma'))
#pheatmap::pheatmap(df_inf[,1],breaks = seq(min(df_non_inf[,1]),max(df_inf[,1]),length.out=100),border_color = NA,cluster_rows = F,cluster_cols = F)->p
#pheatmap::pheatmap(df_inf[,2],breaks = seq(min(df_non_inf[,2]),max(df_inf[,2]),length.out=100),border_color = NA,cluster_rows = F,cluster_cols = F)->p
#t(scale(t(result_sum1_1))) ->result_sum2
df_inf_streng<-data.frame(fib=as.numeric(result_sum2[index_inf,][c(1,2,6)]),endo=as.numeric(result_sum2[index_inf,][c(3,4,5)]),row.names = c('Macrophage','T_cell','B_Plasma'))
df_non_inf_streng<-data.frame(fib=as.numeric(result_sum2[index_non_inf,][c(1,2,6)]),endo=as.numeric(result_sum2[index_non_inf,][c(3,4,5)]),row.names = c('Macrophage','T_cell','B_Plasma'))
pheatmap::pheatmap(df_non_inf[,1],breaks = seq(min(df_non_inf[,1]),max(df_inf[,1]),length.out=100),border_color = NA,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(df_non_inf[,2],breaks = seq(min(df_non_inf[,2]),max(df_inf[,2]),length.out=100),border_color = NA,cluster_rows = F,cluster_cols = F)


#find the signals only significant in LS but not in SSc
res1<-res1_stouffer
labels_F_M<-unlist(lapply(res1, function(x){x['Fibroblasts-Macrophage','combine_ls_adj']<0.05&x['Fibroblasts-Macrophage','combine_ssc_adj']>0.05}))
which(labels_F_M)
labels_F_M<-unlist(lapply(res1, function(x){x['inf_fib-Macrophage','combine_ls_adj']<0.05&x['inf_fib-Macrophage','combine_ssc_adj']>0.05}))
which(labels_F_M)

labels_M_F<-unlist(lapply(res1, function(x){x['Macrophage-Fibroblasts','combine_ls_adj']<0.05&x['Macrophage-Fibroblasts','combine_ssc_adj']>0.05}))
which(labels_M_F)
labels_M_F<-unlist(lapply(res1, function(x){x['Macrophage-inf_fib','combine_ls_adj']<0.05&x['Macrophage-inf_fib','combine_ssc_adj']>0.05}))
which(labels_M_F)
intersect(names(labels_F_M)[which(labels_F_M)],names(labels_M_F)[which(labels_M_F)])

res1<-res1_stouffer
labels_F_T<-unlist(lapply(res1, function(x){x['Fibroblasts-T_cell','combine_ls_adj']<0.05&x['Fibroblasts-T_cell','combine_ssc_adj']>0.05}))
which(labels_F_T)
labels_F_T<-unlist(lapply(res1, function(x){x['inf_fib-T_cell','combine_ls_adj']<0.05&x['inf_fib-T_cell','combine_ssc_adj']>0.05}))
which(labels_F_T)

labels_T_F<-unlist(lapply(res1, function(x){x['T_cell-Fibroblasts','combine_ls_adj']<0.05&x['T_cell-Fibroblasts','combine_ssc_adj']>0.05}))
which(labels_T_F)
labels_T_F<-unlist(lapply(res1, function(x){x['T_cell-inf_fib','combine_ls_adj']<0.05&x['T_cell-inf_fib','combine_ssc_adj']>0.05}))
which(labels_T_F)
intersect(names(labels_F_T)[which(labels_F_T)],names(labels_T_F)[which(labels_T_F)])

#macrophage to fibroblasts: c("commot-cellchat-PDGFB-PDGFRB", "commot-cellchat-CXCL12-CXCR4")
#fibroblasts to macrophage: c("commot-cellchat-COL6A3-CD44", "commot-cellchat-CXCL12-CXCR4","commot-cellchat-COMP-CD47","commot-cellchat-FN1-CD44","commot-cellchat-IL34-CSF1R")
#T_cell to fibroblasts: c("commot-cellchat-PDGFB-PDGFRB", "commot-cellchat-CXCL12-CXCR4")
#fibroblasts to T_cell: c("commot-cellchat-COL1A2-CD44", "commot-cellchat-COL1A1-CD44","commot-cellchat-COMP-CD47","commot-cellchat-CXCL12-CXCR4","commot-cellchat-LAMC1-CD44","commot-cellchat-LAMA4-CD44","commot-cellchat-LAMA2-CD44","commot-cellchat-FN1-CD44","commot-cellchat-COL6A1-CD44")

M_F<-c("commot-cellchat-PDGFB-PDGFRB", "commot-cellchat-CXCL12-CXCR4")
View(all_list_sender_recei1$Macrophage$inf_fib[M_F,])
View(all_list_recei_sender1$inf_fib$Macrophage[M_F,])
F_M<-c("commot-cellchat-COL6A3-CD44", "commot-cellchat-CXCL12-CXCR4","commot-cellchat-COMP-CD47","commot-cellchat-FN1-CD44","commot-cellchat-IL34-CSF1R")
View(all_list_sender_recei1$inf_fib$Macrophage[F_M,])
View(all_list_recei_sender1$Macrophage$inf_fib[F_M,])

T_F<-c("commot-cellchat-PDGFB-PDGFRB", "commot-cellchat-CXCL12-CXCR4")
View(all_list_sender_recei1$T_cell$inf_fib[T_F,])
View(all_list_recei_sender1$inf_fib$T_cell[T_F,])
F_T<-c("commot-cellchat-COL1A2-CD44", "commot-cellchat-COL1A1-CD44","commot-cellchat-COMP-CD47","commot-cellchat-CXCL12-CXCR4","commot-cellchat-LAMC1-CD44","commot-cellchat-LAMA4-CD44","commot-cellchat-LAMA2-CD44","commot-cellchat-FN1-CD44","commot-cellchat-COL6A1-CD44")
View(all_list_sender_recei1$inf_fib$T_cell[F_T,])
View(all_list_recei_sender1$T_cell$inf_fib[F_T,])

pheatmap::pheatmap(t(all_list$inf_fib$Macrophage[M_F,c('ls_scale','ssc_scale')]),cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(t(all_list$Macrophage$inf_fib[F_M,c('ls_scale','ssc_scale')]),cluster_rows = F,cluster_cols = F)

pheatmap::pheatmap(t(all_list$inf_fib$T_cell[T_F,c('ls_scale','ssc_scale')]),cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(t(all_list$T_cell$inf_fib[F_T,c('ls_scale','ssc_scale')]),cluster_rows = F,cluster_cols = F)

unique(c(F_M,M_F,T_F,F_T))->all_lr
result_select<-lapply(result, function(x) {
  lapply(x[all_lr],function(y) y$original)
})
result_select<-result_select[which(names(result_select)!='G_1')]
result_select<-lapply(names(result_select$A_2), function(x) {
  lapply(result_select, function(y) y[[x]])
})
names(result_select)<-all_lr
SSC_index<-c('E_1','H_1','K_1')
LS_index<-c('F_1','I_1','J_1','M_1')

sender<-'inf_fib'
sender_index<-colnames(merged_filter_base_exclude_epi_LS_SSC)[which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==sender)]
list_result_all<-list()
for (receiver in c('Macrophage','T_cell')) {
  receiver_index<-colnames(merged_filter_base_exclude_epi_LS_SSC)[which(merged_filter_base_exclude_epi_LS_SSC$cell_type_l1_all==receiver)]
  lapply(result_select, function(y) {
    lapply(y, function(x) {
      x[which(rownames(x)%in%sender_index),which(colnames(x)%in%receiver_index)]
    })->list_all
    unlist(lapply(list_all[SSC_index], function(m) {as.vector(as.matrix(m))}))->ssc
    unlist(lapply(list_all[LS_index], function(m) {as.vector(as.matrix(m))}))->ls
    if (length(ssc)==0|length(ls)==0) {
      c(stat=NA,P=NA,diff=mean(ssc)-mean(ls),ssc=mean(ssc),ls=mean(ls),log2fc=log2(mean(ssc)/mean(ls))) 
    } else {
      wilcox.test(ssc,ls)->test
      c(stat=test$statistic,P=test$p.value,ssc=mean(ssc),ls=mean(ls),diff=mean(ssc)-mean(ls),log2fc=log2(mean(ssc)/mean(ls))) 
    }
  })->list_result
  do.call(rbind,list_result)->list_result
  list_result->list_result_all[[receiver]]
}


#-asd-asd-asd (trash code)
for (LR_pair in LR_pairs) {
  pair1 <- list()
  for (sample in unique(merged_filter_exclude_epi$orig.ident)) {
    
    M <- recei_list[[sample]][[LR_pair]]
    group <- merged_filter_exclude_epi$cell_type_l1_all[merged_filter_exclude_epi$orig.ident == sample]
    
    # make sure lengths match
    if (length(group) != nrow(M)) {
      stop("Length of group does not match nrow(M) for sample ", sample)
    }
    
    ## unique cell types and precomputed indices
    celltype2 <- unique(group)
    f <- factor(group, levels = celltype2)
    idx_by_g <- split(seq_along(group), f)  # list of indices per cell type
    G <- length(celltype2)
    
    ## observed means matrix (G x G)
    obs <- matrix(NA_real_, nrow = G, ncol = G,
                  dimnames = list(celltype2, celltype2))
    for (i in seq_len(G)) {
      A <- idx_by_g[[i]]
      for (j in seq_len(G)) {
        B <- idx_by_g[[j]]
        obs[i, j] <- mean(M[A, B])
      }
    }
    
    ## permutation test: count how many permuted >= observed
    set.seed(100)
    nperm <- 10000L
    ge_count <- matrix(0L, nrow = G, ncol = G)
    
    n <- length(group)
    
    for (p in seq_len(nperm)) {
      # permute cell indices (equivalent to permuting labels under the null)
      perm <- sample.int(n)
      
      # permuted indices per group:
      # P_g[[k]] = permuted indices of cells in group k
      P_g <- lapply(idx_by_g, function(idx) perm[idx])
      
      # compute permuted means and update counts directly
      for (i in seq_len(G)) {
        A_perm <- P_g[[i]]
        for (j in seq_len(G)) {
          B_perm <- P_g[[j]]
          val <- mean(M[A_perm, B_perm])
          if (!is.na(val) && val >= obs[i, j]) {
            ge_count[i, j] <- ge_count[i, j] + 1L
          }
        }
      }
    }
    
    ## p-value matrix
    pval <- ge_count / nperm
    rownames(pval) <- celltype2
    colnames(pval) <- celltype2
    pair1[[sample]] <- pval
  }
  
  pair_all[[LR_pair]] <- pair1
}
saveRDS(pair_all,'/ix1/wchen/liutianhao/pair_M_F_complement.rds')



-asd-asd-asd
lapply(pair_all, function(x) {
  lapply(x, function(y) {
    df<-as.data.frame(y) %>%
      mutate(row = rownames(y)) %>%        # keep row names
      pivot_longer(
        cols = -row,                        # all other columns become "column names"
        names_to = "col",                   # name of column-name variable
        values_to = "value"                 # name of value variable
      )
    as.data.frame(df)->df
    rownames(df)<-paste(df$row,df$col,sep = '_')
    return(df)
  })
})->result_significance

lapply(result_significance, function(x) {
  Reduce(intersect,lapply(x,rownames))->common_cell
  do.call(cbind,lapply(x, function(y) as.data.frame(y[common_cell,3])))->df
  colnames(df)<-names(x)
  rownames(df)<-common_cell
  return(df)
})->common_cell

combine_fisher <- function(pvals) {
  X <- -2 * sum(log(pvals))
  p <- pchisq(X, df = 2 * length(pvals), lower.tail = FALSE)
  return(p)
}

lapply(common_cell, function(x) {
  x<-x[,which(colnames(x)!='G_1')]
  apply(x, 1, combine_fisher)->x$combined
  x<-x[order(x$combined),]
  return(x)
})->common_cell

lapply(common_cell, function(x) {
  x[c('Fibroblasts_Macrophage','Macrophage_Fibroblasts'),]
})->common_cell_part
do.call(rbind,common_cell_part)->common_cell_part
rownames(common_cell_part)[which(common_cell_part$combined<0.05)]


do.call(cbind,apply(common_cell_part[,grep('1',colnames(common_cell_part))],1,function(x) {
  aggregate(x,by=list(merged_filter_exclude_epi$condition[match(colnames(common_cell_part)[grep('1',colnames(common_cell_part))],merged_filter_exclude_epi$orig.ident)]),mean)
}))->all



