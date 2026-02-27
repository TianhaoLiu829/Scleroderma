merged_filter_exclude_epi <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter_exclude_epi.rds")

FindMarkers(merged_filter_exclude_epi_endo,group.by = 'inf_endo',ident.1 = TRUE,ident.2 = FALSE,logfc.threshold = 0)->marker_inf_endo
res_endo<-DE_wilcox(as.matrix(merged_filter_exclude_epi_endo@assays$RNA@data),merged_filter_exclude_epi_endo$inf_endo,TRUE,FALSE)
res_endo<-res_endo[order(res_endo$pvalue),]
marker_inf_endo<-as.data.frame(cbind(marker_inf_endo,res_endo[match(rownames(marker_inf_endo),rownames(res_endo)),]))
View(marker_inf_endo[which(marker_inf_endo$avg_log2FC<0&marker_inf_endo$p_val<0.05),])

FindMarkers(merged_filter_exclude_epi_fib,group.by = 'inf_fib',ident.1 = TRUE,ident.2 = FALSE,logfc.threshold = 0)->marker_inf_fib
res_fib<-DE_wilcox(as.matrix(merged_filter_exclude_epi_fib@assays$RNA@data),merged_filter_exclude_epi_fib$inf_fib,TRUE,FALSE)
res_fib<-res_fib[order(res_fib$pvalue),]
marker_inf_fib<-as.data.frame(cbind(marker_inf_fib,res_fib[match(rownames(marker_inf_fib),rownames(res_fib)),]))
View(marker_inf_fib[which(marker_inf_fib$avg_log2FC<0-0.5&marker_inf_fib$p_val<0.05),])


# Make sure group is a factor with 2 groups
object<-merged_filter_exclude_epi_fib
object$condition<-object$inf_fib
#object<-merged_filter_exclude_epi
#object<-object[,grep('1',object$orig.ident)]
#object<-object[,which(object$cell_type_l1_all=='Endothelium')]
group <- factor(object$condition)
mat<-as.matrix(object@assays$RNA@data)
g1 <- which(group == levels(group)[1])
g2 <- which(group == levels(group)[2])

# Preallocate results
n_gene <- nrow(mat)
pvals <- numeric(n_gene)
avg1  <- numeric(n_gene)
avg2  <- numeric(n_gene)

# Compute group means (vectorized, fast)
avg1 <- rowMeans(mat[, g1, drop = FALSE])
avg2 <- rowMeans(mat[, g2, drop = FALSE])

pct_1 <- rowMeans(mat[, g1, drop = FALSE]>0)
pct_2 <- rowMeans(mat[, g2, drop = FALSE]>0)
pct_all <- rowMeans(mat>0)

# Compute Wilcoxon p-values (use apply over rows)
pvals <- apply(mat, 1, function(x) {
  wilcox.test(x[g1], x[g2])$p.value
})

# Fold change
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0

# Combine results into data.frame
res <- data.frame(
  gene   = rownames(mat),
  avg_group1 = avg1,
  avg_group2 = avg2,
  pct_1 = pct_1,
  pct_2 = pct_2,
  pct_all = pct_all,
  log2FC = log2fc,
  pvalue = pvals,
  padj   = p.adjust(pvals, method = "fdr")
)

View(res)
res<-res[which(res$pct_all>0.01),]
res$padj2<-p.adjust(res$pvalue,method = 'bonferroni')
res$padj<-p.adjust(res$pvalue,method = 'fdr')
rownames(res)[which(res$log2FC>0&res$padj<0.05)]->features
LS_endo<-FindMarkers(object,group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'LS',ident.2 = 'SSC',only.pos = T,logfc.threshold = 0,min.pct = 0,features = features)

FindMarkers(object,group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'LS',ident.2 = 'SSC')->a
avg1 <- rowMeans(mat[rownames(a), g1, drop = FALSE])
avg2 <- rowMeans(mat[rownames(a), g2, drop = FALSE])
log2fc <- log2((avg1 + 1e-8) / (avg2 + 1e-8))  # avoid dividing by 0
a$avg_log2FC2<-log2fc
Features<-rownames(a)[which(a$avg_log2FC2>0&!grepl('MT-',rownames(a)))]
LS_endo_sc2<-FindMarkers(reference[,which(reference$cell_type_l1=='Endothelium')],group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'LS',ident.2 = 'SSC',only.pos = T,logfc.threshold = 0,min.pct = 0,features = Features)

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", "Fibroblasts" = "Fibroblasts","Myofibroblast"="Smooth_muscle_cell", "Glandular_epithelium" = "Glandular_epithelium", "Keratinocyte" = "Keratinocyte", "Lymphatic_endothelium" = "Endothelium","Mast_cell" = "Myeloid","Melanocyte" = "Melanocyte","Myeloid" = "Myeloid","Neuron" = "Neuron","Pericyte" = "Endothelium","Plasma_cell" = "B_Plasma","Plasma_cell" = "Lymphocyte","prolife_keratinocyte" = "Keratinocyte","Sebaceous_gland"="Keratinocyte","Smooth_muscle_cell" = "Smooth_muscle_cell","T_cell"= "Lymphocyte","Vascular_endothelium"="Endothelium","vascular_SMC"="Smooth_muscle_cell")
reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
LS_endo_sc<-FindMarkers(reference[,which(reference$cell_type_l1=='Endothelium')],group.by="condition",latent.vars = 'orig.ident',test.use = 'MAST',ident.1 = 'LS',ident.2 = 'SSC',only.pos = T,logfc.threshold = 0,min.pct = 0,features = features)
LS_endo_sc$p_val_adj2<-p.adjust(LS_endo_sc$p_val,'fdr')

res_endo_LS<-res[which(res$log2FC>0&res$padj<0.05),]
cbind(res_endo_LS,LS_endo_sc[rownames(res_endo_LS),])->res_endo_LS
res_endo_LS<-res_endo_LS[which(res_endo_LS$p_val_adj<0.05),]
saveRDS(res_endo_LS,'res_endo_LS.rds')
#comparisons between immune-proximal and others#comparisons between immune-proximal and others#comparisons between immune-proximal and others
res_inf_endo <- readRDS("~/res_inf_endo.rds")
dot_plot(merged_filter_exclude_epi_endo,c('CD74','CCL19','LYZ','CXCR4','FOS','FABP4','CSRP1','CAVIN1','COL4A1','COL4A2'),'inf_endo')+coord_flip()+theme(axis.text.x = element_text(size=13))->p
ggsave('p.png',p,width = 6,height = 3,dpi = 300)

res_inf_fib <- readRDS("~/res_inf_fib.rds")
dot_plot(merged_filter_exclude_epi_fib,c('CD74','LYZ','C3','CCL19','CXCL9','COL1A1','COL3A1','COL1A2','SPARC','DSP'),'inf_fib')+coord_flip()+theme(axis.text.x = element_text(size=13))->p
ggsave('p.png',p,width = 6,height = 3,dpi = 300)

dot_plot(merged_filter_exclude_epi_fib,c('CD74','LYZ','C3','CCL19','CXCL9','COL1A1','COL3A1','COL1A2','SPARC','DSP'),'inf_fib')+coord_flip()+theme(axis.text.x = element_text(size=13))->p


DE_sc_fib <- readRDS("~/DE_sc_fib.rds")
DE_sc_endo <- readRDS("~/DE_sc_endo.rds")
res_endo<-res_endo[order(res_endo$pvalue),]
res_endo<-res_endo[which(res_endo$log2FC>0),]
DE_sc_endo[rownames(res_endo)[which(res_endo$pvalue<0.05)],]->inter
res_endo[which(inter$avg_log2FC<0&inter$p_val_adj<0.05),]->res_endo_inter


apply(as.matrix(object@assays$RNA@data[rownames(DE_sc_endo_inter_LS)[1:10],]), 1, aggregate(x,by=list(object$condition),mean))


library(MAST)
library(Matrix)
library(data.table)
mat<-as.matrix(object@assays$RNA@data)
unlist(apply(mat, 1, function(x) length(unique(object$orig.ident[which(x>0)]))))->a
mat<-mat[which(a>1),]


## mat: gene x cell
## cell_meta: rows = cells, columns: group (fixed), sample (random)
cell_meta<-data.frame(group=object$condition,sample=object$orig.ident,row.names = colnames(object))
# Make sure columns/rows line up
stopifnot(all(colnames(mat) == rownames(cell_meta)))

# Feature data (one row per gene)
fdata <- data.frame(
  primerid = rownames(mat),
  row.names = rownames(mat)
)

# Build SingleCellAssay
sca <- FromMatrix(
  exprsArray = mat,
  cData      = cell_meta,
  fData      = fdata
)

# Ensure factors
colData(sca)$group  <- factor(colData(sca)$group)   # variable of interest
colData(sca)$sample <- factor(colData(sca)$sample)  # random effect grouping
zlm_fit <- zlm(
  formula = ~ group  +(1|sample),
  sca     = sca,
  method  = "glmer",   # use glmer to allow random effects
  ebayes  = FALSE,     # empirical Bayes not supported with glmer
  strictConvergence = FALSE
)
