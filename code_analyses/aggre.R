library(Seurat)
list.files('/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1')->sample
spot_to_polygon<-list()
for (i in sample[grep('result',sample)]) {
  i<-substr(i,8,10)
  spot_to_polygon[[i]]<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/',i,'/',i,'_spot_to_polygon.csv'))
}


#lung data
spot_to_polygon<-list()
for (i in 0:4) {
  for (j in 0:4) {
    mtx<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_neo/lung_neo',i,j,'_spot_to_polygon.csv'))
    if (dim(mtx)[1]!=0) {
      mtx$polygon_id<-paste0(mtx$polygon_id,'_',i,j)
    } 
    spot_to_polygon[[paste0(i,j)]]<-mtx
  }
}
do.call(rbind,spot_to_polygon)->spot_to_polygon
aggregate_pro<-function(route){
  library(Matrix)
  library(Matrix.utils)
  mtx<-Read10X_h5(route)
  t(mtx)->mtx
  mtx<-mtx[spot_to_polygon$barcode,]
  aggregate.Matrix(mtx, groupings = setNames(spot_to_polygon$polygon_id, spot_to_polygon$barcode), fun = "sum")->mtx_sc
  mtx_sc<-t(mtx_sc)
  mtx_sc<-CreateSeuratObject(mtx_sc)
  NormalizeData(mtx_sc)->mtx_sc
  FindVariableFeatures(mtx_sc)->mtx_sc
  return(mtx_sc)
}
mtx_sc_pro_neo<-aggregate_pro('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/8743_nA1/square_002um/filtered_feature_bc_matrix.h5')
saveRDS(mtx_sc_pro,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_neo/result/st_neo_raw_all.rds')
#sc_pro(mtx_sc_pro,FALSE)->mtx_sc_pro
mtx_sc_pro_filter<-mtx_sc_pro[,which(mtx_sc_pro$nCount_RNA>100)]
saveRDS(mtx_sc_pro_filter,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_neo/st_neo_raw_filtered.rds')
sc_pro(mtx_sc_pro_filter,FALSE)->mtx_sc_pro_filter



asdasdasd
#crc data
spot_to_polygon<-list()
for (i in 0:4) {
  for (j in 0:4) {
    mtx<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/crc/crc',i,j,'_spot_to_polygon.csv'))
    if (dim(mtx)[1]!=0) {
      mtx$polygon_id<-paste0(mtx$polygon_id,'_',i,j)
    } 
    spot_to_polygon[[paste0(i,j)]]<-mtx
  }
}
do.call(rbind,spot_to_polygon)->spot_to_polygon
aggregate_pro<-function(route){
  library(Matrix)
  library(Matrix.utils)
  mtx<-Read10X_h5(route)
  t(mtx)->mtx
  mtx<-mtx[spot_to_polygon$barcode,]
  aggregate.Matrix(mtx, groupings = setNames(spot_to_polygon$polygon_id, spot_to_polygon$barcode), fun = "sum")->mtx_sc
  mtx_sc<-t(mtx_sc)
  mtx_sc<-CreateSeuratObject(mtx_sc)
  NormalizeData(mtx_sc)->mtx_sc
  FindVariableFeatures(mtx_sc)->mtx_sc
  return(mtx_sc)
}
mtx_sc_pro<-aggregate_pro('/ix1/wchen/liutianhao/data/public_lung_spatial/CRC_10X/HD/P1_CRC/binned_outputs/square_002um/filtered_feature_bc_matrix.h5')
sc_pro(mtx_sc_pro,FALSE)->mtx_sc_pro
mtx_sc_pro_filter<-mtx_sc_pro[,which(mtx_sc_pro$nCount_RNA>50)]
sc_pro(mtx_sc_pro_filter,FALSE)->mtx_sc_pro_filter
saveRDS(mtx_sc_pro_filter,'/ix1/wchen/liutianhao/mtx_sc_pro_filter.rds')
DimPlot(mtx_sc_pro,group.by='')
marker<-c(    "LYZ", "STAB1", "C1QC", "MSR1", "C1QB", "IL7R", "CD40LG", "SPOCK2", "CCR4",
              "TRBC1", "CCL5", "KLRK1", "NKG7", "GZMA", "CD8A", "MS4A1", "CD79A", "IGHM",
              "BLK", "BANK1", "THBS2", "MFAP2", "CTHRC1", "AEBP1", "COL12A1")

DoHeatmap(mtx_sc_pro_filter,marker)




#skin data
for (i in sample[grep('result',sample)]) {
  i<-substr(i,8,10)
  spot_to_polygon[[i]]<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/',i,'/',i,'_spot_to_polygon.csv'))
}

#lung neo
spot_to_polygon<-read.csv(paste0('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_neo/lung_neo22_spot_to_polygon.csv'))
library(Matrix)
library(Matrix.utils)
mtx<-Read10X_h5('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/8743_nA1/square_002um/filtered_feature_bc_matrix.h5')
t(mtx)->mtx
mtx<-mtx[spot_to_polygon$barcode,]
aggregate.Matrix(mtx, groupings = setNames(spot_to_polygon$polygon_id, spot_to_polygon$barcode), fun = "sum")->mtx_sc
mtx_sc<-t(mtx_sc)
mtx_sc<-CreateSeuratObject(mtx_sc)
NormalizeData(mtx_sc)->mtx_sc
FindVariableFeatures(mtx_sc)->mtx_sc
mtx_sc<-mtx_sc[,which(mtx_sc$nCount_RNA>100)]
sc_pro(mtx_sc,FALSE)->mtx_sc
DimPlot(mtx_sc)
#normal analysis
aggregate_pro<-function(i){
  library(Matrix)
  library(Matrix.utils)
  mtx<-Read10X_h5(paste0('/ix1/wchen/liutianhao/result/pathology_ST/count_mtx/2um/filtered_feature_bc_matrix_',gsub('_','',i),'.h5') )
  t(mtx)->mtx
  mtx<-mtx[spot_to_polygon[[i]]$barcode,]
  aggregate.Matrix(mtx, groupings = setNames(spot_to_polygon[[i]]$polygon_id, spot_to_polygon[[i]]$barcode), fun = "sum")->mtx_sc
  mtx_sc<-t(mtx_sc)
  mtx_sc<-CreateSeuratObject(mtx_sc)
  NormalizeData(mtx_sc)->mtx_sc
  FindVariableFeatures(mtx_sc)->mtx_sc
  return(mtx_sc)
}
lapply(names(spot_to_polygon), function(x) aggregate_pro(x)) ->aggregated
names(aggregated)<-names(spot_to_polygon)
library(SeuratData)
library(SeuratDisk)
merged<-merge(aggregated[[1]],aggregated[2:length(aggregated)],add.cell.ids = names(aggregated))
merged$orig.ident<-substr(colnames(merged),1,3)
saveRDS(merged,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_unfilter.rds')

merged_filter<-merged_unfilter[,which(merged_unfilter$nCount_RNA>50)]

sc_pro<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda = 4,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  return(object)
}

merged_filter_list<-list()
for (thresh in c(10,20,30,40,50,60,70,80,90,100)) {
  sc_pro(merged_unfilter[,which(merged_unfilter$nCount_RNA>thresh)],harmony = FALSE)->merged_filter_list[[as.character(thresh)]]
}
saveRDS(merged_filter_list,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter_list.rds')
merged_filter_base<-sc_pro(merged_filter_base,FALSE)
merged_filter_base_exclude_epi<-merged_filter_base[,which(!merged_filter_base$epidermis)]
merged_filter_base_exclude_epi<-sc_pro(merged_filter_base_exclude_epi,FALSE)
DimPlot(merged_filter_base_exclude_epi,group.by = 'cell_type_l1_all')
FeaturePlot(merged_filter_base_exclude_epi,'COL1A1')
summary(merged$nCount_RNA[which(merged$seurat_clusters==9)])
FindMarkers(merged,ident.1 = 9,group.by = 'seurat_clusters')->marker_circle
library(dplyr)
merged@meta.data %>% group_by(seurat_clusters) %>%summarize(sum=mean(nCount_RNA))



FeaturePlot(all_unprocessed_8um_spot2,'ACTA2')
FeaturePlot(all_unprocessed_8um_spot2,'CD19')
FeaturePlot(all_unprocessed_8um_spot2,'CD8A')
FeaturePlot(merged,'ACTA2')
FeaturePlot(merged,'IGHG1')
FeaturePlot(merged,'CD19')
FeaturePlot(merged,'CD8A')
DimPlot(merged,group.by = 'orig.ident',pt.size = 0.002)
DimPlot(merged)

