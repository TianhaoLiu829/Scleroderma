library(Seurat)
library(harmony)
library(ggplot2)
library(parallel)
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
  "T_cell"="#66A61E",
  "Mast_cell"="#AD74EA")

reference <- readRDS('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/combine_all.rds')
array<-c("B_cell" = "B_Plasma", 
         "Fibroblasts" = "Fibroblasts",
         "Myofibroblast"="Smooth_muscle_cell", 
         "Glandular_epithelium" = "Glandular_epithelium", 
         "Keratinocyte" = "Keratinocyte", 
         "Lymphatic_endothelium" = "Endothelium",
         "Mast_cell" = "Mast_cell",
         "Melanocyte" = "Melanocyte",
         "Myeloid" = "Myeloid",
         "Neuron" = "Neuron",
         "Pericyte" = "Endothelium",
         "Plasma_cell" = "B_Plasma",
         "prolife_keratinocyte" = "Keratinocyte",
         "Sebaceous_gland"="Keratinocyte",
         "Smooth_muscle_cell" = "Smooth_muscle_cell",
         "T_cell"= "T_cell",
         "Vascular_endothelium"="Endothelium",
         "vascular_SMC"="Smooth_muscle_cell")

reference@meta.data$cell_type_l1<-array[reference$cell_type_high]
FeatureScatter(reference, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


I1<-Read10X_h5('/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/all_data/flex/result/outs/per_sample_outs/NRCOS356/count/sample_filtered_feature_bc_matrix.h5')
reference[rownames(I1),]->reference

reference<-reference[,which(reference$nCount_RNA>50&reference$nCount_RNA<10^5&reference$nFeature_RNA<10^4)]
(colSums(reference@assays$RNA@counts[grep('^MT-',rownames(reference)),])/reference$nCount_RNA)*100->reference$percent.mt
reference<-reference[,which(reference$percent.mt<20)]
reference<-reference[,which(reference$orig.ident!='G2')]  

reference->object
NormalizeData(object)->object
FindVariableFeatures(object,nfeatures=2000)->object
object <- ScaleData(object,features = rownames(object))
set.seed(100)
object <- RunPCA(object,npcs=50)
ElbowPlot(object, ndims = 50,reduction = 'pca')

n.neighbors<-25
min.dist<-0.3
spread<-1
sigma<-0.12
theta<-2
lambda<-4
dims<-18
for (sigma in c(0.1,0.105,0.11,0.115,0.12,0.13)) {
  for (theta in c(1,1.5,2,2.5)) {
    for (lambda in c(1,2,3,4)) {
      #sigma<-0.11
      set.seed(100)
      object <- RunHarmony(object,group.by.vars = 'orig.ident',lambda = lambda,theta = theta,sigma=sigma)
      #set.seed(100)
      #object<- FindNeighbors(object, reduction = 'harmony', dims = 1:dims,k.param = 20)
      #set.seed(100)
      #object<- FindClusters(object,resolution = 0.4, verbose = FALSE)
      mclapply(c(15,18,20,23,25,28,30), function(dims) {
        paste('umap',sigma,theta,lambda,dims,sep = '_')->name
        object <- RunUMAP(object, reduction = 'harmony',dims = 1:dims,n.neighbors =n.neighbors,min.dist = min.dist,spread =spread,reduction.name = name)
        DimPlot(object,group.by = 'cell_type_l1',cols = colors,reduction = name)->p
        FeaturePlot(object,'IGHG1',reduction = name)
        #ggsave('p2.png',p,width=10,height=8,dpi=300)
        ggsave(paste0('/ix1/wchen/liutianhao/result/pathology_ST/script/test_umap3/',name,'_p.png'),p,width=10,height=8,dpi=300)
        DimPlot(object[,which(object$cell_type_l1=='B_Plasma')],group.by = 'cell_type_l1',cols = colors,reduction = name)->p
        ggsave(paste0('/ix1/wchen/liutianhao/result/pathology_ST/script/test_umap3/new_set/',name,'_p.png'),p,width=10,height=8,dpi=300)
      },mc.cores = 10)

    }
  }
}

object<- FindNeighbors(object, reduction = 'harmony', dims = 1:dims,k.param = 20,random.seed=100)
FindClusters(object = object,random.seed = 100,resolution = 1.5)->object
FindClusters(object = object,random.seed = 100,resolution = 1)->object
DimPlot(object,group.by = 'RNA_snn_res.1.5',reduction = name,label = T)
DimPlot(object,group.by = 'cell_type_l1',cols = colors,reduction = name,label = T)

length(which(object$seurat_clusters==19))
B_Plasma<-object
DimPlot(object,group.by = 'cell_type_l1',reduction = name,label = T,cols = colors)
"Endothelium"=c(3,20,21)
object$cell_type_l1_new<-'Keratinocyte'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(4,21,22))]<-'Endothelium'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(5,10))]<-'Smooth_muscle_cell'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(15,19))]<-'Glandular_epithelium'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(8))]<-'Myeloid'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(20))]<-'B_Plasma'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(12,2))]<-'Fibroblasts'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(26))]<-'Mast_cell'
object$cell_type_l1_new[which(object$seurat_clusters%in%c(7,13,29))]<-'T_cell'
object$cell_type_l1_new[which(object$RNA_snn_res.1.5%in%c(28))]<-'Neuron'
object$cell_type_l1_new[which(object$RNA_snn_res.1.5%in%c(29))]<-'Melanocyte'

DimPlot(object,group.by = 'cell_type_l1_new',reduction = name,label = F,cols = colors)->p
saveRDS(object,'/ix1/wchen/liutianhao/result/skin_ST/HD/singlecell/final_reference/plus_K_flex/reference_updated.rds')
ggsave('p.png',p,width = 11,height = 9,dpi = 300)

#for the spatial data UMAP
merged_filter<-readRDS('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds')
merged_filter$cell_type_l1_all[which(merged_filter$cell_type_l1_all%in%c('Macrophage','Neutrophil'))]<-'Myeloid'
merged_filter->object
merged_filter_marker <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/DE_final/single_cell_DE/merged_filter_marker.rds")
merged_filter_marker<-merged_filter_marker[which(merged_filter_marker$p_val_adj<0.05),]

NormalizeData(object)->object
FindVariableFeatures(object,nfeatures=2000)->object
VariableFeatures(object)<-union(VariableFeatures(object),unique(merged_filter_marker$gene))
object <- ScaleData(object,features = VariableFeatures(object))
set.seed(100)
object <- RunPCA(object,npcs=50)
ElbowPlot(object, ndims = 50,reduction = 'pca')
dims<-18
paste('umap',dims,sep = '_')->name
object <- RunUMAP(object, reduction = 'pca',dims = 1:dims,n.neighbors =20,min.dist = 0.18,spread =0.8,reduction.name = name)
DimPlot(object,group.by = 'cell_type_l1_all',cols = colors,reduction = name)

hard_cell<-c('Fibroblasts','T_cell','Myeloid','Endothelium','Smooth_muscle_cell')
object2<-object[,which(object$cell_type_l1_all%in%hard_cell)]
FindVariableFeatures(object2,nfeatures=1200)->object2
VariableFeatures(object2)<-union(VariableFeatures(object2),merged_filter_marker$gene[which(merged_filter_marker$cluster%in%c(hard_cell,'Macrophage','Neutrophil'))])
object2 <- ScaleData(object2,features = rownames(object2))
set.seed(100)
object2 <- RunPCA(object2,npcs=50)
ElbowPlot(object2, ndims = 50,reduction = 'pca')
dims<-30
paste('umap',dims,sep = '_')->name
object2 <- RunUMAP(object2, reduction = 'pca',dims = 1:dims,n.neighbors =20,min.dist = 0.10,spread =0.5,reduction.name = name)
DimPlot(object2,group.by = 'cell_type_l1_all',cols = colors,reduction = name)


-asd-asd
for (dims in c(20,15,25)) {
  for (n.neighbors in c(20,25,30,35,40)) {
    for (spread in c(1,3,5,7,9)) {
      for (min.dist in c(0.2,0.22,0.25,0.28,0.30,0.33)) {
        for (theta in c(1,1.5,2,2.5)) {
          for (lambda in c(1,2,3)) {
            sigma<-0.11
            set.seed(100)
            object <- RunHarmony(object,group.by.vars = 'orig.ident',lambda = lambda,theta = theta,sigma=sigma)
            #set.seed(100)
            #object<- FindNeighbors(object, reduction = 'harmony', dims = 1:dims,k.param = 20)
            #set.seed(100)
            #object<- FindClusters(object,resolution = 0.4, verbose = FALSE)
            paste('umap',dims,n.neighbors,spread,min.dist,theta,lambda,sep = '_')->name
            object <- RunUMAP(object, reduction = 'harmony',dims = 1:dims,n.neighbors =n.neighbors,min.dist = min.dist,spread =spread,reduction.name = name)
            DimPlot(object,group.by = 'cell_type_l1',cols = colors,reduction = name)->p
            FeaturePlot(object,'IGHG1',reduction = name)
            ggsave('p2.png',p,width=10,height=8,dpi=300)
            ggsave(paste0('/ix1/wchen/liutianhao/result/pathology_ST/script/test_umap2/',name,'_p.png'),p,width=10,height=8,dpi=300)
            DimPlot(object[,which(object$cell_type_l1=='B_Plasma')],group.by = 'cell_type_l1',cols = colors,reduction = name)->p
            ggsave(paste0('/ix1/wchen/liutianhao/result/pathology_ST/script/test_umap2/new_set/',name,'_p.png'),p,width=10,height=8,dpi=300)
            
          }
        }
      }
    }
  }
}



saveRDS(object,'/ix1/wchen/liutianhao/result/pathology_ST/script/test_umap/test_umap2.rds')
-as-as
set.seed(100)
object <- RunUMAP(object, reduction = 'harmony',dims = 1:20,n.neighbors =30,min.dist = 0.25,n.epochs = 700,negative.sample.rate = 10,spread = 1)
DimPlot(object,group.by = 'cell_type_l1',cols = colors,reduction = name)
DimPlot(object,group.by = 'seurat_clusters')
DimPlot(object,group.by = 'cell_type_l1',cols = colors)
FeaturePlot(object,'IGHG1',reduction =name)
FeaturePlot(object,'MS4A1',reduction =name)



