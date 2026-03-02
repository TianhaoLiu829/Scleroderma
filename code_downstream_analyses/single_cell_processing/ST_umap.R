
merged_filter <- readRDS("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/merged_filter.rds")
merged_filter$cell_type_l1_all[which(merged_filter$cell_type_l1_all%in%c('Neutrophil','Macrophage'))]<-'Myeloid'
DimPlot(merged_filter2,cols = colors,group.by = 'cell_type_l1_all')
sc_pro(merged_filter,T)->merged_filter2


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
  "T_cell"="#66A61E")

merged_filter->object
NormalizeData(object)->object
FindVariableFeatures(object,nfeatures=2500)->object
object <- ScaleData(object)
set.seed(100)
object <- RunPCA(object,npcs=50)
ElbowPlot(object, ndims = 50)
dims<-23
#object<- FindNeighbors(object, reduction = 'pca', dims = 1:dims,k.param = 20)
#object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
object <- RunUMAP(object, reduction = 'pca',dims = 1:dims,n.neighbors =20,min.dist=0.2,spread = 1)
cell_types <- c("Melanocyte","Neuron","B_Plasma","T_cell","Smooth_muscle_cell","Myeloid","Endothelium","Glandular_epithelium","Fibroblasts","Keratinocyte")
object$cell_type_l1_all<-factor(object$cell_type_l1_all,levels = cell_types[seq(length(cell_types),1,-1)])
DimPlot(object,group.by = 'cell_type_l1_all',cols = colors)->p
ggsave('p.png',p,width = 11,height = 9,dpi = 300)
object->object3

sc_pro<-function(object,harmony){
  NormalizeData(object)->object
  FindVariableFeatures(object,nfeatures=2500)->object
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object,npcs=50)
  ElbowPlot(object, ndims = 50)
  if (harmony==FALSE) {
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'pca', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30,n.neighbors =15,min.dist=0.2)
  } else {
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = 'orig.ident',theta = 2,lambda =1,sigma = 0.1)
    set.seed(100)
    object<- FindNeighbors(object, reduction = 'harmony', dims = 1:30,k.param = 25)
    set.seed(100)
    object<- FindClusters(object,resolution = 0.7, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30,n.neighbors =15,min.dist=0.2)
  }
  return(object)
}
