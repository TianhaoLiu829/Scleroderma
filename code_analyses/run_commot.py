import numpy as np
import commot as ct
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os
from scipy import sparse
from scipy.io import mmwrite



adata = sc.read_h5ad("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter.h5ad")
coords = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/merged_filter_loc.csv", index_col=0)
adata2 = adata[adata.obs_names==coords.index,:]
coords2 = coords.iloc[adata.obs_names == coords.index,:]
adata2.var_names_make_unique()
adata_full = AnnData(
    X=adata2.raw.X.copy(),
    obs=adata2.obs.copy(),
    var=adata2.raw.var.copy()
)
adata_full.obsm["spatial"] = coords2.loc[adata_full.obs_names].values

def preprocess(object):
    object.raw = object
    sc.pp.normalize_total(object, inplace=True)
    sc.pp.log1p(object)
    return object
adata_full2 = preprocess(adata_full)


df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type=None, database='CellChat')
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_full2, min_cell_pct=0.01)
print(df_cellchat_filtered.shape)
df_cellchat_filtered.columns = df_cellchat.columns



#exclude epidermis
all_adata = []
epidermis = pd.read_csv('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/epidermis.csv',index_col=0)
adata_full2_filter = adata_full2[epidermis.index,:]

adata_full2_filter.obs['epidermis'] = epidermis.iloc[:,0]
adata_full2_filter_epi = adata_full2_filter[adata_full2_filter.obs['epidermis']==False,:]
all_adata_filter = []
for sample in np.unique(adata_full2_filter_epi.obs['orig.ident'].values):
    target = adata_full2_filter_epi[adata_full2_filter_epi.obs['orig.ident']==sample,:]
    ct.tl.spatial_communication(target,database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
    all_adata.append(target)
    target.write_h5ad(f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/{sample}_commot_filter.h5ad")

#test significance
for pathway in adata.uns['commot-cellchat-info']['df_ligrec']['ligand']+'-'+adata.uns['commot-cellchat-info']['df_ligrec']['receptor']:
  for sample in np.unique(adata_full2.obs['orig.ident'].values):
    adata = sc.read_h5ad(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/{sample}_commot_filter.h5ad')
    ct.tl.communication_direction(adata, database_name='cellchat', pathway_name=pathway, k=5)
    ct.tl.cluster_communication(adata, database_name='cellchat', pathway_name=pathway, clustering='cell_type_l1_all',n_permutations=100,random_seed=100)
    adata.uns[f'commot_cluster-cell_type_l1_all-cellchat-{pathway}']['communication_pvalue'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/{sample}_{pathway}_p.csv')
    adata.uns[f'commot_cluster-cell_type_l1_all-cellchat-{pathway}']['communication_matrix'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_epi/{sample}_{pathway}_effect.csv')

#all data
all_adata = []
for sample in np.unique(adata_full2.obs['orig.ident'].values):
    target = adata_full2[adata_full2.obs['orig.ident']==sample,:]
    ct.tl.spatial_communication(target,database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
    all_adata.append(target)
    target.write_h5ad(f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_all/{sample}_commot.h5ad")


sdsdsd

#near fib regions
all_adata = []
min_dis = pd.read_csv('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/distance_all.csv',index_col=0)
adata_full2_filter = adata_full2[min_dis.index,:]
for center in min_dis.columns:
    adata_full2_filter.obs[center] = min_dis[center].values
    adata_full2_filter_near = adata_full2_filter[adata_full2_filter.obs[center]<80,:]
    all_adata_filter = []
    for sample in np.unique(adata_full2_filter_near.obs['orig.ident'].values):
        target = adata_full2_filter_near[adata_full2_filter_near.obs['orig.ident']==sample,:]
        ct.tl.spatial_communication(target,database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
        all_adata.append(target)
        center_name = center.replace("min_dis_to_", "")
        target.write_h5ad(f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_{center_name}/{sample}_commot_filter.h5ad")



#get the cell chat result and read into R
import glob
for target in ['all','epi','endo','fib','kera','T','Macrophage'] :
  files_all = glob.glob(f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/filter_{target}/*.h5ad")
  cell_chat_all = []
  for file in files_all:
      data = sc.read_h5ad(file)
      cell_chat_all.append(data)
  if target=='all':
    dir_list = [s.replace("_commot.h5ad", "") for s in files_all]
  else:
    dir_list = [s.replace("_commot_filter.h5ad", "") for s in files_all]
    
  for i in range(len(files_all)):
      filename = files_all[i].split("/")[-1]      # "H_2_commot.h5ad"
      filename = filename.split("_commot")[0]   # "H_2"
      cell_chat_all[i].obsm['commot-cellchat-sum-receiver'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/result_filter_{target}/cell_chat_'+filename+'_recei.csv')
      cell_chat_all[i].obsm['commot-cellchat-sum-sender'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/result_filter_{target}/cell_chat_'+filename+'_send.csv')
      outdir = dir_list[i]
      os.makedirs(outdir, exist_ok=True)
      pairwise = cell_chat_all[i].obsp 
      pd.Series(cell_chat_all[i].obs_names).to_csv(f"{outdir}/nodes.tsv", index=False, header=False)
      for key, mat in pairwise.items():
          coo = mat.tocoo()  # mmwrite expects COO
          mmwrite(f"{outdir}/{key}.mtx", coo)

  for i in range(len(files_all)):
      filename = files_all[i].split("/")[-1]      # "H_2_commot.h5ad"
      filename = filename.split("_commot")[0]   # "H_2"
      cell_chat_all[i].obsm['commot-cellchat-sum-receiver'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/result_filter_{target}/cell_chat_'+filename+'_recei.csv')
      cell_chat_all[i].obsm['commot-cellchat-sum-sender'].to_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/single_cell/MT_filter/cell_chat/result_filter_{target}/cell_chat_'+filename+'_send.csv')





