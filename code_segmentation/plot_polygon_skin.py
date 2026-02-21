import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from skimage import io
import colorsys
from matplotlib.colors import to_rgba
import colorsys

Image.MAX_IMAGE_PIXELS = None  # allow very large TIFFs



sample2path = {
    "A_2": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8691_A2.tif",
    "E_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/E1.tif",
    "E_2": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8693_E2.tif",
    "K_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/K1_1323.tif",
    "K_2": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/K2_1113.tif",
    "H_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/H1_40x.tif",
    "H_2": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/H2_40x.tif",
    "F_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/F1.tif",
    "I_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/I1_40x.tif",
    "J_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/J1_40x.tif",
    "G_2": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/4945_G2.tif",
    "M_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/M1_1045.tif",
    "G_1": "/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8594_G1.tif",
}

import colorsys


def brighten_hsv(ct2color,ct, v_increase=0.30, s_decrease=0.00):
    color = ct2color[ct]
    no_brighten = {"Melanocyte", "Smooth_muscle_cell"}
    if ct in no_brighten:
        return color
    else:
        r, g, b, *rest = color
        a = rest[0] if rest else 1.0
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        v = min(1.0, v + v_increase)
        s = max(0.0, s - s_decrease)
        r2, g2, b2 = colorsys.hsv_to_rgb(h, s, v)
        if rest:
            return (r2, g2, b2, a)
        else:
            return (r2, g2, b2)


def plot_cells(img,df,sample,no_pink_palette,out_path_overlay,black):
    cell_types = sorted(df['cell_type'].astype(str).unique())
    cell_types = df["cell_type"].cat.categories
    colors = [to_rgba(c) for c in no_pink_palette[:len(cell_types)]]
    ct2color = dict(zip(cell_types, colors))
    df = df[df.index.str.contains(sample)]
    img = Image.fromarray(img)
    W, H = img.size  # pixel dimensions
    legend_handles = [Patch(facecolor=brighten_hsv(ct2color,ct, v_increase=0), edgecolor='none', alpha=1, label=ct)
                  for ct in cell_types]
    fig_leg, ax_leg = plt.subplots(figsize=(4, len(cell_types) * 0.4))  # auto-size legend
    ax_leg.axis("off")
    ax_leg.set_facecolor("black")                  # axes background
    ax_leg.figure.patch.set_facecolor("black")     # figure background

    legend = ax_leg.legend(
        handles=legend_handles,
        loc='center left',
        frameon=True,
        fontsize=10
    )
    legend.get_frame().set_facecolor("black")
    legend.get_frame().set_edgecolor("black")      # remove border contrast if needed

    # Make legend text white so it's visible
    for text in legend.get_texts():
        text.set_color("white")

    legend.get_frame().set_alpha(1)
    out_path_legend = f"{out_path_overlay}/legend.tif"
    plt.savefig(out_path_legend, dpi=200, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig_leg)
    print(f"Saved legend separately to {out_path_legend}")

    polygon = np.load(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}/{sample}_polygons_expanded_array.npy')
    dpi = 100
    
    if black==False: 
        fig, ax = plt.subplots(figsize=(W / dpi, H / dpi), dpi=dpi)
        ax.set_facecolor('white')
        ax.imshow(img, origin='upper', alpha=0.8)  # <— lighten background by blending to white
        for idx, row in df.iterrows():
            i = int(idx.split("_")[-1])
            ct = str(row['cell_type'])
            pts = polygon[i].T       # (32, 2) = (x, y)
            pts = pts[:, ::-1]            # flip columns → (y, x)
            ax.fill(pts[:, 0], pts[:, 1],
                facecolor=ct2color[ct],
                edgecolor='none',
                alpha=0.8)
        
        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)   # y axis downward like image pixels
        ax.set_axis_off()
        plt.savefig(f'{out_path_overlay}/p_HE_color_{sample}', dpi=dpi, bbox_inches='tight', pad_inches=0,
            facecolor='white', transparent=False)
        plt.close(fig)
        
        fig, ax = plt.subplots(figsize=(W / dpi, H / dpi), dpi=dpi)
        ax.set_facecolor('white')
        ax.imshow(img, origin='upper', alpha=0.8)  # <— lighten background by blending to whit
        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)   # y axis downward like image pixels
        ax.set_axis_off()
        plt.savefig(f'{out_path_overlay}/p_blank_{sample}', dpi=dpi, bbox_inches='tight', pad_inches=0,
            facecolor='white', transparent=False)
        plt.close(fig)
    if black==True:
        fig, ax = plt.subplots(figsize=(W / dpi, H / dpi), dpi=dpi)
        ax.set_facecolor('black')  # ensure white background
        for idx, row in df.iterrows():
            i = int(idx.split("_")[-1])
            ct = str(row['cell_type'])
            pts = polygon[i].T       # (32, 2) = (x, y)
            pts = pts[:, ::-1]       # flip columns → (y, x)
            bright_color = brighten_hsv(ct2color,ct, v_increase=0)
            ax.fill(
                pts[:, 0], pts[:, 1],
                facecolor=bright_color,
                edgecolor='none',
                alpha=1
            )

        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)   # y axis downward
        ax.set_axis_off()
        plt.savefig(f'{out_path_overlay}/p_black_{sample}',dpi=dpi,bbox_inches='tight',pad_inches=0,facecolor='black', transparent=False)
        plt.close(fig)


palette_all =["#A6CEE3","#D95F02","#00BBFF","#7570B3","#8DD3C7","#66A61E","#B2DF8A","#E6AB02","#E7298A","#A6761D"]
#palette_all =["#A6CEE3","#EC940F","#00BBFF","#7570B3","#8DD3C7","#66A61E","#B2DF8A","#F9BD07","#E7298A","#A6761D"]
#plot immune proximal fibroblasts
palette_inf_fib = ["#A6CEE3","#E6AB02","#CAB2D6","#66A61E","#CC00CC"]
for target in ['infinite','220','132']:
  for sample in list(sample2path.keys()):
    img = io.imread(sample2path[sample])
    df=pd.read_csv(f'~/immune_fib/result_all_{target}.csv',index_col=0)
    plot_cells(img,df,sample,palette_inf_fib,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_fib/{target}',True)

palette_inf_fib = ["#A6CEE3","#E6AB02","#CAB2D6","#66A61E","#CC00CC"]
for target in ['infinite','220','132']:
  for sample in list(sample2path.keys()):
    img = io.imread(sample2path[sample])
    df=pd.read_csv(f'~/immune_fib/result_all_{target}.csv',index_col=0)
    plot_cells(img,df,sample,palette_inf_fib,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_fib/{target}/H_E',False)

palette_subtype_fib = ["#984EA3","#FFFF00","#00BBFF","#1B9E77"]
sample = 'E_1'
img = io.imread(sample2path[sample])
df=pd.read_csv('~/macrophage_fibroblast_colo/Fibroblast_subtype_E_1.csv',index_col=0)
df['cell_type'] = pd.Categorical(df['cell_type'],categories=["ANGPTL7+", "SFRP2+/SFRP4+/PCOLCE2+","STC1+/FZD9+","Other"])
plot_cells(img,df,sample,palette_subtype_fib,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_fib/fib_subtype_E1_2/',True)
plot_cells(img,df,sample,palette_subtype_fib,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_fib/fib_subtype_E1_2/',False)

-asd-asd-asd
#plot all cells with H&E background and blank H&E
for sample in list(sample2path.keys()):
  img = io.imread(sample2path[sample])
  df=pd.read_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/result_{sample}_2.csv',index_col=1)
  plot_cells(img,df,sample,palette_all,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/HE_color',False)

-asd-asd-asd
#plot all cells in black background
for sample in list(sample2path.keys()):
  img = io.imread(sample2path[sample])
  df=pd.read_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/result_{sample}_2.csv',index_col=1)
  plot_cells(img,df,sample,palette_all,'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black',True)


#plot immune proximal endothelium
palette_inf_endo = ["#A6CEE3","#E6AB02","#CAB2D6","#66A61E","#E31A1C"]
for target in ['infinite','220','132']:
  for sample in list(sample2path.keys()):
    img = io.imread(sample2path[sample])
    df=pd.read_csv(f'~/immune_endo/result_all_{target}.csv',index_col=0)
    plot_cells(img,df,sample,palette_inf_endo,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_endo/{target}',True)


palette_inf_endo = ["#A6CEE3","#E6AB02","#CAB2D6","#66A61E","#E31A1C"]
for target in ['infinite','220','132']:
  for sample in list(sample2path.keys()):
    img = io.imread(sample2path[sample])
    df=pd.read_csv(f'~/immune_endo/result_all_{target}.csv',index_col=0)
    plot_cells(img,df,sample,palette_inf_endo,f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/black/inf_endo/{target}/H_E',False)

#plor endo fib and immune
palette_all = ["#A6CEE3","#E31A1C","#CC00CC","#E6AB02","#CAB2D6","#66A61E"]
for sample in list(sample2path.keys()):
  img = io.imread(sample2path[sample])
  df=pd.read_csv(f'/ix1/wchen/liutianhao/result/pathology_ST/deconvolution/single_cell_level/l1/result_{sample}_2.csv',index_col=1)
  df=pd.read_csv('a.csv',index_col=1)
  df = df[df['cell_type'].isin(['B_Plasma','T_cell','Macrophage','Neutrophil','Fibroblasts','Endothelium'])]
  plot_cells(img,df,sample,palette_all,'/ihome/wchen/tianhao',True)


no_pink_palette =["#A6CEE3","#D95F02","#1B9E77","#7570B3","#8DD3C7","#66A61E","#B2DF8A","#E6AB02","#E7298A","#A6761D"]
no_pink_palette =["#A6CEE3","#D95F02","#7D0CEA","#7570B3","#8DD3C7","#66A61E","#B2DF8A","#E6AB02","#E7298A","#A6761D"]
no_pink_palette =['#D95F02', '#B2DF8A', '#E7298A','tab:blue','tab:red']
no_pink_palette =['#D95F02', '#B2DF8A', '#E7298A','tab:blue','tab:purple']
no_pink_palette =["#D95F02","#1B9E77"]
colors = [to_rgba(c) for c in no_pink_palette[:len(cell_types)]]
ct2color = dict(zip(cell_types, colors))
polygon = np.load(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}/{sample}_polygons_expanded_array.npy')

dpi = 100
fig_w_in = W / dpi
fig_h_in = H / dpi
fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in), dpi=dpi)


#plot polygons on HE
ax.set_facecolor('white')
ax.imshow(img, origin='upper', alpha=0.8)  # <— lighten background by blending to white

for idx, row in df.iterrows():
    i = int(idx.split("_")[-1])
    ct = str(row['cell_type'])
    pts = polygon[i].T       # (32, 2) = (x, y)
    pts = pts[:, ::-1]            # flip columns → (y, x)

    ax.fill(pts[:, 0], pts[:, 1],
            facecolor=ct2color[ct],
            edgecolor='none',
            alpha=0.8)

ax.set_xlim(0, W)
ax.set_ylim(H, 0)   # y axis downward like image pixels
ax.set_axis_off()

out_path_overlay = f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}/{sample}_polygon_color_inf_endo.tif"
plt.savefig(out_path_overlay, dpi=dpi, bbox_inches='tight', pad_inches=0,
            facecolor='white', transparent=False)
plt.close(fig)

#plot legend


#plot polygons on white background
img = io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/E1.tif')
img = io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/F1.tif')
img = io.imread("/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/M1_1045.tif")
img = io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/I1_40x.tif')
img = io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8594_G1.tif')
img = io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/4945_G2.tif')
img = Image.fromarray(img)
W, H = img.size  # pixel dimensions

dpi = 100
fig_w_in = W / dpi
fig_h_in = H / dpi
fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in), dpi=dpi)


fig, ax = plt.subplots(figsize=(W/dpi, H/dpi), dpi=dpi)
ax.set_facecolor('black')  # ensure white background

for idx, row in df.iterrows():
    i = int(idx.split("_")[-1])
    ct = str(row['cell_type'])
    pts = polygon[i].T       # (32, 2) = (x, y)
    pts = pts[:, ::-1]       # flip columns → (y, x)

    ax.fill(
        pts[:, 0], pts[:, 1],
        facecolor=ct2color[ct],
        edgecolor='none',
        alpha=0.8
    )

ax.set_xlim(0, W)
ax.set_ylim(H, 0)   # y axis downward
ax.set_axis_off()

out_path_overlay = f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}/{sample}_polygon_color_inf_fib_{color_back}.tif"
plt.savefig(
    out_path_overlay,
    dpi=dpi,
    bbox_inches='tight',
    pad_inches=0,
    facecolor=color_back,   # keep white background in output
    transparent=False
)
plt.close(fig)




