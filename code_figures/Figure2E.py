import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.colors import to_rgba
from skimage import io

def overlay_cells_union_per_type_allow_cross_type_blend(
    he_img,
    cells_df,
    out_png,
    x_col="x",
    y_col="y",
    label_col="cell_type",
    radius_um=33,
    microns_per_pixel=None,
    radius_px=None,
    celltype_colors=None,
    radius_alpha=0.25,
    dot_size=3,
    dot_color="black",
    dpi=100,
    show=False,
    paint_order=None,  # list of cell types to control draw order (optional)
):
    """
    - Same-type overlap: NO extra darkening (union mask per type)
    - Different-type overlap: allowed (alpha blending between layers)
    """

    # --- radius conversion ---
    if radius_px is None:
        if microns_per_pixel is None:
            raise ValueError("Provide microns_per_pixel or radius_px")
        radius_px = radius_um / microns_per_pixel
    radius_px = float(radius_px)

    H, W = he_img.shape[:2]

    df = cells_df[[x_col, y_col, label_col]].dropna().copy()
    df[x_col] = df[x_col].astype(float)
    df[y_col] = df[y_col].astype(float)

    # keep points inside image
    df = df[
        (df[x_col] >= 0) & (df[x_col] < W) &
        (df[y_col] >= 0) & (df[y_col] < H)
    ].reset_index(drop=True)

    cell_types = sorted(df[label_col].unique())

    # color map only for radius layers
    if celltype_colors is None:
        cmap = plt.get_cmap("tab20")
        celltype_colors = {ct: cmap(i % cmap.N) for i, ct in enumerate(cell_types)}

    # choose draw order (matters for how blending looks)
    if paint_order is None:
        paint_order = cell_types
    else:
        # keep only those present; append missing ones
        paint_order = [ct for ct in paint_order if ct in cell_types] + [ct for ct in cell_types if ct not in paint_order]

    # Precompute a disk footprint (integer offsets)
    r = int(np.ceil(radius_px))
    yy, xx = np.ogrid[-r:r+1, -r:r+1]
    disk = (xx*xx + yy*yy) <= radius_px*radius_px  # boolean (2r+1, 2r+1)

    # --- figure with exact pixel size ---
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(he_img, interpolation="nearest")
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    # --- draw one RGBA layer per type (union mask -> constant alpha) ---
    for ct in paint_order:
        sub = df[df[label_col] == ct]
        if sub.empty:
            continue

        mask = np.zeros((H, W), dtype=bool)

        xs = np.rint(sub[x_col].to_numpy()).astype(int)
        ys = np.rint(sub[y_col].to_numpy()).astype(int)

        for x0, y0 in zip(xs, ys):
            x1, x2 = max(0, x0 - r), min(W, x0 + r + 1)
            y1, y2 = max(0, y0 - r), min(H, y0 + r + 1)

            dx1, dx2 = x1 - (x0 - r), x2 - (x0 - r)
            dy1, dy2 = y1 - (y0 - r), y2 - (y0 - r)

            mask[y1:y2, x1:x2] |= disk[dy1:dy2, dx1:dx2]

        # Build an RGBA image layer for this type (constant alpha on union mask)
        layer = np.zeros((H, W, 4), dtype=np.float32)
        rgba = np.array(to_rgba(celltype_colors[ct], radius_alpha), dtype=np.float32)
        layer[mask, :] = rgba

        # Plot this layer ONCE -> same-type overlap won't darken; cross-type overlaps will blend
        ax.imshow(layer, interpolation="nearest")

    # --- center dots (single color) ---
    ax.scatter(df[x_col], df[y_col], s=dot_size, c=dot_color, linewidths=0, zorder=10)

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", pad_inches=0)
    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, ax


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

def lighten_he_image(img, factor=0.3):
    """
    Lighten an RGB H&E image by blending with white.
    factor ∈ [0,1]: 0 = original, 1 = white
    """
    img = img.astype(np.float32)
    white = np.full_like(img, 255.0)
    out = img * (1 - factor) + white * factor
    return np.clip(out, 0, 255).astype(np.uint8)

for sample in list(sample2path.keys()):
    he = io.imread(sample2path[sample])
    df = pd.read_csv("/ihome/wchen/tianhao/inf_fib.csv",index_col=0)  
    df = df[df.index.str.contains(sample)]
    df.columns = ['y','x','cell_type']
    df = df[df["cell_type"] == "inf_fib"]
    celltype_colors = {
        "Fibroblasts": "#377eb8",
        "inf_fib": "#F3B02E"
    }
    he = lighten_he_image(he, factor=0.25)
    overlay_cells_union_per_type_allow_cross_type_blend(
      he_img= he,
      cells_df=df,
      out_png=f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/radius_no_overlap/inf_only/{sample}_fib_0.png",
      microns_per_pixel=0.25,
      radius_um=33,
      celltype_colors=celltype_colors,
      radius_alpha=0.6,
      paint_order=["inf_fib","Fibroblasts"],
      dot_size=0)
      

for sample in list(sample2path.keys()):
    he = io.imread(sample2path[sample])
    df = pd.read_csv("/ihome/wchen/tianhao/inf_endo.csv",index_col=0)  
    df = df[df.index.str.contains(sample)]
    df.columns = ['y','x','cell_type']
    df = df[df["cell_type"] == "inf_endo"]
    celltype_colors = {
        "Endothelium": "#377eb8",
        "inf_endo": "#F3B02E"
    }
    he = lighten_he_image(he, factor=0.25)
    overlay_cells_union_per_type_allow_cross_type_blend(
      he_img=he,
      cells_df=df,
      out_png=f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/radius_no_overlap/inf_only/{sample}_endo_0.png",
      microns_per_pixel=0.25,
      radius_um=33,
      celltype_colors=celltype_colors,
      radius_alpha=0.6,
      paint_order=["inf_endo","Endothelium"],
      dot_size=0)



for sample in list(sample2path.keys()):
    he = io.imread(sample2path[sample])
    df = pd.read_csv("/ihome/wchen/tianhao/inf_fib.csv",index_col=0)  
    df = df[df.index.str.contains(sample)]
    df.columns = ['y','x','cell_type']
    celltype_colors = {
        "Fibroblasts": "#377eb8",
        "inf_fib": "#F3B02E"
    }
    he = lighten_he_image(he, factor=0.25)
    overlay_cells_union_per_type_allow_cross_type_blend(
      he_img=he,
      cells_df=df,
      out_png=f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/radius_no_overlap/{sample}_fib_0.png",
      microns_per_pixel=0.25,
      radius_um=33,
      celltype_colors=celltype_colors,
      radius_alpha=0.6,
      paint_order=["inf_fib","Fibroblasts"],
      dot_size=0)
      

for sample in list(sample2path.keys()):
    he = io.imread(sample2path[sample])
    df = pd.read_csv("/ihome/wchen/tianhao/inf_endo.csv",index_col=0)  
    df = df[df.index.str.contains(sample)]
    df.columns = ['y','x','cell_type']
    celltype_colors = {
        "Endothelium": "#377eb8",
        "inf_endo": "#F3B02E"
    }
    he = lighten_he_image(he, factor=0.25)
    overlay_cells_union_per_type_allow_cross_type_blend(
      he_img=he,
      cells_df=df,
      out_png=f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/image_out/radius_no_overlap/{sample}_endo_0.png",
      microns_per_pixel=0.25,
      radius_um=33,
      celltype_colors=celltype_colors,
      radius_alpha=0.6,
      paint_order=["inf_endo","Endothelium"],
      dot_size=0)

