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
import colorsys
work_path = '/ihome/wchen/tianhao/'

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
    cell_types = df["cell_type"].astype("category").cat.categories
    colors = [to_rgba(c) for c in no_pink_palette[:len(cell_types)]]
    ct2color = dict(zip(cell_types, colors))
    df = df[df.index.str.contains(sample[:3])]
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

    polygon = np.load(f'{work_path}/Aggregation_result_visualization/data/{sample}_polygons_expanded_array.npy')
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


palette_all = ["#A6CEE3","#D95F02","#00BBFF","#7570B3","#8DD3C7","#E6AB02","#B2DF8A","#E7298A","#CAB2D6","#A6761D","#66A61E"]
img = io.imread(f'{work_path}/Aggregation_result_visualization/data/M1_1045.tif')
df=pd.read_csv(f'{work_path}/Aggregation_result_visualization/data/M_1_celltype.csv',index_col=0)
plot_cells(img,df,'M_1',palette_all,f'{work_path}/Aggregation_result_visualization/result/',True)


