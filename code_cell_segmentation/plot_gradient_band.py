import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap
import numpy as np

def marker_size_pixels_to_s(side_pixels, dpi):
    """
    Convert desired square side in *pixels* (in the saved figure)
    to matplotlib scatter 's' (area in points^2).
    """
    points_per_pixel = 72.0 / dpi
    side_points = side_pixels * points_per_pixel
    return side_points ** 2
# df has columns: 'x', 'y', 'label'
# Convert labels to categorical to get integer codes

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

target = 'band'
for sample in list(sample2path.keys()):
  df = pd.read_csv(f'/ihome/wchen/tianhao/{sample}.csv')
  fig, ax = plt.subplots(figsize=(10, 10), dpi=600)
  dpi = 600          # or 300, etc.
  side_pixels = 8    # desired marker side in output PNG
  x = df['x'].to_numpy()
  y = df['y'].to_numpy()
  val = df[target].to_numpy()
  width_units  = x.max() - x.min()   # span of x in your units
  height_units = y.max() - y.min()   # span of y in your units
  fig_width_in  = width_units  / dpi
  fig_height_in = height_units / dpi

  s = marker_size_pixels_to_s(side_pixels, dpi)

  fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), dpi=dpi)

  sc = ax.scatter(
    x, y,
    c=val,
    cmap='viridis',
    marker='s',       # squares
    s=s,              # ~8x8 pixels in the PNG
    linewidths=0      # no boundary lines
  )

  ax.set_aspect('equal', 'box')
  ax.set_xlim(x.min(), x.max())
  ax.set_ylim(y.min(), y.max())

  ax.axis('off')
  cbar = fig.colorbar(sc, ax=ax)
  cbar.set_label(target)
  fig.savefig(f"{sample}_density.png", dpi=dpi, bbox_inches='tight', pad_inches=0)
  plt.close(fig)
