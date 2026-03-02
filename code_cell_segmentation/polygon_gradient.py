import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point
from shapely.strtree import STRtree
from shapely.ops import unary_union
from collections import OrderedDict

def band_areas_for_polygon(poly, edges):
    areas = OrderedDict()
    areas["inside"] = poly.area
    for lo, hi in zip(edges[:-1], edges[1:]):
        outer = poly.buffer(hi)      # polygon grown by hi pixels
        inner = poly.buffer(lo)      # polygon grown by lo pixels
        ring = outer.difference(inner)
        areas[f"{int(lo)}-{int(hi)}"] = ring.area

    return areas


def classify_band(distance, inside, edges):
    if inside:
        return "inside"
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo <= distance < hi:
            return f"{lo}-{hi}"
    return f">={edges[-1]}"


# --- load your polygons ---
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

for sample in list(sample2path.keys()):
  seg = np.load(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}_2/{sample}_2_polygons_nuclear.npy',allow_pickle=True)

  poly_arrays = []
  for poly in seg:
      coords = np.array(poly.exterior.coords)  # shape (n_points+1, 2)
      coords = coords[:-1]                     # shape (n_points, 2)
      coords = coords.T                        # shape (2, n_points)
      poly_arrays.append(coords)
  seg = np.stack(poly_arrays, axis=0)   # shape (n_polygons, 2, n_points)

  polygons = []
  for i in range(seg.shape[0]):
      xs = seg[i, 0, :]
      ys = seg[i, 1, :]
      coords = list(zip(xs, ys))     # [(x1,y1), (x2,y2), ...]
      poly = Polygon(coords)
      if not poly.is_valid:
          poly = poly.buffer(0)      # fix self-intersections if any
      polygons.append(poly)

  tree = STRtree(polygons)
  geom_to_idx = {id(g): i for i, g in enumerate(polygons)}

  max_dist = 80
  edges = np.arange(0, max_dist + 4, 4)  # [0,1,2,...,10]
  
  all_poly = unary_union(polygons)           # union of all polygons
  global_band_areas = band_areas_for_polygon(all_poly, edges)
  global_areas_df = pd.DataFrame([global_band_areas])
  px_size_um = 0.25
  area_cols = ["inside"] + [f"{int(lo)}-{int(hi)}" for lo, hi in zip(edges[:-1], edges[1:])]
  global_areas_um2 = global_areas_df.copy()
  global_areas_um2[area_cols] = global_areas_df[area_cols] * (px_size_um ** 2)
  global_areas_um2.to_csv(f"global_band_areas_um_{sample}_80.csv", index=False)
  
  spots_df = pd.read_csv(f"/ix1/wchen/liutianhao/result/pathology_ST/position_filter/{sample}_filtered_positions.csv")
  spots_df = spots_df[spots_df['in_tissue']==True]
  spots_df["point"] = spots_df.apply(lambda r: Point(r["pxl_row_in_fullres"], r["pxl_col_in_fullres"]), axis=1)

  results = []
  for _, row in spots_df.iterrows():
      pt = row["point"]

    # In your environment this returns an integer index:
      nearest_idx = tree.nearest(pt)          # numpy.int64
      poly_idx = int(nearest_idx)             # convert to plain int
      poly = polygons[poly_idx]               # get actual Polygon

      inside = poly.contains(pt) or poly.touches(pt)
      d = poly.distance(pt)
      band = classify_band(d, inside, edges)

      results.append({
        "barcode": row["barcode"],
        "x": row["pxl_row_in_fullres"],
        "y": row["pxl_col_in_fullres"],
        "polygon_id": poly_idx,
        "distance_px": d,
        "band": band,
      })

  spots_result = pd.DataFrame(results)
  #spots_result = spots_result[spots_result['band']!='>=16']
  spots_result.to_csv(f"/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/spots_gradient/spots_gradient_{sample}_80.csv", index=False)

-a-sd-asd-asd
for sample in list(sample2path.keys()):
  seg = np.load(f'/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/{sample}_2/{sample}_2_polygons_nuclear.npy',allow_pickle=True)

  poly_arrays = []
  for poly in seg:
      coords = np.array(poly.exterior.coords)  # shape (n_points+1, 2)
      coords = coords[:-1]                     # shape (n_points, 2)
      coords = coords.T                        # shape (2, n_points)
      poly_arrays.append(coords)
  seg = np.stack(poly_arrays, axis=0)   # shape (n_polygons, 2, n_points)

  polygons = []
  for i in range(seg.shape[0]):
      xs = seg[i, 0, :]
      ys = seg[i, 1, :]
      coords = list(zip(xs, ys))     # [(x1,y1), (x2,y2), ...]
      poly = Polygon(coords)
      if not poly.is_valid:
          poly = poly.buffer(0)      # fix self-intersections if any
      polygons.append(poly)

  tree = STRtree(polygons)
  geom_to_idx = {id(g): i for i, g in enumerate(polygons)}

  max_dist = 80
  edges = np.arange(0, max_dist + 4, 4)  # [0,1,2,...,10]
