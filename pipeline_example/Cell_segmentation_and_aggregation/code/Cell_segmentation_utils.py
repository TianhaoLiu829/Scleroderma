import matplotlib.pyplot as plt
from skimage.measure import regionprops
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
from stardist import random_label_cmap
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from csbdeep.utils import normalize
from skimage.draw import circle_perimeter
import pandas as pd
from PIL import Image, ImageDraw
from shapely.geometry import Polygon
from shapely import affinity
from scipy.spatial import KDTree
from tqdm import tqdm
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
import pyarrow.parquet as pa
from scipy.sparse import csr_matrix

import scanpy as sc
from collections import defaultdict
from scipy.io import mmwrite


#required functions
def convert_to_polygons(polygon_coords):
    polygons = []
    for i in range(polygon_coords.shape[0]):
        x = polygon_coords[i, 0, :]
        y = polygon_coords[i, 1, :]
        coords = list(zip(x, y))
        polygons.append(Polygon(coords))
    return polygons

def density_to_expansion(d, max_density,max_exp=20, min_exp=0.5):
    return np.clip(max_exp * (1 - d / max_density), min_exp, max_exp)


def expand_polygons_array(polygons, distances):
    resampled_polygons = []
    n_points = 32
    for poly, dist in zip(polygons, distances):
        poly_buffered = poly.buffer(dist)
        boundary = poly_buffered.exterior
        total_len = boundary.length
        distances2 = np.linspace(0, total_len, n_points, endpoint=False)
        sampled_points = [boundary.interpolate(d) for d in distances2]
        coords = np.array([[pt.x, pt.y] for pt in sampled_points]).T
        resampled_polygons.append(coords)
    return resampled_polygons

def expand_polygons(polygons, distances):
    return [poly.buffer(dist) for poly, dist in zip(polygons, distances)]

def draw_polygons_on_image(image, polygons, color=(255, 0, 0), width=1, fill=False):
    draw = ImageDraw.Draw(image, "RGBA")

    for poly in polygons:
        coords = [(x, y) for y, x in poly.exterior.coords]
        if fill:
            draw.polygon(coords, fill=color + (100,), outline=color + (255,))
        else:
            draw.line(coords + [coords[0]], fill=color + (255,), width=width)  # close the loop

    return image

def reduce_overlaps(expanded_polygons, original_polygons, expansion_distances, step=1.0, min_expansion=0.0, max_iter=10):

    n = len(expanded_polygons)
    for _ in range(max_iter):
        # Build spatial index
        tree = STRtree(expanded_polygons)
        overlap_indices = set()

        # Detect overlaps
        for i, poly in enumerate(expanded_polygons):
            candidates = tree.query(poly)
            for j in candidates:
                if i >= j:
                    continue  # avoid double checking
                if poly.intersects(expanded_polygons[j]):
                    overlap_indices.update([i, j])
        
        if not overlap_indices:
            break  # no more overlaps

        # Reduce expansion for overlapping polygons
        for idx in overlap_indices:
            if expansion_distances[idx] > min_expansion:
                expansion_distances[idx] = max(min_expansion, expansion_distances[idx] - step)
                expanded_polygons[idx] = original_polygons[idx].buffer(expansion_distances[idx])
                
    
    return expanded_polygons, expansion_distances, overlap_indices

def sample_polygon_boundaries(polygons, n_points=32):
    sampled_coords = []
    for poly in polygons:
        boundary = poly.exterior
        total_length = boundary.length
        distances = np.linspace(0, total_length, n_points, endpoint=False)
        points = [boundary.interpolate(d) for d in distances]
        coords = np.array([[pt.x, pt.y] for pt in points]).T  # shape (2, 32)
        sampled_coords.append(coords)
    return np.stack(sampled_coords)  # shape (N, 2, 32)


def combine_draw_expand_step1(image,sample,spot_path,outpath,distance_radius):
    #run star-dist
    model = StarDist2D.from_pretrained('2D_versatile_he')
    img_normalized = normalize(image, 1, 99.8, axis=(0,1))
    labels, details = model.predict_instances(img_normalized)

    #get polygon
    polygons = details['coord']

    N, _, P = polygons.shape
    centroids = polygons.mean(axis=2)

    #estimate density
    tree = KDTree(centroids)
    radius_px = distance_radius  # e.g., radius to define local neighborhood
    neighbors = tree.query_ball_point(centroids, r=radius_px)
    local_density = np.array([len(n) for n in neighbors])
    maxdis = local_density.max()

    #expansion of polygons

    expansions = np.array([density_to_expansion(d,maxdis,max_exp=16, min_exp=0.5) for d in local_density])
    

    polygons_c = convert_to_polygons(polygons)
    polygons_expanded = expand_polygons(polygons_c,expansions)
    polygons_expanded_array = expand_polygons_array(polygons_c,expansions)
    polygons_expanded_array = np.stack(polygons_expanded_array)
    
    #restrict the expansion of polygons
    updated_polygon, updated_expansions, overlap = reduce_overlaps(polygons_expanded,polygons_c,expansions,step=0.5, min_expansion=2.0, max_iter=20)
    
    polygons_expanded_array = sample_polygon_boundaries(updated_polygon)  # shape (N, 2, 32)
    
    centroid_coords = np.array([(poly.centroid.x, poly.centroid.y) for poly in polygons_expanded])
    centroid_df = pd.DataFrame(centroid_coords, columns=['x', 'y'])
    centroid_df.to_csv(outpath+'/'+sample+"_center.csv",index=True)
    
    np.save(outpath+'/'+sample+'_polygons_expanded_array.npy', polygons_expanded_array)
    np.save(outpath+'/'+sample+'_expandions.npy',updated_expansions)
    np.save(outpath+'/'+sample+'_overlap.npy',np.array(overlap))
    #Identify and plot squares
    spots=pd.read_csv(spot_path)
    #spots = spots.to_pandas()
    spots.index = spots.loc[:,'barcode']
    spots = spots.loc[:,['pxl_row_in_fullres', 'pxl_col_in_fullres']]
    from matplotlib.path import Path
    spot_to_polygon = -np.ones(len(spots), dtype=int)
    
    for i, poly in enumerate(polygons_expanded_array.transpose(0, 2, 1)):
        path = Path(poly)
        inside = path.contains_points(spots)
        spot_to_polygon[inside] = i  # this assigns the polygon ID to each spot inside

    spots['polygon_id'] = spot_to_polygon
    spots[spots["polygon_id"] >= 0].to_csv(outpath+"/"+sample+"_spot_to_polygon.csv", index=True)

    # Draw polygons on the image
    image_2 = Image.fromarray(image)
    img_with_polygons = draw_polygons_on_image(image_2, polygons_expanded, color=(0, 255, 127), width=4, fill=False)

    img_with_polygons.save(outpath+'/'+sample+'_polygon.tif', format="TIFF")
    
    # Draw squares on the image
    #spots = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/"+sample+"/"+sample+"_spot_to_polygon.csv")
    spots = spots[spots["polygon_id"] >= 0]
    Image.MAX_IMAGE_PIXELS = None
    image = Image.open(outpath+'/'+sample+'_polygon.tif')
    draw = ImageDraw.Draw(image)


    spots_se = spots[spots["polygon_id"] >=0 ]
    square_centers = list(zip(spots_se['pxl_col_in_fullres'], spots_se['pxl_row_in_fullres']))

    square_size = 8

    for (cx, cy) in square_centers:
        top_left = (cx - square_size // 2, cy - square_size // 2)
        bottom_right = (cx + square_size // 2, cy + square_size // 2)
        draw.rectangle([top_left, bottom_right], outline=(173, 216, 230), width=1)  # Blue outline

    image.save(outpath+'/'+sample+'_spot_polygon.tif')



def aggregate_pro_sparse(outpath, h5path, sample):
    adata = sc.read_10x_h5(h5path)
    mtx = adata.X.tocsr()

    df = pd.read_csv(f"{outpath}{sample}_spot_to_polygon.csv")
    df = df.set_index("barcode")
    adata_barcodes = np.array(adata.obs_names)
    barcodes = np.array(df.index)
    barcode_to_idx = {b: i for i, b in enumerate(adata_barcodes)}
    idx = [barcode_to_idx[b] for b in barcodes]
    mtx = mtx[idx]
    polygon_ids = df["polygon_id"].values
    unique_polygons = np.unique(polygon_ids)
    poly_map = {p: i for i, p in enumerate(unique_polygons)}
    rows = [poly_map[p] for p in polygon_ids]
    cols = np.arange(len(polygon_ids))
    data = np.ones(len(polygon_ids))

    agg_mat = csr_matrix((data, (rows, cols)),
                         shape=(len(unique_polygons), len(polygon_ids)))
    mtx_sc = agg_mat @ mtx   # (polygon × genes)
    mtx_sc = mtx_sc.T
    mmwrite(f"{outpath}{sample}_mtx_sc.mtx", mtx_sc)
    pd.Series(adata.var_names).to_csv(f"{outpath}{sample}_genes.csv", index=False)
    pd.Series(unique_polygons).to_csv(f"{outpath}{sample}_polygons.csv", index=False)
    return mtx_sc, unique_polygons
  
  
  


