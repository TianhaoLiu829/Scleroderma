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





def combine_draw_expand_step1(distance_radius,image,sample):
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
    centroid_df.to_csv("/ix1/wchen/liutianhao/result/pathology_ST/spot_to_polygon"+sample+".csv",index=True)
    
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'_polygons_expanded_array.npy', polygons_expanded_array)
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'_expandions.npy',updated_expansions)
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'_overlap.npy',np.array(overlap))
    #Identify and plot squares
    spots = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/position_filter/"+sample+"_filtered_positions.csv")
    spots.index = spots.loc[:,'barcode']
    spots = spots.loc[:,['pxl_row_in_fullres', 'pxl_col_in_fullres']]
    from matplotlib.path import Path
    spot_to_polygon = -np.ones(len(spots), dtype=int)
    
    for i, poly in enumerate(polygons_expanded_array.transpose(0, 2, 1)):
        path = Path(poly)
        inside = path.contains_points(spots)
        spot_to_polygon[inside] = i  # this assigns the polygon ID to each spot inside

    spots['polygon_id'] = spot_to_polygon
    spots[spots["polygon_id"] >= 0].to_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/"+sample+"/"+sample+"_spot_to_polygon.csv", index=True)

    # Draw polygons on the image
    image_2 = Image.fromarray(image)
    img_with_polygons = draw_polygons_on_image(image_2, polygons_expanded, color=(0, 255, 127), width=4, fill=False)

    img_with_polygons.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'.tif', format="TIFF")
    
    # Draw squares on the image
    spots = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/"+sample+"/"+sample+"_spot_to_polygon.csv")
    Image.MAX_IMAGE_PIXELS = None
    image = Image.open('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'.tif')
    draw = ImageDraw.Draw(image)


    spots_se = spots[spots["polygon_id"] >=0 ]
    square_centers = list(zip(spots_se['pxl_col_in_fullres'], spots_se['pxl_row_in_fullres']))

    square_size = 8

    for (cx, cy) in square_centers:
        top_left = (cx - square_size // 2, cy - square_size // 2)
        bottom_right = (cx + square_size // 2, cy + square_size // 2)
        draw.rectangle([top_left, bottom_right], outline=(173, 216, 230), width=1)  # Blue outline

    image.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'_spot_polygon.tif')
    
    image = Image.open('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'.tif')
    draw = ImageDraw.Draw(image)

    spots_8 = pd.read_csv("/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/SpaceRangerOUT/M1/binned_outputs/square_008um/spatial/tissue_positions.csv")
    spots_8.index = spots_8.loc[:,'barcode']
    spots_8 = spots_8.loc[:,['pxl_row_in_fullres', 'pxl_col_in_fullres']]
    square_centers = list(zip(spots_8['pxl_col_in_fullres'], spots_8['pxl_row_in_fullres']))

    square_size = 32

    for (cx, cy) in square_centers:
        top_left = (cx - square_size // 2, cy - square_size // 2)
        bottom_right = (cx + square_size // 2, cy + square_size // 2)
        draw.rectangle([top_left, bottom_right], outline=(173, 216, 230), width=2)  # Blue outline

    image.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+'_spot_polygon_8.tif')



def combine_draw_expand_step1(distance_radius,image,sample,csv,index=''):

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
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'_polygons_nuclear.npy', polygons_c)
    polygons_expanded = expand_polygons(polygons_c,expansions)
    polygons_expanded_array = expand_polygons_array(polygons_c,expansions)
    polygons_expanded_array = np.stack(polygons_expanded_array)
    
    #restrict the expansion of polygons
    updated_polygon, updated_expansions, overlap = reduce_overlaps(polygons_expanded,polygons_c,expansions,step=0.5, min_expansion=2.0, max_iter=20)
    
    polygons_expanded_array = sample_polygon_boundaries(updated_polygon)  # shape (N, 2, 32)
    
    centroid_coords = np.array([(poly.centroid.x, poly.centroid.y) for poly in polygons_expanded])
    centroid_df = pd.DataFrame(centroid_coords, columns=['x', 'y'])
    centroid_df.to_csv("/ix1/wchen/liutianhao/result/pathology_ST/spot_to_polygon"+sample+".csv",index=True)
    
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'_polygons_expanded_array.npy', polygons_expanded_array)
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'_expandions.npy',updated_expansions)
    np.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'_overlap.npy',np.array(overlap))
    #Identify and plot squares
    spots = csv
    spots.index = spots.loc[:,'barcode']
    spots = spots.loc[:,['pxl_row_in_fullres', 'pxl_col_in_fullres']]
    from matplotlib.path import Path
    spot_to_polygon = -np.ones(len(spots), dtype=int)
    
    for i, poly in enumerate(polygons_expanded_array.transpose(0, 2, 1)):
        path = Path(poly)
        inside = path.contains_points(spots)
        spot_to_polygon[inside] = i  # this assigns the polygon ID to each spot inside

    spots['polygon_id'] = spot_to_polygon
    spots[spots["polygon_id"] >= 0].to_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/"+sample+"/"+sample+index+"_spot_to_polygon.csv", index=True)

    # Draw polygons on the image
    image_2 = Image.fromarray(image)
    img_with_polygons = draw_polygons_on_image(image_2, polygons_expanded, color=(0, 255, 127), width=4, fill=False)

    img_with_polygons.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'.tif', format="TIFF")
    
    # Draw squares on the image
    spots = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/"+sample+"/"+sample+index+"_spot_to_polygon.csv")
    Image.MAX_IMAGE_PIXELS = None
    image = Image.open('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'.tif')
    draw = ImageDraw.Draw(image)


    spots_se = spots[spots["polygon_id"] >=0 ]
    square_centers = list(zip(spots_se['pxl_col_in_fullres'], spots_se['pxl_row_in_fullres']))

    square_size = 8

    for (cx, cy) in square_centers:
        top_left = (cx - square_size // 2, cy - square_size // 2)
        bottom_right = (cx + square_size // 2, cy + square_size // 2)
        draw.rectangle([top_left, bottom_right], outline=(173, 216, 230), width=4)  # Blue outline

    image.save('/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/'+sample+'/'+sample+index+'_spot_polygon.tif')




#skin samples
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
    image =io.imread(sample2path[sample])
    csv = pd.read_csv(f"/ix1/wchen/liutianhao/result/pathology_ST/position_filter/{sample}_filtered_positions.csv")
    combine_draw_expand_step1(40,image,f'{sample}_2',csv)

-asd-asd-asd


#neo infection asthma
neo_coord = pd.read_csv('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/9066_aA7/square_002um/spatial/tissue_positions.csv')
image_neo =io.imread('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/Visium_Images/9066_aA7.tif')

spots = np.array([
    [min(neo_coord.loc[:,'pxl_col_in_fullres']), max(neo_coord.loc[:,'pxl_row_in_fullres'])],
    [min(neo_coord.loc[:,'pxl_col_in_fullres']), min(neo_coord.loc[:,'pxl_row_in_fullres'])],
    [max(neo_coord.loc[:,'pxl_col_in_fullres']), max(neo_coord.loc[:,'pxl_row_in_fullres'])],
    [max(neo_coord.loc[:,'pxl_col_in_fullres']), min(neo_coord.loc[:,'pxl_row_in_fullres'])]
])
# Step 1: Get bounding box from the spots
x_coords = spots[:, 0]
y_coords = spots[:, 1]

x_min = int(np.floor(x_coords.min()))
x_max = int(np.ceil(x_coords.max()))
y_min = int(np.floor(y_coords.min()))
y_max = int(np.ceil(y_coords.max()))

inte_x = (x_max-x_min)/10
inte_y = (y_max-y_min)/10

cropped_image_neo = image_neo[int(y_min):int(y_max),int(x_min):int(x_max), :]
cropped_image_neo = Image.fromarray(cropped_image_neo)

# Step 2: Crop the image using bounding box
# Remember: image is indexed as [row (y), col (x)]
for x in range(0, 5):
  for y in range(0, 5):
    cropped_image_neo = image_neo[int(y_min+inte_y*y*2):int(y_min+inte_y*(y+1)*2),int(x_min+inte_x*x*2):int(x_min+inte_x*(x+1)*2), :]
    pos_neo = pd.read_csv("/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/9066_aA7/square_002um/spatial/tissue_positions_filtered.csv")
    cropped_pos_neo = pos_neo[(pos_neo['pxl_row_in_fullres']<=int(y_min+inte_y*(y+1)*2))&(pos_neo['pxl_row_in_fullres']>=int(y_min+inte_y*y*2))&(pos_neo['pxl_col_in_fullres']<=int(x_min+inte_x*(x+1)*2))&(pos_neo['pxl_col_in_fullres']>=int(x_min+inte_x*x*2))]
    cropped_pos_neo.loc[:,'pxl_row_in_fullres'] = cropped_pos_neo['pxl_row_in_fullres']-int(y_min+inte_y*y*2)
    cropped_pos_neo.loc[:,'pxl_col_in_fullres'] = cropped_pos_neo['pxl_col_in_fullres']-int(x_min+inte_x*x*2)
    combine_draw_expand_step1(40,cropped_image_neo,'lung_adult',cropped_pos_neo,index=(str(x)+str(y)))

#find centers
def find_center(polygon_coords,path):
    polygon_coords_reshape = np.transpose(polygon_coords, (0, 2, 1))
    centroids = polygon_coords_reshape.mean(axis=1)
    centroids = pd.DataFrame(centroids, columns=['x', 'y'])
    centroids.to_csv(path, index=False)
    return centroids


centers = []
path = '/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_adult/lung_adult'
for i in range(5):
    for j in range(5):
        polygongs = np.load(f"{path}{i}{j}_polygons_expanded_array.npy")
        polygongs[:,0,:]+=int(y_min+inte_y*j*2)-y_min
        polygongs[:,1,:]+=int(x_min+inte_x*i*2)-x_min
        np.save(f"{path}{i}{j}_polygons_expanded_array_adj.npy", polygongs)
        centers_df = find_center( polygongs ,f"{path}{i}{j}_center.csv")
        centers_df.index = [f"{index}_{i}{j}" for index in centers_df.index]
        centers.append(centers_df)


centers_all = pd.concat(centers)
centers_all = pd.concat(centers)
centers_all.to_csv("/ix1/wchen/liutianhao/result/pathology_ST/single_cell_label/lung_adult/center_lung_adult.csv",index=True)

-asd-asd-asd
#lung asthma
adult_coord = pd.read_csv('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/9066_aA7/square_002um/spatial/tissue_positions.csv')
image_adult =io.imread('/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/Visium_Images/9066_aA7.tif')

spots = np.array([
    [min(adult_coord.loc[:,'pxl_col_in_fullres']), max(adult_coord.loc[:,'pxl_row_in_fullres'])],
    [min(adult_coord.loc[:,'pxl_col_in_fullres']), min(adult_coord.loc[:,'pxl_row_in_fullres'])],
    [max(adult_coord.loc[:,'pxl_col_in_fullres']), max(adult_coord.loc[:,'pxl_row_in_fullres'])],
    [max(adult_coord.loc[:,'pxl_col_in_fullres']), min(adult_coord.loc[:,'pxl_row_in_fullres'])]
])
# Step 1: Get bounding box from the spots
x_coords = spots[:, 0]
y_coords = spots[:, 1]

x_min = int(np.floor(x_coords.min()))
x_max = int(np.ceil(x_coords.max()))
y_min = int(np.floor(y_coords.min()))
y_max = int(np.ceil(y_coords.max()))

inte_x = (x_max-x_min)/10
inte_y = (y_max-y_min)/10

# Step 2: Crop the image using bounding box
# Remember: image is indexed as [row (y), col (x)]
for x in range(0, 5):
  for y in range(0, 5):
    cropped_image_adult = image_adult[int(y_min+inte_y*y*2):int(y_min+inte_y*(y+1)*2),int(x_min+inte_x*x*2):int(x_min+inte_x*(x+1)*2), :]
    pos_adult = pd.read_csv("/ix1/wchen/liutianhao/Eddens/Eddens.T_Visium_HD_09_24_2024/SpaceRanger_OUT/9066_aA7/square_002um/spatial/tissue_positions_filtered.csv")
    cropped_pos_adult = pos_adult[(pos_adult['pxl_row_in_fullres']<=int(y_min+inte_y*(y+1)*2))&(pos_adult['pxl_row_in_fullres']>=int(y_min+inte_y*y*2))&(pos_adult['pxl_col_in_fullres']<=int(x_min+inte_x*(x+1)*2))&(pos_adult['pxl_col_in_fullres']>=int(x_min+inte_x*x*2))]
    cropped_pos_adult.loc[:,'pxl_row_in_fullres'] = cropped_pos_adult['pxl_row_in_fullres']-int(y_min+inte_y*y*2)
    cropped_pos_adult.loc[:,'pxl_col_in_fullres'] = cropped_pos_adult['pxl_col_in_fullres']-int(x_min+inte_x*x*2)
    combine_draw_expand_step1(40,cropped_image_adult,'lung_adult',cropped_pos_adult,index=(str(x)+str(y)))





-a-dcv-adsc

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/I1_40x.tif')
csv = pd.read_csv("/ix1/wchen/liutianhao/result/pathology_ST/position_filter/I_1_filtered_positions.csv")
combine_draw_expand_step1(40,image,f'{sample}_2',csv)

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/J1_40x.tif')
combine_draw_expand_step1(40,image,'J_1')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/H1_40x.tif')
combine_draw_expand_step1(40,image,'H_1')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_02_27_2025/Visium_Images/H2_40x.tif')
combine_draw_expand_step1(40,image,'H_2')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8594_G1.tif')
combine_draw_expand_step1(40,image,'G_1')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/4945_G2.tif')
combine_draw_expand_step1(40,image,'G_2')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8691_A2.tif')
combine_draw_expand_step1(40,image,'A_2')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/E1.tif')
combine_draw_expand_step1(40,image,'E_1')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_01_28_2025/Visium_Image/8693_E2.tif')
combine_draw_expand_step1(40,image,'E_2')

image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/C.Torok_VisiumHD_08_2024/Visium_image/F1.tif')
combine_draw_expand_step1(40,image,'F_1')


#new_batch
image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/M1_1045.tif')
combine_draw_expand_step1(40,image,'M_1')
image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/K1_1323.tif')
combine_draw_expand_step1(40,image,'K_1')
image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/K2_1113.tif')
combine_draw_expand_step1(40,image,'K_2')
image =io.imread('/ix1/ktorok/shared/ktorok_rfe4_wam30/Torok.C_Visium_HD_05_14_2025/Visium Images/L1_1017.tif')
combine_draw_expand_step1(40,image,'L_1')


#CRC chongyue's method paper
crc_coord = pd.read_csv('/ix1/wchen/liutianhao/data/public_lung_spatial/CRC_10X/HD/P1_CRC/binned_outputs/square_002um/spatial/tissue_positions.csv')
image_crc =io.imread('/ix1/wchen/liutianhao/data/public_lung_spatial/CRC_10X/HD/P1_CRC/Visium_HD_Human_Colon_Cancer_P1_tissue_image.btf',plugin='tifffile')

spots = np.array([
    [min(crc_coord.loc[:,'pxl_col_in_fullres']), max(crc_coord.loc[:,'pxl_row_in_fullres'])],
    [min(crc_coord.loc[:,'pxl_col_in_fullres']), min(crc_coord.loc[:,'pxl_row_in_fullres'])],
    [max(crc_coord.loc[:,'pxl_col_in_fullres']), max(crc_coord.loc[:,'pxl_row_in_fullres'])],
    [max(crc_coord.loc[:,'pxl_col_in_fullres']), min(crc_coord.loc[:,'pxl_row_in_fullres'])]
])
# Step 1: Get bounding box from the spots
x_coords = spots[:, 0]
y_coords = spots[:, 1]

x_min = max(int(np.floor(x_coords.min())),0)
x_max = int(np.ceil(x_coords.max()))
y_min = int(np.floor(y_coords.min()))
y_max = int(np.ceil(y_coords.max()))

inte_x = (x_max-x_min)/10
inte_y = (y_max-y_min)/10

# Step 2: Crop the image using bounding box
# Remember: image is indexed as [row (y), col (x)]
for x in range(0, 5):
  for y in range(0, 5):
    cropped_image_crc = image_crc[int(y_min+inte_y*y*2):int(y_min+inte_y*(y+1)*2),int(x_min+inte_x*x*2):int(x_min+inte_x*(x+1)*2), :]
    pos_crc = pd.read_csv("/ix1/wchen/liutianhao/data/public_lung_spatial/CRC_10X/HD/P1_CRC/binned_outputs/square_002um/spatial/tissue_positions_filtered.csv")
    cropped_pos_crc = pos_crc[(pos_crc['pxl_row_in_fullres']<=int(y_min+inte_y*(y+1)*2))&(pos_crc['pxl_row_in_fullres']>=int(y_min+inte_y*y*2))&(pos_crc['pxl_col_in_fullres']<=int(x_min+inte_x*(x+1)*2))&(pos_crc['pxl_col_in_fullres']>=int(x_min+inte_x*x*2))]
    cropped_pos_crc.loc[:,'pxl_row_in_fullres'] = cropped_pos_crc['pxl_row_in_fullres']-int(y_min+inte_y*y*2)
    cropped_pos_crc.loc[:,'pxl_col_in_fullres'] = cropped_pos_crc['pxl_col_in_fullres']-int(x_min+inte_x*x*2)
    combine_draw_expand_step1(40,cropped_image_crc,'crc',cropped_pos_crc,index=(str(x)+str(y)))

-asd-asd-asd










