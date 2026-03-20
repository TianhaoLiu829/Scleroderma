from skimage import io
from Cell_segmentation_utils import combine_draw_expand_step1, aggregate_pro_sparse

#your path
work_path='/ihome/wchen/tianhao/'

#run the pipeline
img =io.imread(f'{work_path}/Cell_segmentation_and_aggregation/data/M1_1045.tif')
combine_draw_expand_step1(img,'M_1',f'{work_path}/Cell_segmentation_and_aggregation/data/M_1_filtered_positions.csv',f'{work_path}/Cell_segmentation_and_aggregation/result/',40)
out_mtx = aggregate_pro_sparse(f'{work_path}/Cell_segmentation_and_aggregation/result/mtx/',f'{work_path}/Cell_segmentation_and_aggregation/data/filtered_feature_bc_matrix_M1.h5','M_1')
