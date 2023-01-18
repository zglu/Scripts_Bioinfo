# bin the gem file and save to h5ad

import warnings
warnings.filterwarnings('ignore')
import stereo as st

# load the bin1 matrix file and select bin
data_path = './E9.5_E2S4_GEM_bin1.tsv'
data = st.io.read_gem(file_path=data_path, bin_size=50)

# save to h5ad file 
# write into anndata for Seurat/scanpy
st.io.stereo_to_anndata(data, flavor='seurat', output='mouse_E9_5_bin50_seurat.h5ad')
#st.io.stereo_to_anndata(data, flavor='scanpy', output='mouse_E9_5_bin50_scanpy.h5ad')
