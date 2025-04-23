# apt-get install ffmpeg libsm6 libxext6  -y

import sys
sys.path.append('../../../../../SCS')
from scs import preprocessing, transformer, postprocessing

import sys
sys.path.append('../../../scripts/')

from paths import get_data_paths

DATA_FOLDER = get_data_paths('../../../../data_mapping.yml')['human_ovarian_cancer']
OUT_FOLDER = DATA_FOLDER / 'seg_method_results' / 'scs'
INPUT_FOLDER = OUT_FOLDER / 'input'
IMG_PATH = INPUT_FOLDER / 'dapi.tif'
BIN_PATH = INPUT_FOLDER / 'bin_df.tsv'

preprocessing.preprocess(
    BIN_PATH, IMG_PATH, False, None, 0, 0, patchsize=0, bin_size=3,
    n_neighbor=50, watershed_labels=None, out_path=OUT_FOLDER
)

transformer.train(0, 0, 0, 100, 0.0625, out_path=OUT_FOLDER)

cell_labels = postprocessing.postprocess(0, 0, 0, 3, 15, out_path=OUT_FOLDER, figsize=(10, 10))
