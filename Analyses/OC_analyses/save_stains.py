import tifffile
import skimage.io as sio

import sys
sys.path.append('../../scripts')
from paths import get_data_paths

dapi = tifffile.imread(get_data_paths()['human_ovarian_cancer'] / "morphology_focus/morphology_focus_0000.ome.tif", is_ome=False, level=1)

membrane = tifffile.imread(get_data_paths()['human_ovarian_cancer'] / "morphology_focus/morphology_focus_0001.ome.tif", is_ome=False, level=1)

sio.imsave(get_data_paths()['human_ovarian_cancer'] / 'processed' / 'dapi_lowres.tif', dapi)
sio.imsave(get_data_paths()['human_ovarian_cancer'] / 'processed' / 'membrane_lowres.tif', membrane)