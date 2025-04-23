from bidcell import BIDCellModel
import pandas as pd
import numpy as np
import skimage.io as sio

import torch # Otherwise BIDCell fails version check with conda torch
torch.__version__ = '.'.join(torch.__version__.split('.')[:3])

model = BIDCellModel("xenium_config.yaml")

print("Preprocess...")
model.preprocess()

print("Train...")
model.train()

print("Predict...")
model.predict()