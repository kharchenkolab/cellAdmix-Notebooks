import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import lightning as L
from lightning.pytorch.loggers import TensorBoardLogger
from sklearn.metrics import average_precision_score

import data, models

import sys
sys.path.append('../../../scripts')
from paths import get_data_paths

class NeighborSupervisedModule(L.LightningModule):

    def __init__(self, params):
        super().__init__()
        self.save_hyperparameters()
        self.params = params
        self.model = models.get_model(self.params['model_name'], params)
        self.loss_module = torch.nn.BCEWithLogitsLoss(pos_weight=torch.FloatTensor(params['class_weights']))
        self.example_input_array = torch.zeros((1, self.params['in_dim']), dtype=torch.float32)
        self.train_step_preds = []
        self.train_step_targs = []
        self.validation_step_preds = []
        self.validation_step_targs = []

    def forward(self, inputs):
        return self.model(inputs)

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.params['base_lr'])
        return optimizer

    def training_step(self, batch, batch_idx):
        inputs, targs, ids = batch
        preds = self.model(inputs)
        loss = self.loss_module(preds, targs)
        # add L1 penalty?
        self.log('train_loss', loss, batch_size=len(ids))
        self.train_step_preds.append(preds.detach())
        self.train_step_targs.append(targs)
        return loss

    def validation_step(self, batch, batch_idx):
        inputs, targs, ids = batch
        preds = self.model(inputs)
        loss = self.loss_module(preds, targs)
        self.log('val_loss', loss, batch_size=len(ids))
        self.validation_step_preds.append(preds)
        self.validation_step_targs.append(targs)

    def on_train_epoch_end(self):
        # compile predictions:
        all_preds = np.array(torch.squeeze(torch.concatenate(self.train_step_preds, dim=0)).cpu())
        all_targs = np.array(torch.squeeze(torch.concatenate(self.train_step_targs, dim=0)).cpu())
        if all_preds.ndim == 1:
            all_preds = all_preds[:, np.newaxis]
            all_targs = all_targs[:, np.newaxis]
        # compute metrics:
        ap = [average_precision_score(all_targs[:, i], all_preds[:, i]) for i in range(all_preds.shape[1])]
        self.log('train_map', np.mean(ap))
        for i in range(len(ap)):
            self.log(f'train_ap/{self.params["class_names"][i]}', ap[i])
        # clear memory:
        self.train_step_preds.clear()
        self.train_step_targs.clear()

    def on_validation_epoch_end(self):
        # compile predictions:
        all_preds = np.array(torch.squeeze(torch.concatenate(self.validation_step_preds, dim=0)).cpu())
        all_targs = np.array(torch.squeeze(torch.concatenate(self.validation_step_targs, dim=0)).cpu())
        if all_preds.ndim == 1:
            all_preds = all_preds[:, np.newaxis]
            all_targs = all_targs[:, np.newaxis]
        # compute metrics:
        ap = np.array([average_precision_score(all_targs[:, i], all_preds[:, i]) for i in range(all_preds.shape[1])])
        ap_chance = np.array([np.mean(all_targs[:, i]) for i in range(all_preds.shape[1])])
        self.log('val_map', np.mean(ap))
        self.log('val_map_over_chance', np.mean(ap - ap_chance))
        for i in range(len(ap)):
            self.log(f'val_ap/{self.params["class_names"][i]}', ap[i])
            self.log(f'val_ap_over_chance/{self.params["class_names"][i]}', ap[i] - ap_chance[i])

        # clear memory:
        self.validation_step_preds.clear()
        self.validation_step_targs.clear()

    def test_step(self, batch, batch_idx):
        inputs, targs, ids = batch
        preds = self.model(inputs)
        loss = self.loss_module(preds, targs)
        self.log('test_loss', loss, batch_size=len(ids))

    def predict_step(self, batch, batch_idx):
        inputs, targs, ids = batch
        preds = self.model(inputs)
        return preds, targs, ids

if __name__ == '__main__':
    data_mapping = get_data_paths('../../../data_mapping.yml')

    torch.set_float32_matmul_precision('high') # high, medium
    # Log train / val / test counts?
    params = {}

    # data params:
    params['run_name'] = '2024_04_11_gut_in_enterocyte_out_goblet'
    params['dataset_name'] = 'gut_v1' # gut_v1, nsclc
    params['randomize'] = False # whether to shuffle transcriptomes within each cell type
    if params['dataset_name'] == 'gut_v1':
        params['load_base'] = data_mapping['mouse_gut']
        # params['class_names_whitelist'] = ['Enterocyte']
        params['class_names_whitelist'] = ['Goblet']
        # params['class_names_whitelist'] = None # full multi-label
        # params['input_cell_type'] = None # which cell type to use as input; None -> use all
        params['input_cell_type'] = 'Enterocyte' # which cell type to use as input; None -> use all
        # params['input_cell_type'] = 'Goblet'
        params['class_merges'] = [(['Enterocyte (Bottom Villus)', 'Enterocyte (Mid Villus)', 'Enterocyte (Top Villus)'], 'Enterocyte')]
        params['gene_names_blacklist'] = None
    elif params['dataset_name'] == 'nsclc':
        params['load_base'] = data_mapping['nsclc_lung5_1']
        # params['class_names_whitelist'] = ['fibroblast']
        params['class_names_whitelist'] = ['tumor']
        params['class_merges'] = [(['tumor 5', 'tumor 6', 'tumor 9', 'tumor 12', 'tumor 13'], 'tumor')]
        # params['input_cell_type'] = None # which cell type to use as input; None -> use all
        # params['input_cell_type'] = 'fibroblast'
        # params['input_cell_type'] = 'tumor'
        params['input_cell_type'] = 'neutrophil'
        # params['input_cell_type'] = 'macrophage'
        params['gene_names_blacklist'] = None
    else:
        params['class_name_whitelist'] = None
    params['split_type'] = 'checkerboard' # horizontal, vertical, checkerboard
    params['num_nb'] = 6 # number of neighbors used to define prediction targets
    params['split_intervals'] = {
        'train': (0.0, 0.5),
        'val': (0.5, 0.75),
        'test': (0.75, 1.0)
    }
    params['checkerboard_params'] = {
        'num_x': 10,
        'num_y': 10,
        'train': 0.5,
        'val': 0.25,
        'test': 0.25
    }
    params['input_proc'] = 'norm_scale_log1p' # norm_scale_log1p, binarize
    params['targ_proc'] = 'binarize' # binarize
    params['standardize_inputs'] = True

    # model params:
    params['enc_depth'] = 0
    params['enc_width'] = 32
    params['p_drop'] = 0.0

    # optimizer params:
    params['base_lr'] = 1e-4
    params['batch_size'] = 64
    params['max_epochs'] = 300
    params['num_workers'] = 12

    # misc params:
    params['log_every_n_steps'] = 10

    # build dataset:
    loaders, in_dim, out_dim, class_counts, class_names, gene_list = data.get_loaders(params)
    if 0:
        # positive weight balance relative to other classes
        params['class_weights'] = 1.0 / class_counts
        params['class_weights'] /= np.sum(params['class_weights'])
    else:
        # positive weight balance each class for presence/absence
        neg_counts = len(loaders['train'].dataset) - class_counts
        params['class_weights'] = neg_counts / class_counts
    assert np.all(params['class_weights'] > 0)
    params['class_names'] = class_names
    params['gene_list'] = gene_list
    assert not np.any(np.isnan(params['class_weights']))

    # model params:
    params['model_name'] = 'multilabel_mlp'
    params['in_dim'] = in_dim
    params['out_dim'] = out_dim

    logger = TensorBoardLogger(save_dir='./', version=params['run_name'])

    supervised_model = NeighborSupervisedModule(params)
    lr_monitor = L.pytorch.callbacks.LearningRateMonitor(logging_interval='step')
    checkpoint_callback = L.pytorch.callbacks.ModelCheckpoint(
        save_top_k=3,
        monitor='val_map',
        mode='max',
        save_last=True, # Note: This refers to latest checkpoint, not final after training.
        save_on_train_epoch_end=False, # save after val, not after train
        )
    trainer = L.Trainer(
        devices=[1],
        accelerator='gpu',
        max_epochs=params['max_epochs'],
        log_every_n_steps=params['log_every_n_steps'],
        callbacks=[checkpoint_callback, lr_monitor],
        logger=logger,
        )

    trainer.fit(model=supervised_model, train_dataloaders=loaders['train'], val_dataloaders=loaders['val'])