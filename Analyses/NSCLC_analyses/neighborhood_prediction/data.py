import os
import numpy as np
import pandas as pd
import torch
from sklearn.neighbors import NearestNeighbors

class SimpleVectorDataset(torch.utils.data.Dataset):

    def __init__(self, inputs, targs, ids):
        self.inputs = inputs
        self.targs = targs
        self.ids = ids
    
    def __len__(self):
        return self.inputs.shape[0]
    
    def __getitem__(self, idx):
        
        cur_inputs = torch.FloatTensor(self.inputs[idx, :])
        cur_targs = torch.FloatTensor(self.targs[idx, :])
        cur_ids = self.ids[idx]
        
        return cur_inputs, cur_targs, cur_ids

# TODO: Decide whether we should use z in defining nearest neighbors.

def build_targets(cell_coords, cell_type, num_nb):
    num_cells = cell_coords.shape[0]
    # generate integer labels:
    cell_type_list = np.sort(np.unique(cell_type))
    cell_type_to_index = {cell_type_list[i]: i for i in range(len(cell_type_list))}
    cell_type_int = np.array([cell_type_to_index[x] for x in cell_type])
    # NN lookups:
    nn_obj = NearestNeighbors(n_neighbors=num_nb+1, n_jobs=32)
    nn_obj.fit(cell_coords)
    nb_inds = nn_obj.kneighbors(cell_coords, return_distance=False)
    nb_inds = nb_inds[:, 1:]
    nb_classes = cell_type_int[nb_inds]
    nb_counts = np.zeros((num_cells, len(np.unique(cell_type_int))))
    for i in range(num_cells): # TODO: Faster implementation.
        for j in nb_classes[i, :]:
            nb_counts[i, j] += 1
    return nb_counts, cell_type_list

def parse_data(params):
    if params['dataset_name'] == 'gut_v1':
        # load raw data:
        df_seg = pd.read_csv(os.path.join(params['load_base'], 'segmentation/segmentation.csv'))
        df_cell_stats = pd.read_csv(os.path.join(params['load_base'], 'segmentation/segmentation_cell_stats.csv'))
        df_cell_assignments = pd.read_csv(os.path.join(params['load_base'], 'clustering/cell_assignment.csv'))
        df_cell_stats['leiden_final'] = df_cell_assignments['leiden_final'] # assuming the rows line up
        # clean up:
        cells_drop = df_cell_stats[df_cell_stats['leiden_final'] == 'Removed']['cell'].to_numpy()
        cells_drop = list(cells_drop)
        cells_drop.append(0) # Skip the "background" cell.
        cells_drop = np.array(cells_drop)
        df_seg = df_seg[~df_seg['cell'].isin(cells_drop)]
        df_cell_stats = df_cell_stats[df_cell_stats['leiden_final'] != 'Removed']
        # merge some classes:
        for i in range(len(params['class_merges'])):
            df_cell_stats.loc[df_cell_stats['leiden_final'].isin(params['class_merges'][i][0]), 'leiden_final'] = params['class_merges'][i][1]

        gene_list = sorted(np.unique(df_seg['gene'].to_numpy()))
        gene_name_to_gene_id = {gene_list[i]: i for i in range(len(gene_list))}

        num_cells = df_cell_stats.shape[0]
        cell_list = np.sort(df_cell_stats['cell'].to_numpy()) # defines order of outputs
        cell_id_to_index = {cell_list[i]: i for i in range(len(cell_list))}

        # format neighbor cell type data:
        cell_coords = np.stack(np.array([df_cell_stats['x'], df_cell_stats['y']]), axis=0).T
        cell_type = df_cell_stats['leiden_final'].to_numpy()
        nb_counts, cell_type_list = build_targets(cell_coords, cell_type, params['num_nb'])

        nb_counts_cell_ids = df_cell_stats['cell'].to_numpy()
        nb_counts_cell_id_to_index = {nb_counts_cell_ids[i]: i for i in range(len(nb_counts_cell_ids))}
        idx_nb = np.array([nb_counts_cell_id_to_index[cell_list[i]] for i in range(len(cell_list))])
        nb_counts = nb_counts[idx_nb, :]
        cell_type = cell_type[idx_nb]

        # format expression data:
        X = np.zeros((len(cell_list), len(gene_list)))
        for cell_id, cur_df in df_seg[['cell', 'gene']].groupby('cell'):
            cell_index = cell_id_to_index[cell_id]
            vc = cur_df['gene'].value_counts().reset_index().to_numpy()
            indices = np.array([gene_name_to_gene_id[vc[i, 0]] for i in range(vc.shape[0])])
            X[cell_index, indices] = vc[:, 1]
        
    elif params['dataset_name'] == 'nsclc':
        '''
        Data prep here is somewhat slow. Consider doing offline. 
        '''
        # load coordinates:
        df_meta = pd.read_csv(os.path.join(params['load_base'], 'Lung5_Rep1_metadata_file.csv'))
        df_meta['cell_ID_unique'] = df_meta['fov'].astype(str) + '_' + df_meta['cell_ID'].astype(str)
        assert len(np.unique(df_meta['cell_ID_unique'])) == df_meta.shape[0]
        # load cell type labels:
        df_cell_types = pd.read_csv(os.path.join(params['load_base'], 'annotation_adj.csv'))
        df_cell_types['cell_ID_unique'] = df_cell_types['fov'].astype(str) + '_' + df_cell_types['cell_id'].astype(str)
        # merge some classes:
        for i in range(len(params['class_merges'])):
            df_cell_types.loc[df_cell_types['cell_type'].isin(params['class_merges'][i][0]), 'cell_type'] = params['class_merges'][i][1]
        assert len(np.unique(df_cell_types['cell_ID_unique'])) == df_cell_types.shape[0]
        # merge dataframes:
        df = df_meta.merge(df_cell_types, on='cell_ID_unique', how='inner')

        cell_list = np.sort(df['cell_ID_unique'].to_numpy()) # defines order of outputs
        cell_id_to_index = {cell_list[i]: i for i in range(len(cell_list))}

        # load expression data:
        df_expr = pd.read_csv(os.path.join(params['load_base'], 'Lung5_Rep1_tx_file.csv'))
        df_expr = df_expr.loc[df_expr['CellComp'] != '0', :] # drop background "cell"
        df_expr['cell_ID_unique'] = df_expr['fov'].astype(str) + '_' + df_expr['cell_ID'].astype(str)
        df_expr = df_expr.loc[df_expr['cell_ID_unique'].isin(cell_list), :] # drop cells without annotations

        gene_list = np.unique(df_expr['target'].to_numpy())
        gene_name_to_gene_id = {gene_list[i]: i for i in range(len(gene_list))}

        # format neighbor cell type data:
        cell_coords = df[['CenterX_global_px', 'CenterY_global_px']].to_numpy().astype(float) # Note: not reversing y-axis, which is required for plotting on image. 
        cell_type = df['cell_type'].to_numpy().astype(str)
        nb_counts, cell_type_list = build_targets(cell_coords, cell_type, params['num_nb'])
        nb_counts_cell_ids = df['cell_ID_unique'].to_numpy()
        nb_counts_cell_id_to_index = {nb_counts_cell_ids[i]: i for i in range(len(nb_counts_cell_ids))}
        idx_nb = np.array([nb_counts_cell_id_to_index[cell_list[i]] for i in range(len(cell_list))])
        nb_counts = nb_counts[idx_nb, :]
        cell_type = cell_type[idx_nb]

        # format expression data:
        X = np.zeros((len(cell_list), len(gene_list)))
        for cell_id, cur_df in df_expr[['cell_ID_unique', 'target']].groupby('cell_ID_unique'):
            cell_index = cell_id_to_index[cell_id]
            vc = cur_df['target'].value_counts().reset_index().to_numpy()
            indices = np.array([gene_name_to_gene_id[vc[i, 0]] for i in range(vc.shape[0])])
            X[cell_index, indices] = vc[:, 1]

    else:
        raise NotImplementedError
    
    # preprocess targets:
    if params['targ_proc'] == 'binarize':
        nb_counts = np.array(nb_counts > 0).astype(float)
    else:
        raise NotImplementedError
    
    # preprocess inputs:
    if params['input_proc'] == 'norm_scale_log1p':
        X /= np.sum(X, axis=1, keepdims=True)
        X *= 1e4
        X += 1
        X = np.log(X)
    elif params['input_proc'] == 'binarize':
        X = np.array(X > 0).astype(np.float)
    else:
        raise NotImplementedError
    
    return cell_coords, X, nb_counts, cell_list, cell_type_list, gene_list, cell_type

def get_loaders(params):
    coords, inputs, targs, ids, class_names, gene_list, focal_cell_type = parse_data(params)
    if params['class_names_whitelist'] is not None:
        class_mask = np.array([class_names[i] in params['class_names_whitelist'] for i in range(len(class_names))])
        assert np.any(class_mask) # ensure at least one class is included
        targs = targs[:, class_mask]
        class_names = class_names[class_mask]
    if params['gene_names_blacklist'] is not None:
        gene_mask = np.array([gene_list[i] not in params['gene_names_blacklist'] for i in range(len(gene_list))])
        assert np.any(gene_mask) # ensure at least one gene is included
        inputs = inputs[:, gene_mask]
        gene_list = np.array(gene_list)[gene_mask].tolist()
    if params['input_cell_type'] is not None:
        focal_mask = focal_cell_type == params['input_cell_type']
        coords = coords[focal_mask, :]
        inputs = inputs[focal_mask, :]
        targs = targs[focal_mask, :]
        ids = ids[focal_mask]
        focal_cell_type = focal_cell_type[focal_mask]
    in_dim = inputs.shape[1]
    out_dim = targs.shape[1]
    loaders = {split: build_dataset(coords, inputs, targs, ids, focal_cell_type, split, params) for split in ['train', 'val', 'test']}
    class_counts = loaders['train'].dataset.targs.sum(axis=0)
    return loaders, in_dim, out_dim, class_counts, class_names, gene_list

def build_dataset(coords, inputs, targs, ids, focal_cell_type, split, params):
    
    # generate splits:
    if params['split_type'] == 'horizontal':
        x_min, x_max = params['split_intervals'][split]
        x = coords[:, 0]
        is_valid = (x > np.quantile(x, x_min)) & (x < np.quantile(x, x_max))
    elif params['split_type'] == 'vertical':
        y_min, y_max = params['split_intervals'][split]
        y = coords[:, 1]
        is_valid = (y > np.quantile(y, y_min)) & (y < np.quantile(y, y_max))
    elif params['split_type'] == 'checkerboard':
        # generate a checkerboard ID for each data point:
        x = coords[:, 0]
        y = coords[:, 1]
        xmin = np.min(x)
        xmax = np.max(x)
        check_size_x = (xmax - xmin) / params['checkerboard_params']['num_x']
        ymin = np.min(y)
        ymax = np.max(y)
        check_size_y = (ymax - ymin) / params['checkerboard_params']['num_y']
        x_assign = np.floor((x-xmin) / check_size_x).astype(int)
        y_assign = np.floor((y-ymin) / check_size_y).astype(int)
        check_ids = np.array([str(x_assign[i]) + '_' + str(y_assign[i]) for i in range(len(x_assign))])
        unique_check_ids = np.unique(check_ids)
        rng = np.random.default_rng(seed=6372)
        idx_rand = rng.permutation(len(unique_check_ids))
        unique_check_ids_shuffled = unique_check_ids[idx_rand]
        num_train = np.round(params['checkerboard_params']['train'] * len(unique_check_ids)).astype(int)
        num_val = np.round(params['checkerboard_params']['val'] * len(unique_check_ids)).astype(int)
        check_ids_train = unique_check_ids_shuffled[:num_train]
        check_ids_val = unique_check_ids_shuffled[num_train:(num_train+num_val)]
        check_ids_test = unique_check_ids_shuffled[(num_train+num_val):]
        is_valid_train = np.array([cur_check_id in check_ids_train for cur_check_id in check_ids])
        if split == 'train':
            is_valid = is_valid_train
        elif split == 'val':
            is_valid = np.array([cur_check_id in check_ids_val for cur_check_id in check_ids])
        elif split == 'test':
            is_valid = np.array([cur_check_id in check_ids_test for cur_check_id in check_ids])
        # Note: Inefficient; need to compute all three splits every time we ask for one. 
    else:
        raise NotImplementedError
    if params['standardize_inputs']:
        train_means = np.mean(inputs[is_valid_train, :], axis=0, keepdims=True)
        train_stds = np.std(inputs[is_valid_train, :], axis=0, keepdims=True)
        train_stds[train_stds == 0] = 1.0
    inputs = inputs[is_valid, :]
    targs = targs[is_valid, :]
    ids = ids[is_valid]
    if params['standardize_inputs']:
        inputs -= train_means
        inputs /= train_stds


    # randomization baseline: 
    if params['randomize'] and (split == 'train'):
        print('randomizing assignment of transcriptomes to cells (within types) in training set')
        focal_cell_type_masked = focal_cell_type[is_valid] # subset list of focal cell types to current split
        rng = np.random.default_rng(seed=9753)
        for cur_cell_type in sorted(np.unique(focal_cell_type_masked)):
            cur_idx = np.ravel(np.argwhere(cur_cell_type == focal_cell_type_masked))
            if len(cur_idx) == 0:
                continue
            cur_idx_shuffled = cur_idx[rng.permutation(len(cur_idx))]
            inputs[cur_idx, :] = inputs[cur_idx_shuffled, :]
    
    # define dataset:
    ds = SimpleVectorDataset(inputs, targs, ids)
    loader = torch.utils.data.DataLoader(
        ds, 
        batch_size=params['batch_size'], 
        shuffle=split=='train', 
        num_workers=params['num_workers'], 
        pin_memory=True
    )
    
    return loader

