from typing import Optional
import numpy as np
import pandas as pd
import skimage.io as sio
from skimage import draw


def process_clustermap_tile(
        out, tile_num, min_spot_per_cell: int = 5,
        gene_list: list = None, num_dims: int = 2, xy_radius: int = 0
    ):
    from ClusterMap.clustermap import ClusterMap

    spots_tile=out.loc[tile_num,'spots']
    dapi_tile=out.loc[tile_num,'img']

    ### instantiate model
    model_tile = ClusterMap(
        spots=spots_tile, dapi=dapi_tile, gene_list=gene_list, num_dims=num_dims,
        xy_radius=xy_radius,z_radius=0,fast_preprocess=False
    )

    ###preprocessing
    model_tile.preprocess(dapi_grid_interval=3, pct_filter=0)

    # ### segmentation
    model_tile.min_spot_per_cell = min_spot_per_cell
    model_tile.segmentation(cell_num_threshold=0.1,dapi_grid_interval=3,add_dapi=True,use_genedis=True)

    # ### plot cell segmentation results in spots (colored by cells)
    # model_tile.plot_segmentation(figsize=(4,4),s=0.01,plot_with_dapi=True,plot_dapi=True, show=False)
    return model_tile


def create_xenium_dapi_mask(
    bundle_path: str, shape: Optional[tuple] = None, scale: float = 1.0, return_mask: bool = True,
    labels: Optional[dict[str, int]] = None
):
    if shape is None:
        dapi = sio.imread(f'{bundle_path}/morphology_focus/morphology_focus_0000.ome.tif')[:,:,0]
        shape = dapi.shape

    nuclei_boundaries = pd.read_parquet(f'{bundle_path}/nucleus_boundaries.parquet')
    nuclei_polygons = (
        nuclei_boundaries.
        groupby('cell_id').
        apply(lambda x: [x[['vertex_x', 'vertex_y']].values / scale]).
        map(lambda x: x[0])
    )

    cell_labels = np.zeros(shape, dtype=np.uint32)

    for i, cell_id in enumerate(sorted(nuclei_polygons.keys())):
        vertices = nuclei_polygons[cell_id]
        lab = labels.get(cell_id, 0) if labels is not None else i + 1

        rr, cc = draw.polygon(vertices[:, 1], vertices[:, 0], shape)
        cell_labels[rr, cc] = lab

    if return_mask:
        return (cell_labels > 0)

    return cell_labels


def bin_2d_data(x_vals: np.ndarray, y_vals: np.ndarray, step=1):
    """
    Bin 2D data from a DataFrame into a grid with specified step size.

    Parameters:
    -----------
    x_vals : np.ndarray
        Input array containing x coordinates
    y_vals : np.ndarray
        Input array containing y coordinates
    step : float
        Size of the bins (default: 1)

    Returns:
    --------
    tuple
        (bin_counts, x_edges, y_edges, x_centers, y_centers)
        - bin_counts: 2D array with counts in each bin
        - x_edges: array of bin edges in x direction
        - y_edges: array of bin edges in y direction
        - x_centers: array of bin centers in x direction
        - y_centers: array of bin centers in y direction
    """
    # Calculate bin edges
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()

    # Add small offset to max to include the maximum values
    x_bins = np.arange(x_min, x_max + step, step)
    y_bins = np.arange(y_min, y_max + step, step)

    # Calculate histogram
    hist, x_edges, y_edges = np.histogram2d(x_vals, y_vals, bins=[x_bins, y_bins])

    # Calculate bin centers
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    x_bin_ids = np.digitize(x_vals, x_edges) - 1
    y_bin_ids = np.digitize(y_vals, y_edges) - 1

    return hist.T, x_edges, y_edges, x_centers, y_centers, x_bin_ids, y_bin_ids
