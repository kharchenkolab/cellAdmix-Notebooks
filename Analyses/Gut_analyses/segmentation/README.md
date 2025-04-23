# Segmentation of the MERFISH Mouse Ileum dataset

This folder contains the code for the following methods:

- [BIDCell](./bidcell_mouse_ileum.ipynb). Please, note that the method relies on [merfish_ileum_config.yaml](./merfish_ileum_config.yaml), which should be updated with your local paths.
- [ClusterMap](./clustermap_mouse_ileum.ipynb)
- [ComSeg](./comseg_mouse_ileum.ipynb). It also requires the [preprocessing notebook](./comseg_prepare_data.ipynb).
- [DeepCell](./deepcell_mouse_ileum.ipynb) requires registration or [web UI](https://www.deepcell.org/predict). Registration failed, so we used the web UI. The notebook exports the image and imports the results.
- [SCS](./scs_mouse_ileum.ipynb). By default the method works using Watershed segmentation, which has notoriously low quality. So, we used [a pathced version](https://github.com/VPetukhov/SCS), which allows passing an arbitrary segmentation and enables some utility configurations.

We also tested ProSeg using the following CLI command:

```bash
proseg -x x -y y -z z --gene-column gene \
    --cell-id-column prior_segmentation --cell-id-unassigned 0 \
    -t 60 --enforce-connectivity --qv-column qc_score --detect-layers \
    ../../../../data/mouse_gut/merfish/segmentation/segmentation.csv
```